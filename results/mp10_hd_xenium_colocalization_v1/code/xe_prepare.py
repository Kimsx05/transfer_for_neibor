from base import *
import gzip,hashlib,matplotlib
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from matplotlib.colors import to_hex
from datetime import date
source=P/'highres_ST/highres_spatial_audit_v2_external/00_local_input_inventory/GSE326226_h5'
refpath=next((P/'subtype/v2/reference_training_v2/04_gene_harmonization').glob('*h5ad'))
sets={'Chol_biosynthesis_8':'HMGCR MVD FDPS FDFT1 SQLE EBP DHCR7 DHCR24'.split(),
      'MP10_covered':(P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split()}
for p in sorted((O/'shared/msigdb').glob('*.json')):sets[p.stem]=json.loads(p.read_text())[p.stem]['geneSymbols']
common=None;coverage=[];inputrows=[]
for s in XE:
    d=O/'xenium'/s; d.mkdir(exist_ok=True)
    p=next(source.glob('*_'+s+'_*h5'));x,ids,g=tenx(p);keep=g.feature_type=='Gene Expression';x=x[:,keep.to_numpy()];g=g[keep].copy()
    assert g.name.is_unique and g.id.is_unique and pd.Index(ids).is_unique
    g['cells_detected']=x.getnnz(axis=0);g['transcripts']=np.asarray(x.sum(0)).ravel()
    put(g,d/'panel.tsv.gz'); common=set(g.name) if common is None else common&set(g.name)
    cf=next((O/'xenium/geo').glob('*_'+s+'_cells*'));bf=next((O/'xenium/geo').glob('*_'+s+'_cell_boundaries*'))
    cells=pd.read_csv(cf).set_index('cell_id');assert cells.index.is_unique and set(cells.index)==set(ids);cells=cells.loc[ids]
    with gzip.open(bf,'rb') as f:b=pd.read_parquet(f)
    vertex=b.groupby('cell_id').size();assert set(vertex.index)==set(ids);cells['boundary_vertices']=vertex
    cells['boundary_degenerate']=cells.boundary_vertices<4
    cells['area_outlier']=(cells.cell_area<=cells.cell_area.quantile(.005))|(cells.cell_area>=cells.cell_area.quantile(.995))
    cells['nFeature']=x.getnnz(axis=1);cells['nGeneTranscripts']=np.asarray(x.sum(1)).ravel()
    assert (cells.nGeneTranscripts==cells.transcript_counts).all() and np.isfinite(cells[['x_centroid','y_centroid']]).all().all()
    # Some boundaries approximate the original mask with a small number of vertices.
    cells.to_parquet(d/'cells.parquet')
    cells[["area_outlier","boundary_degenerate"]].to_csv(d/"flags.tsv.gz",sep="\t",index_label="cell_id")
    for key,genes in sets.items():
        for gene in genes:
            hit=g[g.name==gene];coverage.append((s,key,gene,not hit.empty,int(hit.cells_detected.iloc[0]) if len(hit) else None,int(hit.transcripts.iloc[0]) if len(hit) else None,len(ids)))
    inputrows.extend([(s,'expression',str(p)),(s,'cell centroids + QC',str(cf)),(s,'original cell polygons',str(bf))])
cols=['sample','program','gene','measured','cells_detected','transcripts','total_cells']
cov=pd.DataFrame(coverage,columns=cols);put(cov,O/'shared/panel_coverage.tsv')
summary=cov.groupby(['sample','program']).agg(set_size=('gene','size'),measured=('measured','sum'),detected=('cells_detected',lambda x:(x>0).sum()),cells=('total_cells','first')).reset_index()
put(summary,O/'shared/panel_summary.tsv')
pd.Series(sorted(common)).to_csv(O/'shared/common_panel.txt',index=False,header=False)
(O/'shared/programs.json').write_text(json.dumps({k:[g for g in v if g in common] for k,v in sets.items()},indent=2))
# Materialize only the balanced reference and the common panel for the transfer.
cells=pd.read_csv(P/'mp10_sc_c2l_final_annotation_v1/tables/analysis_cells.tsv.gz',sep='\t')
q=pd.read_csv(P/'mp10_cellchat_high_definition_sensitivity_v1/tables/specimen_quantiles_coverage.tsv',sep='\t').set_index('sample_uid')
from frozen_truth import frozen_truth
frozen=frozen_truth()
cells['H05_truth']=cells.barcode.map(frozen.H05_truth)
epi=cells.H05_truth.notna()
assert cells.loc[epi,'analysis_celltype_c2l_final'].eq('Epithelial').all() and epi.sum()==len(frozen)
arch=pd.read_csv(P/'mp10_cellchat_high_definition_sensitivity_v1/tables/membership/C30__ALL__H05.tsv.gz',sep='\t')
assert set(arch.loc[arch.candidate_pool=='H05','barcode'])<=set(cells.loc[cells.H05_truth.eq('H05').fillna(False),'barcode'])
rng=np.random.default_rng(123)
def balanced(df,strata,n):
    y=df.copy();y['random']=rng.random(len(y));y=y.sort_values('random');y['within']=y.groupby(strata).cumcount()
    return y.sort_values(['within','random']).head(n).index
ids1=np.concatenate([balanced(gr,'sample_uid',450 if label!='Epithelial' else 1200) for label,gr in cells.groupby('analysis_celltype_c2l_final')])
# Representative per-specimen validation includes both states at natural prevalence.
ids2=np.concatenate([rng.choice(gr.index,min(len(gr),240),False) for _,gr in cells[epi].groupby('sample_uid')])
# Reproduce the original balanced sampling frame at its serialized tail-count boundary.
# Labels used by the transfer itself are taken from explicit archived barcode membership.
sampling_frame=cells.loc[epi].copy()
sampling_frame['sampling_state']='nonH05'
for sample,group in sampling_frame.groupby('sample_uid'):
    n=int(q.loc[sample,'n_H05'])
    if n:
        hi=group.MP10_cell.nlargest(n).min()
        sampling_frame.loc[group.index[group.MP10_cell>=hi],'sampling_state']='H05'
ids3=np.concatenate([balanced(gr,'sample_uid',1800) for _,gr in sampling_frame.groupby('sampling_state')])
cells['identity_ref']=cells.index.isin(ids1);cells['validate']=cells.index.isin(ids2);cells['state_ref']=cells.index.isin(ids3)
meta=cells[cells.identity_ref|cells.validate|cells.state_ref].copy().set_index('barcode')
with h5py.File(refpath) as f:
    ids=f['obs/_index/values'].asstr()[:];genes=f['var/_index/values'].asstr()[:]
    jj=np.flatnonzero(np.isin(ids,meta.index));assert len(jj)==len(meta)
    sel=np.flatnonzero(np.isin(genes,list(common))); ptr=f['X/indptr'][:];parts=[]
    for a in range(0,len(ids),1024):
        rows=jj[(jj>=a)&(jj<a+1024)]-a
        if not len(rows):continue
        offsets=ptr[a:a+1025];lo,hi=offsets[[0,-1]]
        z=sparse.csr_matrix((f['X/data'][lo:hi],f['X/indices'][lo:hi],offsets-lo),shape=(len(offsets)-1,len(genes)))
        parts.append(z[rows][:,sel])
    x=sparse.vstack(parts,format='csr'); assert np.all(x.data==np.floor(x.data))
    d=O/'shared/reference';d.mkdir(exist_ok=True);export(x,genes[sel],ids[jj],d);put(meta.loc[ids[jj]].reset_index(),d/'membership.tsv.gz')
rules=runpy.run_path(str(A/'cell_annotation_v2/code/marker_rules.py'))
major=runpy.run_path(str(A/'cell_annotation_v2_1/code/programs.py'))['MAJOR']; sub=rules['SUBTYPE'].copy();anchors=rules['ANCHORS'].copy()
sub['C1QC TAM']+=sub['FOLR2 TAM'];anchors['C1QC TAM']+=anchors['FOLR2 TAM']
sub['Naive/Memory B']=sub['Naive B']+sub['Memory B']
anchors.update({'Classical Mono':['FCN1','S100A12','CCR2'],'LAMP3 DC':['LAMP3','FSCN1'],'cDC1':['CLEC9A','XCR1'],'cDC2':['CD1C','CLEC10A','FCER1A'],'Capillary EC':['CA4','RGCC'],'Lymphatic EC':['PROX1','LYVE1'],'rCAF':['CCL19','CCL21']})
m=cells.groupby('analysis_celltype_c2l_final').analysis_major_c2l_final.first().to_dict()
tax=pd.read_csv(P/'subtype/v2/reference_training_v2/11_final_spatial_reference/00_32_to_30_label_mapping.csv').rename(columns={'reference_celltype_v2':'original','spatial_reference_celltype_v2_final':'label'})
tax['major']=tax.label.map(m);tax['markers']=[' '.join(sub.get(k,major[m[k]])) for k in tax.label];tax['anchors']=[' '.join(anchors.get(k,sub.get(k,major[m[k]]))) for k in tax.label]
tax['status']='included in frozen reference; query support evaluated separately';put(tax,O/'shared/taxonomy.tsv')
(O/'shared/major_markers.json').write_text(json.dumps(major,indent=2))
palette=json.loads((A/'cell_annotation_v2_1/metadata/annotation_palette.json').read_text())
for i,k in enumerate(sorted(m)):palette.setdefault(k,to_hex(matplotlib.colormaps['turbo']((i*0.61803398875)%1)))
palette.update({'Mixed':'#ee7733','Mural':'#987450','Chol-high Epi':'#ac1645','other Epi':'#e4b286'})
(O/'shared/colors.json').write_text(json.dumps(palette,indent=2))
inputrows.append(('shared','reference raw counts',str(refpath)))
put(pd.DataFrame(inputrows,columns=['sample','role','path']),O/'shared/xenium_inputs.tsv')
print(summary.to_string(index=False),flush=True)

