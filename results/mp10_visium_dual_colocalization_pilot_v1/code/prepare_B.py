from shared import *
from sparse_io import strings,pack,slice_h5
from scipy import sparse,io
import h5py
ann=read(A/'analysis_cell_metadata.tsv.gz')
fine='analysis_celltype_c2l_final'; major='analysis_major_c2l_final'
mapping=ann[[fine,major]].drop_duplicates().sort_values([major,fine]);assert mapping[fine].is_unique and len(mapping)==30
mapping.columns=['subtype','major'];save(mapping,'annotation_mapping.tsv')
old=read(P/'reference_sampling_ids.tsv.gz'); epi=old[old.label.isin(['Epi_H05','Epi_nonH05'])].copy()
pool=read(R/'mp10_cellchat_high_definition_sensitivity_v1/tables/epi_full_scores_pools.tsv.gz')
audit=epi.merge(pool[['barcode','sample_uid','H05']],on=['barcode','sample_uid'],validate='one_to_one');assert len(audit)==len(epi) and ((audit.label=='Epi_H05')==audit.H05).all()
cand=ann[ann.sample_uid.isin(old.sample_uid.unique())&(ann[major]!='Epithelial')]
assert len(old.sample_uid.unique())==37 and not cand.sample_uid.str.contains('CNP0000460_P01|T4463209').any()
rng=np.random.default_rng(20261125)
tme=cand.loc[np.concatenate([rng.choice(g.index,min(len(g),35),replace=False) for _,g in cand.groupby(['sample_uid',fine])])].sort_values(['sample_uid',fine,'barcode']).copy()
tme['label']=tme[fine];epi['model']='Epi';tme['model']='TME'
ref=pd.concat([epi,tme],ignore_index=True)[['model','barcode','sample_uid','dataset','patient_uid_recorded',fine,major,'label']]
assert ref.barcode.is_unique and tme.label.nunique()==29
save(ref,'reference_manifest.tsv.gz')
counts=ref.groupby(['model','sample_uid','label']).size().rename('selected').reset_index()
orig=pd.concat([pool.assign(model='Epi',label=np.where(pool.H05,'Epi_H05','Epi_nonH05')),cand.assign(model='TME',label=cand[fine])])
orig=orig[orig.sample_uid.isin(old.sample_uid.unique())]
counts=counts.merge(orig.groupby(['model','sample_uid','label']).size().rename('eligible').reset_index(),validate='one_to_one')
counts['reference_share']=counts.selected/counts.groupby('model').selected.transform('sum');save(counts,'reference_sampling.tsv')
source=read(GENES);gs=[]
for name,key in PATHWAYS.items():
 g=source[source.gs_name==key].drop_duplicates('gene_symbol');assert len(g)>10
 for variant in ['full','minus20']:
  sub=g if variant=='full' else g[~g.gene_symbol.isin(TOP)]
  gs.append(sub.assign(module=name,variant=variant,Top20_overlap=sub.gene_symbol.isin(TOP)))
save(pd.concat(gs),'gene_sets.tsv')
niche=R/'neighborhood_ring_harmony_v3/03_main/SAMPLE_R12_LEIDEN/labels.parquet'
m=read(P/'input_manifest.tsv')
files=[]; queries=[]; brows=[];joins=[]
def add(p,role,sec='shared'):
 p=Path(p);files.append(dict(section_id=sec,role=role,path=str(p.resolve()),size=p.stat().st_size))
for p,role in [(GENES,'MSigDB 2026.1.Hs'),(P/'input_manifest.tsv','pilot sources'),(P/'reference_sampling_ids.tsv.gz','frozen reference'),(A/'analysis_cell_metadata.tsv.gz','fine labels'),(REF,'raw scRNA CSR'),(R/'mp10_cellchat_high_definition_sensitivity_v1/tables/epi_full_scores_pools.tsv.gz','original H05'),(niche,'original accepted labels'),(F/'tables/input_manifest.tsv','spatial manifest'),(R/'neighborhood_ring_harmony_v3/01_features/feature_manifest.json','original Ring1/Ring2 definition'),(R/'mp10_paper_distance_full_v1/tables/source_R12_profiles.tsv','original peripheral profiles'),(R/'neighborhood_niche_visualization_v1/data_tables/niche_palette.tsv','palette')]:add(p,role)
for s,r in M.iterrows():
 p=P/'spots'/f'{s}.tsv.gz';add(p,'original spots and query',s)
 d=read(p);assert set(d.section_id)=={s} and d['query'].sum()==S[s]
 d=d[KEY+['query',MP,'Epithelial_fraction','nCount_RNA','x','y','array_row','array_col','lesion_location','H05_support','H05_conditional_support','Epi_support','top20_UMI']+[c for c in d if c.startswith('UMI_')]].copy()
 queries.append(d.loc[d['query'],KEY])
 for k in ['counts_path','genes_path','barcodes_path','c2l_mean_path','image_path']:add(r[k],k,s)
 ab=pd.read_csv(r.c2l_mean_path,index_col=0);assert set(ab)==set(mapping.subtype) and ab.index.is_unique and set(ab.index)==set(d.barcode) and (ab>=0).all().all()
 ab=ab.loc[d.barcode];total=ab.sum(axis=1,min_count=30);assert np.allclose(d.Epithelial_fraction,ab.Epithelial/total)
 parts=ab.rename(columns=ccol).reset_index(drop=True);parts['c2l_mean_all30']=total.to_numpy();parts['c2l_mean_TME29']=(total-ab.Epithelial).to_numpy()
 for group,g in mapping[mapping.major!='Epithelial'].groupby('major'):
  parts[ccol('major_'+group)]=ab[g.subtype].sum(axis=1,min_count=len(g)).to_numpy()
 for t in [*ab.columns,*['major_'+k for k in mapping.major.unique() if k!='Epithelial']]:
  parts[lcol(t)]=np.log1p(parts[ccol(t)])
 d=pd.concat([d,parts],axis=1)
 nodep=R/'neighborhood_spottedpy_v8/01_geometry'/s/'nodes.parquet';add(nodep,'original spot_uid and rings',s)
 node=pd.read_parquet(nodep,columns=KEY+['spot_uid','ring1_valid_n','ring2_valid_n','complete_rings'])
 n_extra=len(node)-len(d);d=join(d,node[node.barcode.isin(d.barcode)]);lab=pd.read_parquet(niche,filters=[('section_id','==',s)],columns=['spot_uid','accepted_label','assignment_status'])
 d=d.merge(lab,on='spot_uid',how='left',validate='one_to_one');d['niche']=d.accepted_label.where(d.accepted_label>=0).add(1).astype('Int64')
 d.loc[d.assignment_status.isna(),'assignment_status']='outside_original_anchor_table'
 save(d,f'spots/{s}.base.tsv.gz')
 e=d[d['query']];brows.extend(square(e,s,'B',[MP]+[lcol(t) for t in mapping.subtype if t!='Epithelial']))
 joins.append(dict(section_id=s,full=len(d),query=len(e),outside_frozen_full_nodes=n_extra,assigned_query=e.niche.notna().sum(),unassigned_query=e.niche.isna().sum()))
save(brows,'tables/B_pairwise.tsv.gz');save(pd.concat(queries),'frozen_query.tsv.gz');save(joins,'tables/join_audit.tsv');save(files,'input_manifest.tsv')
print('B complete',joins,flush=True)
with h5py.File(REF) as f:
 ids=pd.Index(strings(f['obs/_index'])); genes=strings(f['var/_index'])
 for model,g in ref.groupby('model'):
  idx=ids.get_indexer(g.barcode);assert min(idx)>=0 and np.array_equal(np.char.add(np.char.add(strings(f['obs/dataset'])[idx].astype(str),'::'),strings(f['obs/orig.ident'])[idx].astype(str)),g.sample_uid)
  pack(model,slice_h5(f,idx),genes,g.barcode)
coverage=[]
for s,r in M.iterrows():
 x=io.mmread(r.counts_path).tocsr();genes=pd.Index(Path(r.genes_path).read_text().splitlines());ids=pd.Index(Path(r.barcodes_path).read_text().splitlines())
 if not genes.is_unique:
  ix,genes=pd.factorize(genes);x=sparse.csr_matrix((np.ones(len(ix),dtype=int),(ix,np.arange(len(ix)))),shape=(len(genes),len(ix)))@x
 d=read(O/f'spots/{s}.base.tsv.gz').set_index('barcode');assert set(d.index)==set(ids) and np.array_equal(np.asarray(x.sum(axis=0)).ravel(),d.loc[ids,'nCount_RNA'])
 pack(s,x.T,genes,ids)
 det=set(genes[np.asarray((x>0).sum(axis=1)).ravel()>0])
 for g in gs:
  members=set(g.gene_symbol);coverage.append(dict(section_id=s,module=g.module.iat[0],variant=g.variant.iat[0],n_total=len(members),n_assayed=len(members&set(genes)),n_detected=len(members&det),detected_genes=','.join(sorted(members&det)),n_background=len(genes),aucMaxRank=int(np.ceil(.05*len(genes))),member_top20=','.join(sorted(members&set(TOP)))))
 # Individual raw marker counts for discrepant regions, including genes outside MP10.
 targets=['ABCA1','ABCG1','DHCR7','DHCR24','LSS','CYP51A1','FASN','CPT1A','ACADVL']
 dex=pd.DataFrame(dict(section_id=s,barcode=ids))
 for t in targets:dex['UMI_'+t]=x[genes.get_loc(t)].toarray().ravel() if t in genes else np.nan
 save(dex,f'spots/{s}.markers.tsv.gz')
save(coverage,'tables/gene_coverage_base.tsv')

