from common import *
import anndata as ad
from scipy.io import mmwrite
from scipy.spatial import cKDTree
from datetime import datetime
(R/'temporary').mkdir(exist_ok=True)
# Freeze endpoints and selection rules before any association is calculated.
cfg=dict(version=1,created=datetime.now().isoformat(),sections={'Cervilla':'full','D1_Johnson':'light_target_evaluability'},seed=123,validation_seeds=[123,124,125],threads=4,night_window='23:00–07:00 Asia/Shanghai',ROI_selection='5 inherited; 256um squares centred on strongest unsmoothed T-panel and endothelial-panel block, plus nonimmune measured control; no MP10 relationship',reference='final30 exact membership, complete lineage background',reference_cap_per_type=600,epi_reference_cap=1800,holdout='every fifth sorted orig.ident; sample not patient holdout',rctd=dict(mode='doublet',UMI_min=100,UMI_min_sigma=300,counts_MIN=10,max_cores=4),cell_annotation='SingleR wilcox aggregate reference, no TIGIT; all genes shared before DE selection; lineage support and fine subtype uncertainty separate',radii_um=[50,200],niche_K=[2,4,6,8],niche_seeds=[123,456,789],niche_selection='stable per-cluster match and interpretable composition; K2 retained; choose smallest finer K with seedARI>=0.8 and blockARI>=0.6 and minimum matched-cluster Jaccard>=0.4; otherwise coarse',MP10=dict(genes=TOP,original_aucMaxRank_fraction=.05,zero_ties='random seed audit; keepZeroesAsNA sensitivity',quantile_type=7,groups=['L50','H50','H10','H05'],eligibility='credible Epi; >=5 detected Top20; rank/depth validation; if poor single-unit precision use fixed Epi-only 64um grid pooled raw counts, >=5 credible Epi units, not crossing input components'),association=dict(primary='continuous MP10',block_um=400,bootstrap=499,min_blocks=12,continuous_BH_family='Treg support, Treg/CD4 or Treg/T, Treg distance, Endothelial support and distance at 50/200; reverse comparison',niche_family='niche MP10 contrasts at both radii',conditions='log1p nCount, Epi support/purity, observed support; additional vessel distance; no collinear nFeature together with nCount'),frozen_before_association=True)
if not (R/'frozen_config.yaml').exists():(R/'frozen_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False,allow_unicode=True))
refpath=P/'subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad';a=ad.read_h5ad(refpath,backed='r');md=a.obs.copy();md['cell_id']=md.index;md['label']=md.spatial_reference_celltype_v2_final.astype(str);md['sample']=md['orig.ident'].astype(str);hold=set(sorted(md['sample'].unique())[::5]);rng=np.random.default_rng(123);ids=[];valid=[]
for typ,g in md.groupby('label',observed=True):
 for sample,v in g.groupby('sample',observed=True):
  if sample in hold:valid.extend(rng.choice(v.index,min(8,len(v)),replace=False))
 train=g[~g['sample'].isin(hold)];cap=1800 if typ=='Epithelial' else 600
 # round robin samples; capped type contribution independent of abundance
 buckets=[list(rng.permutation(v.index)) for _,v in train.groupby('sample',observed=True)];chosen=[]
 while len(chosen)<cap and any(buckets):
  for b in buckets:
   if b and len(chosen)<cap:chosen.append(b.pop())
 ids.extend(chosen)
sel=ids+valid;z=md.loc[sel].copy();z['partition']=['train']*len(ids)+['sample_holdout']*len(valid);full=a.X.to_memory();xx=full[md.index.get_indexer(sel)].tocsr();del full;print('REFERENCE_SLICE',xx.shape,flush=True);assert (xx.data>=0).all() and np.allclose(xx.data,np.round(xx.data));mmwrite(R/'temporary/reference.mtx',xx.T);put(z.reset_index(drop=True),'reference_membership.tsv.gz');put(a.var.reset_index(names='symbol'),'reference_genes.tsv');a.file.close()
manifest=[dict(role='fine_reference_raw',path=str(refpath),size_bytes=refpath.stat().st_size)]
roi=pd.read_csv(A/'hd_epi_niche_multimethod_v2/tables/fixed_ROI_index.tsv',sep='\t');print('ROI cols',roi.columns.tolist(),flush=True)
allrois=[];geneaudit=[]
for s in ['Cervilla','D1_Johnson']:
 x,d,ge,sy,path=read8(s);manifest.append(dict(role=s+'_native8_raw',path=str(path),size_bytes=path.stat().st_size));d=panel(x,d,sy);d.to_parquet(R/'tables'/f'{s}_bin8_measurement.parquet',index=False)
 geneaudit.extend(dict(section=s,panel=k,gene=g,present=bool(np.any(sy==g)),detected_bins=int((x[:,sy==g].sum(1)>0).sum())) for k,v in {**PANELS,'MP10':TOP,'readout':['TIGIT']}.items() for g in v)
 if s=='Cervilla':
  rr=roi[roi.section.eq(s)].copy();print(rr.to_string(index=False),flush=True)
  # names below resolved from actual recorded columns
  for row in rr.to_dict('records'):
   allrois.append(dict(section=s,roi=row.get('roi',row.get('name')),x0=row['x'],y0=row['y'],x1=row['X'],y1=row['Y'],reason='inherited pre-MP10 ROI'))
  block=d.assign(bx=np.floor(d.x_um/256).astype(int),by=np.floor(d.y_um/256).astype(int));g=block[block.analysis_qc].groupby(['bx','by']).agg(n=('bin_id','size'),T=('T_UMI','sum'),Endo=('Endo_UMI','sum'),UMI=('nCount_RNA','sum'))
  for name,col in [('T_lineage_focus','T'),('vascular_focus','Endo')]:
   bx,by=g[g.n>=250][col].idxmax();allrois.append(dict(section=s,roi=name,x0=bx*256,y0=by*256,x1=(bx+1)*256,y1=(by+1)*256,reason='maximum raw panel UMI in measured 256um block with >=250 QC bins; no MP10 used'))
  control=g[(g.n>=500)&(g.T/g.UMI<.002)] if False else g[(g.n>=500)&(g['T']/g.UMI<.002)]
  bx,by=control.UMI.idxmax();allrois.append(dict(section=s,roi='nonimmune_control',x0=bx*256,y0=by*256,x1=(bx+1)*256,y1=(by+1)*256,reason='high measured coverage, raw T-panel UMI fraction<0.002; no MP10 used'))
  # Fixed origin row/col aggregation from all native supplied 8um, not QC subset.
  keys=list(zip((d.array_row//2).astype(int),(d.array_col//2).astype(int)));lev=sorted(set(keys));idx=pd.Index(lev).get_indexer(keys);G=sparse.csr_matrix((np.ones(len(d),dtype=np.int64),(idx,np.arange(len(d)))),shape=(len(lev),len(d)));y=G@x;assert np.array_equal(np.asarray(y.sum(0)),np.asarray(x.sum(0)));sparse.save_npz(R/'models/Cervilla_bin16_counts.npz',y)
  # Verified affine grid transform to fixed full-square centers, no observed centroid shift.
  X=np.c_[np.ones(len(d)),d.array_row,d.array_col];coef=np.linalg.lstsq(X,d[['x_um','y_um']].to_numpy(),rcond=None)[0];rc=np.array(lev);xy=np.c_[np.ones(len(lev)),2*rc[:,0]+.5,2*rc[:,1]+.5]@coef
  m=pd.DataFrame(dict(bin_id=['s_016um_%05d_%05d-1'%q for q in lev],array_row=rc[:,0],array_col=rc[:,1],x_um=xy[:,0],y_um=xy[:,1],supplied8=np.asarray(G.sum(1)).ravel(),nCount_RNA=np.asarray(y.sum(1)).ravel(),nFeature_RNA=np.asarray((y>0).sum(1)).ravel(),in_tissue_supplied8=np.asarray(G@d.in_tissue.to_numpy()).ravel()))
  m['qc_keep']=(m.in_tissue_supplied8>0)&(m.nCount_RNA>5);m['partial_supplied']=m.supplied8<4;m=panel(y,m,sy);m.to_parquet(R/'tables/Cervilla_bin16_measurement.parquet',index=False);put(pd.DataFrame(dict(gene_id=ge,symbol=sy,input_UMI=np.asarray(x.sum(0)).ravel(),aggregated_UMI=np.asarray(y.sum(0)).ravel())),'bin16_gene_conservation.tsv');put(pd.DataFrame(dict(gene_id=ge,symbol=sy)),'hd_genes.tsv')
  # R temporary sparse input is restricted to fixed ROIs and is removed after conversion.
  r=pd.DataFrame(allrois);assert r[['x0','y0','x1','y1']].notna().all().all()
  for unit,mat,meta in [('bin8',x,d),('bin16',y,m)]:
   flags=np.zeros(len(meta),bool);meta['roi']=''
   for v in r.itertuples():
    hit=meta.x_um.between(v.x0,v.x1,inclusive='left')&meta.y_um.between(v.y0,v.y1,inclusive='left');flags|=hit;meta.loc[hit,'roi']=v.roi
   take=np.flatnonzero(flags);mmwrite(R/'temporary'/f'{unit}_ROI.mtx',mat[take].T);put(meta.iloc[take],f'{unit}_ROI.tsv.gz')
 del x,d
put(pd.DataFrame(allrois),'roi_manifest.tsv');put(pd.DataFrame(geneaudit),'gene_panel_coverage.tsv');put(pd.DataFrame(manifest),'input_manifest.tsv');print('PREPARED',flush=True)
