from pathlib import Path
import sys,json,gzip
import numpy as np,pandas as pd,h5py,pyarrow as pa,anndata as ad
from scipy import sparse
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path(__file__).resolve().parents[1];ROOT=O.parents[1];RETRY=O.parent/'hd_xenium_annotation_niche_retry_v1'
S={'hd':['Cervilla','D1_Johnson'],'xenium':['A1R2','A1R4','A1R5','A1R6','A2R2']}
R=['0.2','0.5','1.0']
def cname(r):return 'leiden_r'+r.replace('.','_')
def write(frame,path):frame.to_csv(path,sep='\t',index=False,na_rep='NA')
def read_sample(p,s,save=False):
 if p=='hd':
  d=RETRY/'hd'/s;f=d/'cells.h5ad';a=ad.read_h5ad(f);x=a.X.tocsr();q=a.obs.copy();ids=q.index
  g=pd.DataFrame({'feature_id':a.var.feature_id.to_numpy(),'gene_symbol':a.var.symbol.to_numpy(),'unique_name':a.var_names});controls={}
  supp=pd.read_parquet(d/'secondary_image_support.parquet')
  assert supp.index.is_unique and set(supp.index)==set(ids[q.source=='GEX_secondary'])
  q['secondary_image_suspect']=supp.image_suspect.reindex(ids).eq(True);q['base_geometry_flag']=q.geometry_flag
  rules='n_counts>=10 and nFeature>=10 and not base_geometry_flag and not secondary_image_suspect'
  coords=str(f)+' /obs/x_um,y_um';store='/X CSR; integer raw sum of unique assigned 2um bins';idfield='/obs/_index';q['area_outlier']=False;q['boundary_degenerate']=False
 else:
  paths=pd.read_csv(O.parent/'mp10_hd_xenium_colocalization_v1/shared/xenium_inputs.tsv',sep='\t');ps=paths[paths['sample']==s].set_index('role').path
  f=Path(ps['expression']);coords=ps['cell centroids + QC'];store='/matrix CSC features x cells; Gene Expression only';idfield='/matrix/barcodes'
  with h5py.File(f) as h:
   m=h['matrix'];fg=m['features'];ids=pd.Index(m['barcodes'].asstr()[:]);ft=fg['feature_type'].asstr()[:];sel=ft=='Gene Expression'
   allx=sparse.csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:][::-1]));controls={t:{'n_features':int((ft==t).sum()),'counts':int(allx[:,ft==t].sum())} for t in np.unique(ft[~sel])}
   x=allx[:,sel].tocsr();g=pd.DataFrame({'feature_id':fg['id'].asstr()[:][sel],'gene_symbol':fg['name'].asstr()[:][sel]});g['unique_name']=g.feature_id
   if save:
    cf=pd.DataFrame({'platform':p,'sample':s,'feature_id':fg['id'].asstr()[:][~sel],'name':fg['name'].asstr()[:][~sel],'feature_type':ft[~sel],'n_counts':allx[:,~sel].sum(axis=0).A1})
    write(cf,O/p/(s+'_controls.tsv'))
  q=pd.read_csv(coords).set_index('cell_id');assert q.index.is_unique and set(q.index)==set(ids);q=q.loc[ids].rename(columns={'x_centroid':'x_um','y_centroid':'y_um','cell_area':'area_um2'})
  with gzip.open(ps['original cell polygons'],'rb') as stream:
   vertices=pd.read_parquet(stream,columns=['cell_id']).groupby('cell_id').size()
  q['boundary_vertices']=vertices.reindex(ids).fillna(0);q['boundary_degenerate']=q.boundary_vertices<4
  q['area_outlier']=(q.area_um2<=q.area_um2.quantile(.005))|(q.area_um2>=q.area_um2.quantile(.995))
  q['base_geometry_flag']=False;q['secondary_image_suspect']=False;q['source']='original_segmentation';rules='nFeature>=10; vertices>=4; area strictly between per-section P0.5 and P99.5'
 assert ids.is_unique and g.feature_id.is_unique and (x.data>=0).all() and np.equal(x.data,np.floor(x.data)).all()
 x.sum_duplicates();x.eliminate_zeros();lib=x.sum(1).A1;nf=x.getnnz(1)
 if p=='hd':assert np.array_equal(q.n_counts,lib) and np.array_equal(q.nFeature,nf)
 else:assert np.array_equal(q.transcript_counts,lib)
 q['n_counts']=lib;q['nFeature']=nf;q['low_features']=nf<10;q['low_counts']=(lib<10)&(p=='hd')
 flags=['low_features','low_counts','base_geometry_flag','secondary_image_suspect','area_outlier','boundary_degenerate']
 q['qc_pass']=~q[flags].any(axis=1);q['qc_reasons']=['|'.join(k for k,test in zip(flags,row) if test) for row in q[flags].to_numpy()]
 q['platform']=p;q['sample']=s;q['original_cell_id']=ids;q['cell_key']=s+'::'+ids;q['source_row']=np.arange(len(q))
 info=dict(platform=p,sample=s,counts_file=str(f),counts_location=store,objects=len(q),genes=len(g),input_id=idfield,coordinates=coords,units='um',qc=rules,included=int(q.qc_pass.sum()),excluded=int((~q.qc_pass).sum()),counts_dtype=str(x.dtype),input_total_GEX=int(lib.sum()),controls=json.dumps(controls),**{k:int(q[k].sum()) for k in flags})
 return x,g,q,info
def combine(p,save=False):
 xx=[];qq=[];genes=None;manifest=[]
 for s in S[p]:
  x,g,q,m=read_sample(p,s,save=save);manifest.append(m)
  if genes is None:genes=g.copy()
  else:
   order=pd.Index(g.feature_id).get_indexer(genes.feature_id)
   if len(g)!=len(genes) or min(order)<0:raise ValueError('Measured feature sets differ: implement measured-only strata before pooling')
   assert np.array_equal(g.gene_symbol.to_numpy()[order],genes.gene_symbol)
   x=x[:,order]
  ok=q.qc_pass.to_numpy();xx.append(x[ok]);qq.append(q.loc[ok]);del x
  if save:q.reset_index(drop=True).to_parquet(O/p/(s+'_qc.parquet'),index=False)
 q=pd.concat(qq,ignore_index=True);assert q.cell_key.is_unique
 if save:
  write(pd.DataFrame(manifest),O/p/'input_manifest.tsv');write(genes,O/p/'features.tsv');q.to_parquet(O/p/'qc_eligible.parquet',index=False)
 return sparse.vstack(xx,format='csr'),genes,q

