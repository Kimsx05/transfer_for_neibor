from common import *
import anndata as ad
from scipy.stats import spearmanr
from scipy.io import mmwrite
from numba import njit
@njit
def auc_sparse(indptr,indices,data,n_genes,top_idx,seed,keep_zero):
 np.random.seed(seed);n=len(indptr)-1;K=int(np.ceil(.05*n_genes));den=sum([K-i for i in range(1,len(top_idx)+1)])
 out=np.zeros(n);lookup=np.full(n_genes,-1,np.int32)
 for i in range(len(top_idx)):lookup[top_idx[i]]=i
 for j in range(n):
  lo=indptr[j];hi=indptr[j+1];v=data[lo:hi];ii=indices[lo:hi];nz=len(v);order=np.argsort(-v.astype(np.float64)-np.random.random(nz)*.1);rr=np.full(len(top_idx),-1,np.int32)
  for k in range(nz):
   g=lookup[ii[order[k]]]
   if g>=0:rr[g]=k+1
  sm=0.;drawn=np.empty(len(top_idx),np.int32);nd=0
  for k in range(len(top_idx)):
   rank=rr[k]
   if rank<0 and not keep_zero:
    rank=np.random.randint(nz+1,n_genes+1)
    while rank in drawn[:nd]:rank=np.random.randint(nz+1,n_genes+1)
    drawn[nd]=rank;nd+=1
   if rank>0 and rank<K:sm+=K-rank
  out[j]=sm/den
 return out
@njit
def auc_historical(indptr,indices,data,n_genes,top_idx):
 n=len(indptr)-1;K=int(np.ceil(.05*n_genes));den=sum([K-i for i in range(1,len(top_idx)+1)]);out=np.zeros(n);lookup=np.zeros(n_genes,np.int8);lookup[top_idx]=1
 for j in range(n):
  lo=indptr[j];hi=indptr[j+1];v=data[lo:hi];ii=indices[lo:hi];nz=len(v)
  if nz==0:out[j]=np.nan;continue
  order=np.argsort(-v.astype(np.float64));k=0;sm=0.
  while k<nz:
   end=k+1
   while end<nz and v[order[end]]==v[order[k]]:end+=1
   mean_rank=(k+1+end)/2.
   if mean_rank<K:
    for a in range(k,end):
     if lookup[ii[order[a]]]:sm+=K-mean_rank
   k=end
  out[j]=sm/den
 return out
if __name__=='__main__':
 gate();audit=[]
 for unit in ['cell','bin8','bin16']:
  if unit=='cell':
   a=ad.read_h5ad(A/'segmentation_assignment_v2/Cervilla/cells.h5ad');keep=a.obs.qc_keep_analysis.to_numpy();x=sparse.csr_matrix(a.X[keep]);d=a.obs.loc[keep].copy();d['cell_id']=d.index;d['bin_id']=d.index;sy=a.var['symbol'].to_numpy() if 'symbol'in a.var else a.var.iloc[:,0].to_numpy();print('CELLVAR',a.var.columns.tolist(),flush=True);del a
  elif unit=='bin8':x,d,ge,sy,_=read8('Cervilla')
  else:x=sparse.load_npz(R/'models/Cervilla_bin16_counts.npz');d=pd.read_parquet(R/'tables/Cervilla_bin16_measurement.parquet');sy=pd.read_csv(R/'tables/hd_genes.tsv',sep='\t').symbol.to_numpy()
  d=panel(x,d,sy);ti=np.flatnonzero(np.isin(sy,TOP));assert len(ti)==20;d['mp10_historical_auc']=auc_historical(x.indptr,x.indices,x.data,x.shape[1],ti);d['mp10_auc_observed']=auc_sparse(x.indptr,x.indices,x.data,x.shape[1],ti,123,True)
  # Small fixed stratified validation subset by RNA detection only, independent of score/association.
  rng=np.random.default_rng(123);strata=pd.cut(d.nFeature_RNA,[-1,20,100,500,np.inf]);ii=np.concatenate([rng.choice(v,min(150,len(v)),False)for _,v in pd.Series(np.arange(len(d))).groupby(strata.to_numpy(),observed=True)])
  xx=x[ii];rows=[]
  for seed in [123,124,125]:
   s0=auc_sparse(xx.indptr,xx.indices,xx.data,xx.shape[1],ti,seed,False);s1=auc_sparse(xx.indptr,xx.indices,xx.data,xx.shape[1],ti,seed,True)
   rows.append(pd.DataFrame(dict(id=d.iloc[ii].bin_id.to_numpy(),unit=unit,seed=seed,nCount=d.iloc[ii].nCount_RNA.to_numpy(),nFeature=d.iloc[ii].nFeature_RNA.to_numpy(),MP10_detected=d.iloc[ii].MP10_detected.to_numpy(),auc_original=s0,auc_observed=s1)))
  vv=pd.concat(rows);put(vv,f'{unit}_score_seed_validation.tsv.gz')
  for frac in [.5]:
   z=xx.copy();z.data=rng.binomial(np.round(z.data).astype(int),frac);z.eliminate_zeros();ss=auc_sparse(z.indptr,z.indices,z.data,z.shape[1],ti,123,True);base=vv[vv.seed==123];put(pd.DataFrame(dict(id=base.id,unit=unit,full=base.auc_observed,half_depth=ss,nFeature=base.nFeature,MP10_detected=base.MP10_detected)),f'{unit}_score_depth_validation.tsv')
  ro=pd.read_csv(R/'tables/roi_manifest.tsv',sep='\t');d['roi']=''
  for r in ro.itertuples():d.loc[d.x_um.between(r.x0,r.x1,inclusive='left')&d.y_um.between(r.y0,r.y1,inclusive='left'),'roi']=r.roi
  d.to_parquet(R/'tables'/f'Cervilla_{unit}_measured.parquet',index=False)
  # Actual AUCell package receives bounded raw examples for source-level equivalence and ties audit.
  jj=ii[:64];mmwrite(R/'temporary'/f'{unit}_AUCell.mtx',x[jj].T);put(pd.DataFrame(dict(id=d.iloc[jj].bin_id.to_numpy())) ,f'{unit}_AUCell_ids.tsv');pd.Series(sy).to_csv(R/'temporary'/f'{unit}_genes.txt',index=False,header=False)
  audit.append(dict(unit=unit,n=len(d),UMI=int(x.sum()),genes=x.shape[1],top20_covered=len(ti),MP10_detected_median=float(d.MP10_detected.median()),aucMaxRank=int(np.ceil(.05*x.shape[1]))));print('MEASURED',unit,audit[-1],flush=True);del x,d
 put(pd.DataFrame(audit),'measurement_summary.tsv')
