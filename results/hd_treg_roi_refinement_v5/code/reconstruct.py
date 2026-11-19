from common5 import *
from numba import njit
from scipy import ndimage
import anndata as ad
@njit
def relax(labels,valid,rgb,logdepth):
 """Five synchronous relaxations suffice for a ten micron budget at two micron pitch."""
 rows,cols=labels.shape
 dist=np.where(labels>0,0.,np.inf).astype(np.float32)
 for iteration in range(5):
  out=labels.copy();cost=dist.copy()
  for i in range(1,rows-1):
   for j in range(1,cols-1):
    if not valid[i,j] or dist[i,j]==0:continue
    for di in (-1,0,1):
     for dj in (-1,0,1):
      if di==dj==0:continue
      u=i+di;v=j+dj;n=labels[u,v]
      if n==0:continue
      colour=0.
      for h in range(3):colour+=abs(float(rgb[i,j,h])-float(rgb[u,v,h]))
      edge=2*np.sqrt(di*di+dj*dj)*(1+2*colour/765+.15*abs(logdepth[i,j]-logdepth[u,v]))
      new=dist[u,v]+edge
      if new<=10 and (new<cost[i,j]-1e-5 or (abs(new-cost[i,j])<=1e-5 and n<out[i,j])):
       cost[i,j]=new;out[i,j]=n
  dist=cost;labels=out
 return labels,dist
def process(s):
 out=R/s;d,audit=inputdata(s);ids,sy=genedata(audit);cells=qc(s);x=d.array_col.to_numpy().astype(int);y=d.array_row.to_numpy().astype(int);off=np.array([x.min()-1,y.min()-1]);x-=off[0];y-=off[1];shape=(y.max()+2,x.max()+2)
 valid=np.zeros(shape,bool);usable=d.usable.to_numpy();valid[y,x]=usable;nuc=np.zeros(shape,np.int32);nuc[y,x]=np.where(usable,d.nucleus_label,0)
 tr=js(N/s/'metadata/transform.json');scale=np.array(tr['local_to_model_scale_xy']);org=np.array(tr['fullres_crop_xyxy'][:2]);xy=d[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()
 pix=np.rint((xy-org+.5)*scale-.5).astype(int);im=np.load(N/s/'inputs/model_image.npy',mmap_mode='r');ok=(pix>=0).all(axis=1)&(pix[:,0]<im.shape[1])&(pix[:,1]<im.shape[0]);rgb=np.zeros((*shape,3),np.uint8);rgb[y[ok],x[ok]]=im[pix[ok,1],pix[ok,0],:3]
 depth=np.zeros(shape,np.float32);depth[y,x]=np.log1p(d.raw_counts)
 # Engineering pilot on one registered historical window per section; full map checks its core against the tiled run.
 q=rd(U/'roi_manifest.tsv').query('section==@s').iloc[0]
 hit=(d.x.between(q.x0,q.x1)&d.y.between(q.y0,q.y1)).to_numpy()
 ya,yb=y[hit].min(),y[hit].max()+1;xa,xb=x[hit].min(),x[hit].max()+1
 sl=np.s_[ya-8:yb+8,xa-8:xb+8];pilot,_=relax(nuc[sl],valid[sl],rgb[sl],depth[sl]);assert np.array_equal(pilot[8:-8,8:-8][nuc[ya:yb,xa:xb]>0],nuc[ya:yb,xa:xb][nuc[ya:yb,xa:xb]>0])
 lab=np.zeros(shape,np.int32);cost=np.full(shape,np.inf,np.float32)
 for r in range(0,shape[0],512):
  for c in range(0,shape[1],512):
   r2=min(r+512,shape[0]);c2=min(c+512,shape[1]);a=max(0,r-8);b=max(0,c-8);sl=np.s_[a:r2+8,b:c2+8];ll,dd=relax(nuc[sl],valid[sl],rgb[sl],depth[sl]);inside=np.s_[r-a:r2-a,c-b:c2-b];lab[r:r2,c:c2]=ll[inside];cost[r:r2,c:c2]=dd[inside]
 assert np.array_equal(lab[ya:yb,xa:xb],pilot[8:-8,8:-8])
 final=lab[y,x];assert np.array_equal(final[d.nucleus_label.gt(0)&usable],d.nucleus_label[d.nucleus_label.gt(0)&usable])
 sparse.save_npz(out/'inferred_cell_grid.npz',sparse.csr_matrix(lab));np.save(out/'grid_origin_col_row.npy',off)
 mp=d[['source_bin_index','array_row','array_col','cell_label','nucleus_label','usable']].copy();mp['new_id']=final;mp['path_cost']=cost[y,x];mp.to_parquet(out/'assignments.parquet',index=False)
 n=len(cells);ng=len(sy);X=sparse.csr_matrix((n,ng),dtype=np.int64);allcounts=np.zeros(ng,np.int64);excluded=allcounts.copy();unassigned=allcounts.copy()
 for a,b,z in slabs(audit):
  zz=final[a:b];keep=zz>0;j=np.flatnonzero(keep)
  ass=sparse.coo_matrix((np.ones(len(j),np.int32),(zz[j]-1,j)),shape=(n,len(keep))).tocsr()
  X=X+ass@z;allcounts+=np.asarray(z.sum(0)).ravel();unassigned+=np.asarray(z[(~keep)&usable[a:b]].sum(0)).ravel();excluded+=np.asarray(z[~usable[a:b]].sum(0)).ravel()
 assert np.array_equal(allcounts,excluded+unassigned+np.asarray(X.sum(0)).ravel())
 cells['nCount']=np.asarray(X.sum(1)).ravel();cells['nFeature']=X.getnnz(1);cells['n_bins']=np.bincount(final,minlength=n+1)[1:];cells['nuclear_bins']=np.bincount(nuc[y,x],minlength=n+1)[1:]
 for g in T+reg+epi+fib+mye:cells[g]=np.asarray(X[:,sy==g].sum(1)).ravel()
 for name,genes in [('T',T),('core',core),('reg',reg),('Epi',epi),('Fib',fib),('Myeloid',mye)]:
  cells[name+'_width']=cells[genes].gt(0).sum(axis=1);cells[name+'_UMI']=cells[genes].sum(axis=1)
 cells=pd.concat([cells,classify(X,sy)],axis=1)
 mixed=(cells.T_width>=2)&((cells.Epi_width>=2)|(cells.Fib_width>=2)|(cells.Myeloid_width>=2))
 bad=cells.image_confirmed_invalid|cells.geometry_review|cells.nuclear_bins.eq(0)
 pred=cells.pred.fillna('');treg=pred.eq('Treg')&(cells.T_width>=2)&(cells.core_width>=1)&((cells.support>=.95)|((cells.FOXP3>0)&(cells.core_width>=2)&(cells.support>=.5)))
 ctrl=pred.str.match('^(CD4|CD8|Tfh)')&(cells.support>=.5)&(cells.T_width>=2)&(cells.FOXP3==0)&(cells.core_width<2)
 cells['status']=np.select([bad,mixed,treg,ctrl,pred.eq('Treg')|cells.FOXP3.gt(0),cells.nCount.lt(25)],['geometry_uncertain','mixed','supported_Treg','supported_other_T','Treg_hint','low_information'],default='other')
 cells['mixed']=mixed;cells['epithelial']=(~bad)&(~mixed)&pred.eq('Epithelial')&(cells.support>=.5)&(cells.Epi_width>=2)
 cells['T_hint']=ist(pred)|cells.T_UMI.gt(0)|cells.reg_UMI.gt(0)
 cells.to_parquet(out/'all_cells.parquet',index=False);sparse.save_npz(out/'counts.npz',X);put(pd.DataFrame(dict(id=ids,symbol=sy)),out/'genes.tsv')
 put(pd.DataFrame(dict(gene=ids,input=allcounts,unassigned=unassigned,excluded=excluded,assigned=np.asarray(X.sum(0)).ravel())),out/'counts_audit.tsv.gz')
 old=ad.read_h5ad(S/s/'cells.h5ad',backed='r');assert np.array_equal(old.obs_names,cells.cell_id)and np.array_equal(old.var_names,ids)
 hist=[]
 for a in range(0,n,10000):hist.append(classify(old.X[a:a+10000],sy))
 old.file.close();p=pd.concat(hist,ignore_index=True);z=cells.loc[ist(p.pred),['id','x','y','cell_id']].copy();z['pred']=p.loc[z.index,'pred'];z['source']='full_v4_model_on_old_cells';z.to_parquet(out/'old_model_candidates.parquet',index=False)
 print(s,'nuclei',n,'assigned',int(X.sum()),'Treg',cells.status.value_counts().to_dict(),'changes',int(np.sum(final!=d.cell_label)),flush=True)
if __name__=='__main__':
 for s in sections:
  if not (R/s/'old_model_candidates.parquet').exists():process(s)
