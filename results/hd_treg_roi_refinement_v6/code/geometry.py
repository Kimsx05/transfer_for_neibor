from base6 import *
from numba import njit
from classifier import predict,to_ref

@njit
def nearest_two(mask,xx,yy,valid):
 n=len(xx);ids=np.zeros(n,np.int32);nuc=ids.copy();second=ids.copy()
 d1=np.full(n,100.,np.float32);d2=d1.copy()
 for i in range(n):
  if not valid[i]:continue
  x=xx[i];y=yy[i];cx=round(x);cy=round(y)
  if 0<=cy<mask.shape[0] and 0<=cx<mask.shape[1]:nuc[i]=mask[cy,cx]
  for v in range(max(0,cy-10),min(mask.shape[0],cy+11)):
   for u in range(max(0,cx-10),min(mask.shape[1],cx+11)):
    lab=mask[v,u]
    if lab==0:continue
    di=(u-x)**2+(v-y)**2
    if lab==ids[i]:
     if di<d1[i]:d1[i]=di
    elif di<d1[i]:
     second[i]=ids[i];d2[i]=d1[i];ids[i]=lab;d1[i]=di
    elif di<d2[i]:
     second[i]=lab;d2[i]=di
 return ids,nuc,second,np.sqrt(d1)*.5,np.sqrt(d2)*.5
def process(s):
 o=R/s;t=js(N/s/'metadata/transform.json');scale=np.array(t['local_to_model_scale_xy']);org=np.array(t['fullres_crop_xyxy'][:2]);mpp=t['native_mpp']
 new=sparse.load_npz(o/'nuclei_v6.npz').toarray().astype('int32',copy=False)
 present=np.unique(new); compact=np.zeros(int(present.max())+1,np.int32); positive=present[present>0]; compact[positive]=np.arange(1,len(positive)+1); new=compact[new]
 old=sparse.load_npz(S/s/'inputs/nuclei.npz').tocsr();m=int(new.max());n=int(old.max())
 oa=np.bincount(old.data,minlength=n+1);na=np.bincount(new.ravel(),minlength=m+1)
 pairs=[];y_sum=np.zeros(m+1);x_sum=y_sum.copy()
 for a in range(0,new.shape[0],256):
  block=new[a:a+256];y,x=np.nonzero(block);ids=block[y,x];y_sum+=np.bincount(ids,weights=y+a,minlength=m+1);x_sum+=np.bincount(ids,weights=x,minlength=m+1)
  oo=old[a:a+256].tocoo();v=block[oo.row,oo.col];ok=v>0;k,count=np.unique(oo.data[ok].astype('int64')*(m+1)+v[ok],return_counts=True);pairs.append(pd.DataFrame(dict(old_id=k//(m+1),new_id=k%(m+1),intersection=count)))
 pair=pd.concat(pairs).groupby(['old_id','new_id']).intersection.sum().reset_index()
 pair['overlap_old']=pair.intersection/oa[pair.old_id];pair['overlap_new']=pair.intersection/na[pair.new_id]
 pair['substantial']=(pair.intersection>=4)&(pair[['overlap_old','overlap_new']].max(axis=1)>=.25)
 sub=pair[pair.substantial];mappedold=sub.groupby('old_id').size();mappednew=sub.groupby('new_id').size()
 od=qc(s);uid=od.id.to_numpy()
 # Only restore a previously accepted old object with no touching new mask and contrast against its local surround.
 no_overlap=~od.id.isin(pair.old_id)&~od.geometry_review&~od.image_confirmed_invalid
 bg=np.load(N/s/'inputs/model_image.npy',mmap_mode='r');px=np.rint((od[['x','y']].to_numpy()/mpp-org+.5)*scale-.5).astype(int)
 rescue=[];dark=[]
 for j in np.flatnonzero(no_overlap):
  cx,cy=px[j];y0=max(0,cy-24);x0=max(0,cx-24);sl=np.s_[y0:cy+25,x0:cx+25];om=old[sl].toarray();hit=om==j+1;ring=ndimage.binary_dilation(hit,iterations=4)&(om==0)&(new[sl]==0)
  if hit.sum()<12 or ring.sum()<8:continue
  # H&E nuclei locally darker in green channel than cytoplasm/background.
  contrast=np.median(bg[sl][:,:,1][ring])-float(np.median(bg[sl][:,:,1][hit]))
  if contrast>=15:rescue.append(j+1);dark.append(contrast)
 resmap=np.zeros(n+1,np.int32);resmap[rescue]=np.arange(len(rescue))+m+1
 if rescue:
  for a in range(0,new.shape[0],256):
   oo=old[a:a+256].tocoo();v=resmap[oo.data];k=v>0;new[a+oo.row[k],oo.col[k]]=v[k]
  p=pd.DataFrame(dict(old_id=rescue,new_id=resmap[rescue],intersection=oa[rescue],overlap_old=1.,overlap_new=1.,substantial=True));pair=pd.concat([pair,p],ignore_index=True)
 oldout=pd.DataFrame(dict(old_id=uid,nuclear_area_pixels=oa[1:],v6_overlaps=od.id.map(mappedold).fillna(0).astype(int)))
 oldout['action']=np.select([resmap[uid]>0,oldout.v6_overlaps>1,oldout.v6_overlaps==1],['retain_old_visible_unmatched','split_or_partial_split','overlap_with_redetection'],default='not_retained_by_redetection')
 oldout.to_parquet(o/'nuclear_audit_old.parquet',index=False);pair.to_parquet(o/'nuclear_overlap.parquet',index=False)
 ids=np.arange(1,m+1)
 centers=(np.c_[x_sum[1:],y_sum[1:]]/na[1:,None]+.5)/scale-.5+org
 centers=centers*mpp
 cells=pd.DataFrame(dict(id=ids,x=centers[:,0],y=centers[:,1],nuclear_area_um2=na[1:]*np.prod(t['actual_model_mpp_xy']),overlap_old=pd.Series(ids).map(mappednew).fillna(0).astype(int)))
 cells['nuclear_action']=np.select([cells.overlap_old==0,cells.overlap_old>1],['new_detection','merge_or_partial_merge'],default='redetected')
 if rescue:
  rr=od[od.id.isin(rescue)][['id','x','y','nuclear_area_um2']].copy();rr.id=resmap[rr.id];rr['nuclear_action']='retain_old_visible_unmatched';rr['overlap_old']=1;cells=pd.concat([cells,rr],ignore_index=True)
 assert np.array_equal(cells.id,np.arange(1,len(cells)+1))
 sparse.save_npz(o/'final_nuclei.npz',sparse.csr_matrix(new))
 bins,audit=inputdata(s);xy=bins[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();p=(xy-org+.5)*scale-.5
 first,nuc,other,dist,dist2=nearest_two(new,p[:,0],p[:,1],bins.usable.to_numpy());del new
 ambiguous=(nuc==0)&(other>0)&(dist<=4)&(dist2-dist<1)
 hard=np.where(nuc>0,nuc,np.where((dist<=4)&~ambiguous,first,0)).astype('int32')
 mp=bins[['source_bin_index','array_row','array_col','usable']].copy();mp['nucleus_id']=nuc;mp['id']=hard;mp['nearest_id']=first;mp['runner_id']=other;mp['distance_um']=dist;mp['ambiguous']=ambiguous;mp.to_parquet(o/'assignment.parquet',index=False)
 n=len(cells);genes,sy=genedata(audit);ng=len(sy);X=sparse.csr_matrix((n,ng),dtype=np.int64);acc=np.zeros((4,ng),np.int64);use=bins.usable.to_numpy()
 for a,b,z in slabs(audit):
  v=hard[a:b];ii=np.flatnonzero(v);mat=sparse.csr_matrix((np.ones(len(ii),int),(v[ii]-1,ii)),shape=(n,len(v)));X+=mat@z
  for j,f in enumerate([np.ones(len(v),bool),~use[a:b],ambiguous[a:b],(v==0)&use[a:b]&~ambiguous[a:b]]):acc[j]+=np.asarray(z[f].sum(0)).ravel()
 assigned=np.asarray(X.sum(0)).ravel();assert np.array_equal(acc[0],acc[1:].sum(0)+assigned)
 put(pd.DataFrame(dict(gene=genes,total=acc[0],excluded=acc[1],ambiguous=acc[2],unassigned=acc[3],assigned=assigned)),o/'count_balance.tsv.gz')
 sparse.save_npz(o/'counts.npz',X)
 cells['nCount']=np.asarray(X.sum(1)).ravel();cells['nFeature']=X.getnnz(1);cells['nuclear_bins']=np.bincount(nuc,minlength=n+1)[1:]
 for g in T+reg+epi+fib+mye:cells[g]=np.asarray(X[:,sy==g].sum(1)).ravel()
 for key,gs in [('direct',direct),('other_core',core[1:]),('T',T),('Epi',epi),('Fib',fib),('Myeloid',mye)]:
  cells[key+'_width']=cells[gs].gt(0).sum(1)
 model=joblib.load(R/'model.joblib');q=predict(model,to_ref(X,sy,model[0]));cells=pd.concat([cells,q],axis=1)
 # Establish all molecular hints including raw unassigned signal before scoring.
 raw=pd.read_parquet(V5/s/'raw_candidates.parquet');raw['id']=hard[raw.source_bin_index];raw['nucleus_id']=nuc[raw.source_bin_index];raw['ambiguous']=ambiguous[raw.source_bin_index]
 raw.to_parquet(o/'raw_signal.parquet',index=False)
 cells.to_parquet(o/'cells_geometry.parquet',index=False)
 print(s,len(cells),len(rescue),int(assigned.sum()),flush=True)
if __name__=='__main__':
 for s in sys.argv[1:] or sections:process(s)
