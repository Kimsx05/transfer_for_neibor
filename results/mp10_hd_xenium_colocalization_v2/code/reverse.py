from common import *
from spatial import groups,ratio
from numba import njit
@njit
def calc(H,grp,ptr,ix,dist,limit,nrep,seed):
 np.random.seed(seed);order=np.argsort(grp);off=np.concatenate((np.zeros(1,np.int64),np.cumsum(np.bincount(grp))));move=np.arange(len(H));N=len(ptr)-1
 out=np.full((nrep+1,H.shape[1],3),np.nan);weight=np.zeros(len(H));ok=0
 for i in range(N):
  n=ptr[i+1]-ptr[i]
  if n:
   ok+=1
   for j in range(ptr[i],ptr[i+1]):weight[ix[j]]+=1./n
 if not ok:return out
 for it in range(nrep+1):
  if it:
   for g in range(len(off)-1):
    for u in range(off[g+1]-1,off[g],-1):
     a=order[u];b=order[np.random.randint(off[g],u+1)];move[a],move[b]=move[b],move[a]
  for k in range(H.shape[1]):
   h=np.zeros(len(H),np.bool_);frac=0.
   for i in range(len(H)):
    if H[i,k]:h[move[i]]=True;frac+=weight[move[i]]
   nh=0;total=0.
   for i in range(N):
    if ptr[i]==ptr[i+1]:continue
    ds=limit
    for j in range(ptr[i],ptr[i+1]):
     if h[ix[j]]:ds=dist[j];nh+=1;break
    total+=ds
   out[it,k,0]=frac/ok;out[it,k,1]=nh/ok;out[it,k,2]=total/ok
 return out
def run(s):
 hd=s in HD;v,d=folders(s);z=labels(s);i=pd.read_parquet(d/'identity.parquet');cov=pd.read_parquet(d/'target_cover.parquet');out=[];loo=[];methods=['MP10'] if hd else STATES
 for res in ([16,32] if hd else ['cell']):
  e=pd.read_parquet(d/f'centers_{res}.parquet');H=e[[k+'_high' for k in methods]].to_numpy(bool)
  for r in ([50,200] if hd else [50,200,1,2,5]):
   if r>=50:G=sparse.load_npz(d/f'epi_{res}_{r}.npz')
   else:
    p=pd.read_parquet(d/'boundary_pairs.parquet');p=p[p.gap_um<=r];G=sparse.csr_matrix((np.ones(len(p),np.int32),(p.row,p.neighbor)),shape=(len(e),len(z)))
   for name,grades in [('supported',['supported']),('broad',['supported','possible'])]:
    tg=i.grade.isin(grades);tgt=np.flatnonzero(tg);keep=cov.loc[z.index[tgt],f'coverage_{max(50,r)}'].values>=.5
    t=tgt[keep];co=G[:,t].T.tocoo();ds=np.linalg.norm(z[['x_um','y_um']].values[t[co.row]]-e[['x_um','y_um']].values[co.col],axis=1);order=np.lexsort((ds,co.row));ix=co.col[order];dist=ds[order];ptr=np.r_[0,np.cumsum(np.bincount(co.row,minlength=len(t)))]
    vals=calc(H,groups(e,400),ptr,ix,dist,r,4999,CFG['seed'])
    for j,key in enumerate(methods):
     for k,metric in enumerate(['high exposure','any high','restricted distance to high um']):
      if r<50 and k==2:continue
      a=vals[:,j,k];nu=a[1:];mn=nu.mean();sd=nu.std();valid=len(tgt)>1 and sd>1e-12
      out.append(dict(sample=s,platform='HD' if hd else 'Xenium',unit_resolution=res,identity_definition=name,state_definition=key,radius_or_contact_um=r,endpoint='reverse '+metric,observed=a[0],null_mean=mn,null_sd=sd,effect=a[0]-mn,OE=a[0]/mn if mn>0 else np.nan,P=(1+(abs(nu-mn)>=abs(a[0]-mn)-1e-12).sum())/5000 if valid else np.nan,targets=len(tgt),valid_centers=(np.diff(ptr)>0).sum(),coverage_pass=len(t),reason='estimable' if valid else 'target or null insufficient',family='reverse',block_um=400,version=CFG['version'],human_review_status=CFG['human_review_status']))
    if name=='supported' and r==50:
     tden=G@i.reliable_T.to_numpy('i4');area=e.area_50.values;valid=(e.coverage_50.values>=.5)&(np.diff(G.indptr)>0);g=G[:,tgt]
     for k,key in enumerate(methods):
      for end,den in [('density',area/1e6),('Treg/T',tden)]:
       ok=valid&(den>0);h=H[:,k];w=np.zeros(len(e));hh=ok&h;ll=ok&~h
       if not hh.any() or not ll.any():continue
       w[hh]=1/hh.sum()/den[hh];w[ll]=-1/ll.sum()/den[ll];con=g.T@w;obs=con.sum();loo.extend((s,res,key,end,z.index[idx],obs,obs-val) for idx,val in zip(tgt,con))
 put(pd.DataFrame(out),d/'reverse_effects.tsv.gz');put(pd.DataFrame(loo,columns=['sample','res','state','endpoint','excluded_id','full','leave_one']),d/'leave_one.tsv.gz')
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:run(s)
