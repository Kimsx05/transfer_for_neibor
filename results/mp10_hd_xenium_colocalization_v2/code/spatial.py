from common import *
from numba import njit
from scipy.stats import spearmanr
from scipy.spatial.distance import cdist
from scipy.special import expit
def ratio(x,n):return np.divide(x,n,out=np.full(np.shape(x),np.nan,dtype=float),where=np.asarray(n)>0)
def groups(y,scale,depth=True):
 a=y[['component','x_um','y_um']].copy();a[['x_um','y_um']]=np.floor(a[['x_um','y_um']]/scale)
 if depth:a['depth']=y.UMI.gt(y.UMI.median()).values
 return pd.factorize(pd.MultiIndex.from_frame(a))[0]
@njit
def shuffle_effect(H,Y,groups,nperm,seed):
 np.random.seed(seed);order=np.argsort(groups);offset=np.concatenate((np.zeros(1,np.int64),np.cumsum(np.bincount(groups))));N,K=H.shape;M=Y.shape[1]
 out=np.full((nperm+1,K,M),np.nan);sums=np.zeros(M);ns=np.zeros(M)
 for i in range(N):
  for j in range(M):
   if np.isfinite(Y[i,j]):sums[j]+=Y[i,j];ns[j]+=1
 idx=np.arange(N)
 for it in range(nperm+1):
  if it:
   for l in range(len(offset)-1):
    for u in range(offset[l+1]-1,offset[l],-1):
     a=order[u];b=order[np.random.randint(offset[l],u+1)];idx[a],idx[b]=idx[b],idx[a]
  for k in range(K):
   sh=np.zeros(M);nh=np.zeros(M)
   for i in range(N):
    if H[i,k]:
     v=Y[idx[i]]
     for j in range(M):
      if np.isfinite(v[j]):sh[j]+=v[j];nh[j]+=1
   for j in range(M):
    if nh[j]>0 and ns[j]>nh[j]:out[it,k,j]=sh[j]/nh[j]-(sums[j]-sh[j])/(ns[j]-nh[j])
 return out
def boot(y,h,blk):
 valid=np.isfinite(y);nb=max(blk)+1
 a=np.column_stack([np.bincount(blk[valid&(h==u)],weights=w[valid&(h==u)],minlength=nb) for u in [1,0] for w in [np.nan_to_num(y),np.ones(len(y))]])
 rng=np.random.default_rng(CFG['seed']);ci=[]
 for i in range(0,999,40):
  chunk=a[rng.integers(nb,size=(min(40,999-i),nb))].sum(1);v=(chunk[:,1]>0)&(chunk[:,3]>0)
  ci.extend(chunk[v,0]/chunk[v,1]-chunk[v,2]/chunk[v,3])
 return np.quantile(ci,[.025,.975]) if ci and nb>=2 else [np.nan,np.nan]
def spatial(s):
 hd=s in HD;v,d=folders(s);z=labels(s);idt=pd.read_parquet(d/'identity.parquet');xy=z[['x_um','y_um']].values
 cov=pd.read_parquet(d/'target_cover.parquet');methods=['MP10'] if hd else STATES
 result=[];reverse=[];allows=[];curves=[]
 for res in ([16,32] if hd else ['cell']):
  ee=pd.read_parquet(d/f'centers_{res}.parquet');Y={};spec=[];H=ee[[k+'_high' for k in methods]].to_numpy(bool);ex=ee[['x_um','y_um']].values
  for r in [50,200]:
   G=sparse.load_npz(d/f'epi_{res}_{r}.npz');n=np.diff(G.indptr);t=G@idt.reliable_T.to_numpy('i4');known=G@(~z.major.isin(['Mixed','Unknown'])).to_numpy('i4');valid=(ee[f'coverage_{r}'].values>=.5)&(n>0)
   support=pd.DataFrame({'id':ee.index,'neighbors':n,'known':known,'reliable_T':t,'coverage':ee[f'coverage_{r}'].values,'valid':valid})
   for name,grade in [('supported',['supported']),('broad',['supported','possible'])]:
    target=np.flatnonzero(idt.grade.isin(grade));num=G@idt.grade.isin(grade).to_numpy('i4');support[name]=num
    def add(k,arr,ok=valid,scale=r):
     ky=(name,scale,k);Y[ky]=np.where(ok,arr,np.nan);spec.append(ky)
    add('density',num*1e6/ee[f'area_{r}'].clip(lower=1));add('Treg/all',ratio(num,n));add('Treg/T',ratio(num,t))
    for cut in [3,5]:add('Treg/T>='+str(cut),ratio(num,t),valid&(t>=cut))
    if hd and name=='supported':
     w=idt.treg_weight;add('Treg weight intensity',ratio(G@w.fillna(0).to_numpy(),G@w.notna().to_numpy('i4')))
    if r==200:
     g=G[:,target].tocoo();ds=np.linalg.norm(ex[g.row]-xy[target[g.col]],axis=1);nearest=np.full(len(ee),np.inf);np.minimum.at(nearest,g.row,ds)
     add('nearest observed um',np.where(np.isfinite(nearest),nearest,np.nan));add('restricted distance um',np.minimum(nearest,200))
     support[name+'_distance_um']=np.where(np.isfinite(nearest),nearest,np.nan)
     for radius in [25,50,100,150,200]:add('CDF_'+str(radius),np.asarray(nearest<=radius,float))
    gg=G[:,target].T.tocsr();ts=cov.loc[z.index[target],f'coverage_{r}'].to_numpy()>=.5
    gg=gg[ts];ii=target[ts];den=np.diff(gg.indptr)
    # Reverse statistic uses all eligible Epi neighbors; area gate applies to Treg centers.
    for k,key in enumerate(methods):
     high=gg@H[:,k].astype(int);rf=ratio(high,den);rev=pd.DataFrame({'id':z.index[ii],'epi':den,'high':high,'fraction':rf})
     if r==200:
      co=gg.tocoo();dst=np.linalg.norm(xy[ii[co.row]]-ex[co.col],axis=1);best=np.full(len(ii),np.inf);keep=H[co.col,k];np.minimum.at(best,co.row[keep],dst[keep]);rev['high_distance_um']=np.where(np.isfinite(best),best,np.nan)
     rev.to_parquet(d/f'reverse_{res}_{r}_{name}_{key}.parquet',index=False)
     reverse.append((s,res,r,name,key,len(target),len(ii),sum(den>0),np.nanmean(rf) if (den>0).any() else np.nan,(high[den>0]>0).mean() if (den>0).any() else np.nan))
    if not hd and r==50:
     pairs=pd.read_parquet(d/'boundary_pairs.parquet')
     for gap in [1,2,5]:
      p=pairs[pairs.gap_um<=gap];cg=sparse.csr_matrix((np.ones(len(p),np.int32),(p.row,p.neighbor)),shape=G.shape)
      ct=cg@idt.grade.isin(grade).to_numpy('i4');add('contact edges',ct,valid,gap);add('any contact',1.*(ct>0),valid,gap)
      cgt=cg[:,target].T;den2=np.asarray(cgt.sum(1)).ravel()
      for k,key in enumerate(methods):
       hh=cgt@H[:,k].astype(int);qq=ts&(den2>0);reverse.append((s,res,'contact'+str(gap),name,key,len(target),ts.sum(),qq.sum(),np.mean(hh[qq]/den2[qq]) if qq.any() else np.nan,np.mean(hh[qq]>0) if qq.any() else np.nan))
   support.to_parquet(d/f'neighborhood_{res}_{r}.parquet',index=False)
  vals=np.ascontiguousarray(np.column_stack([Y[a] for a in spec]));blk=groups(ee,400,False);qc=pd.read_csv(d/'score_QC.tsv',sep='\t',dtype={'res':str})
  for scale in [400,800]:
   grp=groups(ee,scale);P=shuffle_effect(np.ascontiguousarray(H),vals,grp,4999,CFG['seed'])
   assert np.allclose(P[0],np.array([[np.nanmean(vals[h],axis=0)-np.nanmean(vals[~h],axis=0) for h in H.T]])[0],equal_nan=True,atol=1e-8)
   for k,key in enumerate(methods):
    h=H[:,k];swap=pd.DataFrame({'g':grp,'h':h}).groupby('g').h.agg(['sum','size']);exchg=swap.loc[(swap['sum']>0)&(swap['sum']<swap['size']),'size'].sum()
    allows.append((s,res,key,scale,len(swap),exchg,len(h)))
    for j,(name,r,end) in enumerate(spec):
     y=vals[:,j];obs=P[0,k,j];nu=P[1:,k,j];nu=nu[np.isfinite(nu)];targ=idt.grade.isin(['supported'] if name=='supported' else ['supported','possible']).sum()
     valid=np.isfinite(y);vh=valid&h;vl=valid&~h;mean=nu.mean() if len(nu) else np.nan;sd=nu.std() if len(nu) else np.nan
     reason='estimable' if targ>=2 and vh.any() and vl.any() and len(nu)==4999 and sd>1e-10 and exchg>0 else 'sparse target / absent denominator / no exchange or degenerate null'
     p=(1+np.sum(abs(nu-mean)>=abs(obs-mean)-1e-12))/5000 if reason=='estimable' else np.nan
     ci=boot(y,h,blk) if scale==400 and vh.any() and vl.any() else [np.nan]*2
     # Direct standardization to pooled depth/area strata represented in both arms.
     adj=pd.DataFrame({'y':y,'h':h,'depth':ee.UMI.gt(ee.UMI.median()).values,'area':(ee.cell_area>ee.cell_area.median()).values if 'cell_area' in ee else ee.n_bins.to_numpy()}).dropna().groupby(['depth','area','h']).y.agg(['mean','count'])
     tab=adj['mean'].unstack();weight=adj['count'].groupby(level=[0,1]).sum();diff=tab.get(True,np.nan)-tab.get(False,np.nan);balanced=np.average(diff.dropna(),weights=weight.loc[diff.dropna().index]) if isinstance(diff,pd.Series) and diff.notna().any() else np.nan
     family='primary' if scale==400 and name=='supported' and r==50 and end in ['Treg/T','density'] and res!=32 else 'secondary'
     result.append(dict(platform='HD' if hd else 'Xenium',sample=s,state_definition=key,unit_resolution=res,identity_definition=name,radius_or_contact_um=r,endpoint=end,block_um=scale,effect=obs,high_mean=np.nanmean(y[vh]) if vh.any() else np.nan,low_mean=np.nanmean(y[vl]) if vl.any() else np.nan,balanced_effect=balanced,CI_low=ci[0],CI_high=ci[1],null_mean=mean,null_sd=sd,relative_deviation=obs-mean,P=p,family=family,targets=targ,valid_high=vh.sum(),valid_low=vl.sum(),Spearman=ee[key].corr(pd.Series(y,index=ee.index),method='spearman'),reason=reason,score_QC=qc.loc[(qc.res==str(res))&(qc.state==key),'QC'].item(),version=CFG['version'],human_review_status=CFG['human_review_status']))
    for name in ['supported','broad']:
     ix=[j for j,p0 in enumerate(spec) if p0[0]==name and p0[2].startswith('CDF_')];vs=P[:,k,ix];mu=vs[1:].mean(0);sd=vs[1:].std(0);ok=sd>1e-10
     if ok.any():
      mx=np.max(np.abs((vs[:,ok]-mu[ok])/sd[ok]),axis=1);crit=np.quantile(mx[1:],.95)
      for jj,j in enumerate(ix):curves.append((s,res,name,key,scale,spec[j][2],vs[0,jj],mu[jj]-crit*sd[jj],mu[jj]+crit*sd[jj],(1+(mx[1:]>=mx[0]).sum())/5000))
  print(s,res,'spatial effects complete',flush=True)
 put(pd.DataFrame(result),d/'effects.tsv.gz');put(pd.DataFrame(reverse,columns=['sample','res','scale','definition','state','targets','coverage_pass','with_epi','high_exposure','any_high']),d/'reverse_summary.tsv')
 put(pd.DataFrame(curves,columns=['sample','res','definition','state','block','distance','observed','global_lower','global_upper','global_P']),d/'distance_envelope.tsv')
 put(pd.DataFrame(allows,columns=['sample','res','state','block','strata','exchangeable','eligible']),d/'exchange.tsv')
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:spatial(s)
