"""Sparse two-sided Wilcoxon rank sums. All zeros share one rank; positive entries are ranked exactly."""
from data import *
from scipy.special import log_ndtr
from scipy.stats import mannwhitneyu
from contextlib import ExitStack
import warnings
p=sys.argv[1];d=O/p;X,genes,rawq=combine(p);z=pd.read_parquet(d/'objects.parquet').set_index('cell_key')
assert z.index.is_unique and set(z.index)==set(rawq.cell_key);z=z.loc[rawq.cell_key]
n,g=X.shape;lib=rawq.n_counts.to_numpy(dtype=float);X=X.tocsc();X.sort_indices()
samples=S[p];sample_codes=pd.Categorical(rawq['sample'],categories=samples).codes.astype('i4');ns=len(samples)
sample_sizes=np.bincount(sample_codes,minlength=ns);grp={r:z[cname(r)].to_numpy('i4') for r in R}
# Calculate all per-sample normalized sums/detection counts and pooled rank sums in one sparse pass.
ks={r:int(grp[r].max())+1 for r in R}
by={r:sample_codes*ks[r]+grp[r] for r in R}
N={r:np.bincount(by[r],minlength=ns*ks[r]).reshape(ns,ks[r]) for r in R}
sums={r:np.zeros((ns,ks[r],g)) for r in R};det={r:np.zeros((ns,ks[r],g),np.int32) for r in R}
rank_excess={r:np.zeros((ks[r],g)) for r in R}
tievar=np.zeros(g);status=np.full(g,'ok',object);zero_rank=np.zeros(g)
for j in range(g):
 l,h=X.indptr[j:j+2];idx=X.indices[l:h];nz=len(idx);n0=n-nz
 if not nz:status[j]='all_zero';continue
 norm=X.data[l:h].astype('f8')*10000/lib[idx];v=np.log1p(norm)
 unique,inv,t=np.unique(v,return_inverse=True,return_counts=True);t=t.astype(float)
 zr=(n0+1)/2;zero_rank[j]=zr
 avgpos=n0+np.cumsum(t)-(t-1)/2
 delta=avgpos[inv]-zr
 tie_terms=(float(n0)**3-n0)+(t**3-t).sum();tievar[j]=((n+1)-tie_terms/(float(n)*(n-1)))/12
 if n0==0 and len(unique)==1:status[j]='constant'
 for r in R:
  k=ks[r];sums[r][:,:,j]=np.bincount(by[r][idx],weights=norm,minlength=ns*k).reshape(ns,k)
  det[r][:,:,j]=np.bincount(by[r][idx],minlength=ns*k).reshape(ns,k)
  rank_excess[r][:,j]=np.bincount(grp[r][idx],weights=delta,minlength=k)
def bh(lp):
 good=np.flatnonzero(np.isfinite(lp));adj=np.full(len(lp),np.nan)
 if len(good):
  order=good[np.argsort(lp[good],kind='stable')];v=lp[order]+np.log(len(order)/np.arange(1,len(order)+1));adj[order]=np.minimum(0,np.minimum.accumulate(v[::-1])[::-1])
 return adj
def effect(si,so,di,do,ni,no):
 ai=si/ni;ao=np.where(do==0,0.,np.maximum(so,0))/no
 pi=di/ni;po=do/no;return ai,ao,pi,po,pi-po,np.log2((ai+1e-9)/(ao+1e-9))
base=genes.copy();base['measured_samples']='|'.join(samples)
ver=[]
for r in R:
 dr=d/('r'+r);tops=[];bottoms=[];st=sums[r].sum(1);dt=det[r].sum(1)
 with ExitStack() as ctx:
  files=[ctx.enter_context(gzip.open(dr/(key+'.tsv.gz.part'),'wt',compresslevel=3)) for key in ['pooled_all_genes','within_sample_all_genes']]
  for c in range(ks[r]):
   ni=int(N[r][:,c].sum());no=n-ni;di=det[r][:,c].sum(0);do=dt.sum(0)-di;si=sums[r][:,c].sum(0);so=st.sum(0)-si
   av=effect(si,so,di,do,ni,no);frame=base.assign(platform=p,resolution=r,cluster=c,n_in=ni,n_out=no)
   fields=['mean_norm_in','mean_norm_out','pct_in','pct_out','delta_pct','log2FC']
   for col,val in zip(fields,av):frame[col]=val
   centered=rank_excess[r][c]+ni*(zero_rank-(n+1)/2);u=centered+float(ni)*no/2
   state=status.copy()
   if min(ni,no)<2:state[:] = 'insufficient_group'
   good=(state=='ok')&(tievar>0);state[(state=='ok')&(tievar<=0)]='constant'
   zz=np.full(g,np.nan);zz[good]=centered[good]/np.sqrt(float(ni)*no*tievar[good]);lp=np.minimum(0,np.log(2)+log_ndtr(-abs(zz)));la=bh(lp)
   frame['test_statistic']=np.where(good,u,np.nan);frame['z_tie_corrected']=zz;frame['p_value']=np.exp(lp);frame['p_adj_bh']=np.exp(la);frame['neg_log10_p']=-lp/np.log(10);frame['neg_log10_fdr']=-la/np.log(10)
   frame['direction']=np.select([av[-1]>0,av[-1]<0],['positive','negative'],default='zero');frame['status']=state
   effects=np.full((ns,g),np.nan)
   for h,s in enumerate(samples):
    a=int(N[r][h,c]);b=int(sample_sizes[h]-a)
    w=base.assign(platform=p,resolution=r,cluster=c,sample=s,n_in=a,n_out=b)
    if min(a,b)<20:
     for col in fields:w[col]=np.nan
     w['status']='insufficient_n_in' if a<20 else 'insufficient_n_out'
    else:
     ai=det[r][h,c];ao=dt[h]-ai;vs=effect(sums[r][h,c],st[h]-sums[r][h,c],ai,ao,a,b)
     for col,val in zip(fields,vs):w[col]=val
     effects[h]=vs[-1];w['status']=np.where(ai+ao==0,'measured_all_zero','ok')
    w.to_csv(files[1],sep='\t',index=False,header=c==0 and h==0,na_rep='NA',float_format='%.12g')
   evaluable=np.isfinite(effects);ne=evaluable.sum(0);agree=((np.sign(effects)==np.sign(av[-1]))&evaluable&(av[-1]!=0)).sum(0)
   frame['n_evaluable_samples']=ne;frame['n_same_direction']=np.where(av[-1]!=0,agree,np.nan);frame['fraction_same_direction']=np.divide(agree,np.maximum(ne,1),dtype=float);frame.loc[(ne==0)|(av[-1]==0),'fraction_same_direction']=np.nan
   with warnings.catch_warnings():
    warnings.simplefilter('ignore',RuntimeWarning)
    for suf,fn in [('min',np.nanmin),('median',np.nanmedian),('max',np.nanmax)]:frame['within_log2FC_'+suf]=fn(effects,axis=0)
   frame.to_csv(files[0],sep='\t',index=False,header=c==0,na_rep='NA',float_format='%.12g')
   pos=frame[good&(frame.p_adj_bh<=.05)&(frame.log2FC>=.25)&(frame.pct_in>=.1)];neg=frame[good&(frame.p_adj_bh<=.05)&(frame.log2FC<=-.25)&(frame.pct_out>=.1)]
   for f,target in [(pos,tops),(neg,bottoms)]:
    top=f.sort_values(['neg_log10_fdr','neg_log10_p','feature_id'],ascending=[False,False,True]).head(20).copy();top['rank']=np.arange(1,len(top)+1);target.append(top)
   if c==0:
    # Independent tests cover genes with low, intermediate and high detection plus an all-zero if present.
    valid=np.flatnonzero(good)
    for j in valid[np.argsort((di+do)[valid])][np.linspace(0,len(valid)-1,5).astype(int)]:
     v=np.log1p(X[:,j].toarray().ravel().astype(float)*10000/lib);sel=grp[r]==c
     test=mannwhitneyu(v[sel],v[~sel],alternative='two-sided',method='asymptotic',use_continuity=False)
     assert abs(test.statistic-u[j])<1e-4 and np.isclose(test.pvalue,np.exp(lp[j]),rtol=1e-7,atol=1e-300)
     ver.append((p,r,c,genes.feature_id.iloc[j],test.statistic,u[j],test.pvalue,np.exp(lp[j])))
 for k in ['pooled_all_genes','within_sample_all_genes']:(dr/(k+'.tsv.gz.part')).replace(dr/(k+'.tsv.gz'))
 write(pd.concat(tops),dr/'top20_positive.tsv');write(pd.concat(bottoms),dr/'top20_negative.tsv')
 print('full DEG',r,ks[r],g,flush=True)
write(pd.DataFrame(ver,columns=['platform','resolution','cluster','feature_id','scipy_U','our_U','scipy_P','our_P']),d/'wilcoxon_checks.tsv')

