from common import *
from scipy.linalg import qr

def fit_effect(d,predictor,ycol='Epi_MP10',adjusted=True,block_um=800,resamples=499,seed=123,covars=None):
 cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());res={'n':len(d),'n_blocks':0,'effect':np.nan,'low':np.nan,'high':np.nan,'P':np.nan,'status':'insufficient','block_um':block_um,'covariates':''}
 if len(d)<5:return res
 covars= covars if covars is not None else ['T_lineage','T_nCount','Epi_nCount','Epi_area_um2','observed_circle_fraction','T_background_balance']
 covars=[c for c in covars if c!=predictor and c in d]
 req=[predictor,ycol,'x_um','y_um','component']+(covars if adjusted else [])
 v=d.dropna(subset=req).copy();res['n']=len(v)
 if len(v)<5 or v[predictor].nunique()<3 or v[ycol].nunique()<3:return res
 def scale(x):
  z=rankdata(x);return (z-z.mean())/max(z.std(),1e-12)
 x=scale(v[predictor]);y=scale(v[ycol]);C=[np.ones(len(v))];kept=[]
 if adjusted:
  for c in covars:
   if v[c].nunique()<3:continue
   z=scale(v[c]);a=np.column_stack(C);rem=z-a@np.linalg.lstsq(a,z,rcond=None)[0]
   if np.var(rem)>.1:C.append(z);kept.append(c) # fixed-order VIF<10 against prior nuisance columns
 Z=np.column_stack(C);X=np.column_stack([Z,x]);beta=np.linalg.lstsq(X,y,rcond=None)[0][-1];res.update(effect=float(beta),covariates=','.join(kept))
 blocks=pd.factorize(pd.MultiIndex.from_arrays([v.component,np.floor(v.x_um/block_um),np.floor(v.y_um/block_um)]))[0];nb=blocks.max()+1;res['n_blocks']=int(nb)
 if len(v)<cfg['association_min_regions'] or nb<cfg['association_min_blocks']:res['status']='descriptive_insufficient_independent_support';return res
 if np.var(x-Z@np.linalg.lstsq(Z,x,rcond=None)[0])<.1:res['status']='descriptive_predictor_collinearity';return res
 rng=np.random.default_rng(seed);B=int(resamples);bs=[];idx=[np.flatnonzero(blocks==j)for j in range(nb)]
 for _ in range(B):
  ii=np.concatenate([idx[j]for j in rng.integers(0,nb,nb)]);bs.append(np.linalg.lstsq(X[ii],y[ii],rcond=None)[0][-1])
 # Restricted-null residuals: one Rademacher sign per spatial block; statistic uses fixed full X.
 fitted0=Z@np.linalg.lstsq(Z,y,rcond=None)[0];resid=y-fitted0;hat=np.sum(Z*np.linalg.pinv(Z).T,axis=1);resid=resid/np.sqrt(np.maximum(1-hat,1e-6));weight=np.linalg.pinv(X)[-1];group=np.bincount(blocks,weights=weight*resid,minlength=nb);null=rng.choice([-1.,1.],size=(B,nb))@group
 p=(1+np.sum(np.abs(null)>=abs(beta)))/(B+1);lo,hi=np.quantile(bs,[.025,.975]);res.update(low=float(lo),high=float(hi),P=float(p),status='within_section_exploratory_test');return res

def bh(d,groupcols):
 from statsmodels.stats.multitest import multipletests
 d=d.copy();d['BH_q']=np.nan
 for _,ix in d.groupby(groupcols,dropna=False).groups.items():
  ii=[i for i in ix if np.isfinite(d.loc[i,'P'])]
  if ii:d.loc[ii,'BH_q']=multipletests(d.loc[ii,'P'],method='fdr_bh')[1]
 return d
