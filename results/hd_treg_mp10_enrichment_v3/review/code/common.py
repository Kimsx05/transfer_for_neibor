import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','OMP_THREAD_LIMIT']:os.environ[k]='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
from pathlib import Path
import numpy as np,pandas as pd,yaml,json,sys
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from datetime import datetime
from zoneinfo import ZoneInfo
R=Path(__file__).resolve().parents[1];A=R.parent;V=A/'hd_treg_mp10_spatial_all_v1';F=A/'hd_treg_mp10_spatial_followup_v2';H=A.parent
SECTIONS=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];LABELS=['FOXP3_supported','strict'];HIGHS=['H05','H10'];RADII=[50,200]
def put(d,f):d.to_csv(R/'tables'/f,sep='\t',index=False,na_rep='NA')
def rankz(v):
 r=rankdata(v);return (r-r.mean())/max(r.std(),1e-12)
def status(x):return 'NA'if np.isnan(x)else'Inf'if np.isinf(x)else'zero'if x==0 else'finite'
def divide(a,b):
 a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float));o=np.full(a.shape,np.nan);np.divide(a,b,out=o,where=b!=0);o[(b==0)&(a>0)]=np.inf;return o
def effects(a,b,c,d):
 a,b,c,d=[np.asarray(x,float)for x in [a,b,c,d]];N=a+b+c+d;pH=divide(a,a+c);pC=divide(b,b+d);base=divide(a+b,N)
 return {'FE':divide(pH,base),'PR':divide(pH,pC),'OR':divide(a*d,b*c)}
def empirical(obs,null,valid=True):
 null=np.asarray(null,float);sd=np.std(null,ddof=1);mu=np.mean(null)
 if not valid or not np.isfinite(obs)or not np.isfinite(sd)or sd<1e-12:return np.nan,np.nan,mu,sd
 return (1+np.sum(np.abs(null-mu)>=abs(obs-mu)-1e-12))/(len(null)+1),(obs-mu)/sd,mu,sd
def strata(d,degree):
 # No label, MP10 score, or high exposure used to choose cutpoints.
 q=d.copy();q['opportunity']=degree;bits=[]
 for col in ['nCount','T_UMI','opportunity']:
  med=q.groupby('component')[col].transform('median');bits.append((q[col]>med).astype(int))
 bits.append((q.coverage>=.8).astype(int));key=pd.MultiIndex.from_arrays([q.component]+bits);return pd.factorize(key)[0]
def perm_batches(labels,group,B,seed=123,batch=128):
 rng=np.random.default_rng(seed);indices=[np.flatnonzero(group==g)for g in np.unique(group)];labels=np.asarray(labels,np.float32)
 for start in range(0,B,batch):
  n=min(batch,B-start);out=np.broadcast_to(labels,(n,)+labels.shape).copy()
  for ix in indices:
   if len(ix)>1:
    order=np.argsort(rng.random((n,len(ix))),axis=1);out[:,ix,:]=labels[ix[order],:]
  yield start,out

def block_ids(d,size):return pd.factorize(pd.MultiIndex.from_arrays([d.component,np.floor(d.x_um/size),np.floor(d.y_um/size)]))[0]
def pairlist(d,r):
 out=[]
 for _,v in d.groupby('component'):
  pp=cKDTree(v[['x_um','y_um']]).query_pairs(r,output_type='ndarray')
  if len(pp):out.append(v.index.to_numpy()[pp])
 return np.concatenate(out)if out else np.empty((0,2),int)
