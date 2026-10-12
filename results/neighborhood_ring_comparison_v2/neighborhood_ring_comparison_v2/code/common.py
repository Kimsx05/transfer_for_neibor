import os,sys,json,hashlib,time
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism')
OUT=ROOT/'neighborhood_ring_comparison_v2'
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(ROOT/'spatial_ecology_top20_v2/runtime'),str(ROOT/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
sys.path.append(str(ROOT/'neighborhood_v1/.venv/lib/python3.12/site-packages'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(OUT/'logs/mplconfig')
import numpy as np,pandas as pd
from scipy import sparse

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n')
def tab(p,x):pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def log(x):print(time.strftime('%F %T'),x,flush=True)
def allocate(cap,total):
 cap=np.asarray(cap,int);a=np.zeros(len(cap),int)
 while a.sum()<total:
  active=np.flatnonzero(a<cap)
  if not len(active):break
  add=max(1,(total-int(a.sum()))//len(active))
  for j in active:
   d=min(add,int(cap[j]-a[j]),total-int(a.sum()));a[j]+=d
 return a

def sample_hierarchy(meta,n,seed):
 rng=np.random.default_rng(seed);indices=[];rows=[]
 def rec(df,level,budget):
  if level==3:
   ix=rng.choice(df.index.to_numpy(),budget,replace=False);indices.extend(ix.tolist())
   rows.append(dict(dataset_id=df.dataset_id.iloc[0],patient_stratum=df.patient_stratum.iloc[0],section_id=df.section_id.iloc[0],capacity=len(df),sampled=budget,inclusion_probability=budget/len(df)));return
  key=['dataset_id','patient_stratum','section_id'][level];groups=list(df.groupby(key,sort=True,dropna=False));alloc=allocate([len(v) for _,v in groups],budget)
  for (_,g),b in zip(groups,alloc):rec(g,level+1,int(b))
 rec(meta,0,min(n,len(meta)));return np.sort(indices),rows
sys.path.insert(0,str(OUT/'runtime'))
