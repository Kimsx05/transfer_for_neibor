from pathlib import Path
import json,sys,hashlib,os
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
O=Path(__file__).resolve().parents[1];ROOT=O.parents[1];V=O.parent/'mp10_visium_niche_boundary_v1'
from frozen_geometry import tissue_graph,band,CATS,sha,geometry
SEED=20261127;AREA=np.sqrt(3)*.005
B=['all_E','E_without_I_enriched'];G=['infiltrated-like','non-infiltrated-like','mixed']
A=['E','I','T','CD8','Treg','myeloid'];SCORE='EpiMP10_Top20_AUCell'
def read(p):return pd.read_csv(p,sep='\t')
M=read((O if (O/'audit/input_manifest.tsv').exists() else V)/'audit/input_manifest.tsv').set_index('section_id',drop=False)
def save(d,p):
    path=O/p;path.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(d).to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.12g')
def pq(d,p):
    path=O/p;path.parent.mkdir(parents=True,exist_ok=True);d.to_parquet(path,index=False)
def load(p):return pd.read_parquet(O/p,use_threads=False)
def bh(p,total=None):
    v=np.asarray(p,float);ok=np.isfinite(v);out=np.full(len(v),np.nan)
    if ok.any():out[ok]=np.minimum(1,stats.false_discovery_control(v[ok])*(total or ok.sum())/ok.sum())
    return out
def safe(a,b):return a/b if b>0 else np.nan
def avg(x):return np.mean(x) if len(x) else np.nan
def unit_summary(x):
    x=pd.Series(x).dropna();n=len(x);mu=x.mean();se=x.std()/np.sqrt(n) if n>1 else np.nan;ci=stats.t.ppf(.975,n-1)*se if n>=2 else np.nan
    return dict(effect=mu,low=mu-ci,high=mu+ci,p=2*stats.t.sf(abs(mu/se),n-1) if n>=3 and se>0 else np.nan)
def longmeta(sec):return dict(section_id=sec,dataset_id=M.at[sec,'dataset_id'],analysis_unit_uid=M.at[sec,'analysis_unit_uid'],patient_uid=M.at[sec,'patient_uid'])

