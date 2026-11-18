import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','OMP_THREAD_LIMIT']:os.environ[k]='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
from pathlib import Path
import numpy as np,pandas as pd,json,yaml,sys,importlib.util
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from datetime import datetime
from zoneinfo import ZoneInfo
R=Path(__file__).resolve().parents[1];A=R.parent;V=A/'hd_treg_mp10_spatial_all_v1';H=A.parent;P=H.parent;OLD=A/'hd_mp10_treg_endo_precommunication_v1'
sys.path.insert(0,str(V/'code'))
import common as v1
SECTIONS=v1.SECTIONS;cfg=yaml.safe_load((V/'frozen_config.yaml').read_text());SIG=cfg['signatures']['Treg_main'];CORE=['FOXP3','IL2RA','CTLA4','IKZF2'];TL=cfg['signatures']['T_lineage']
def put(d,name):d.to_csv(R/'tables'/name,sep='\t',index=False,na_rep='NA')
def stamp():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def zrank(x):
 z=rankdata(x);return (z-z.mean())/max(z.std(),1e-12)
def wilson(k,n):
 if n==0:return (np.nan,np.nan)
 p=k/n;z=1.96;den=1+z*z/n;cen=(p+z*z/(2*n))/den;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;return max(0,cen-h),min(1,cen+h)
