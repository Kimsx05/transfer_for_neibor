import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
from pathlib import Path
import numpy as np,pandas as pd,json,yaml
import warnings
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
from scipy import sparse
import h5py
R=Path(__file__).resolve().parents[1];H=R.parents[1];P=H.parent;A=H/'analysis'
PANELS={'Epi':['EPCAM','KRT7','KRT19','UPK2','KRT5','KRT14','KRT6A','KRT6B','KRT17','COL17A1'],'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','CD2'],'Treg':['FOXP3','IL2RA','CTLA4','IKZF2','LRRC32'],'Tfh':['CXCL13','PDCD1','ICOS','BCL6','MAF'],'CD4':['CD4','IL7R','CCR7','TCF7','LEF1'],'CD8':['CD8A','CD8B','CCL5','GZMK'],'NK':['NKG7','GNLY','KLRD1'],'Endo':['PECAM1','VWF','CDH5','KDR','EMCN','RAMP2','ESAM'],'Fib':['COL1A1','COL1A2','DCN','LUM'],'Myeloid':['LYZ','LST1','TYROBP','FCER1G'],'B':['MS4A1','CD79A','CD79B'],'Plasma':['MZB1','JCHAIN','SDC1'],'Mural':['RGS5','ACTA2','MYH11'],'Mast':['TPSAB1','TPSB2','KIT']}
TOP=(P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split()
def put(d,n):d.to_csv(R/'tables'/n,sep='\t',index=False,na_rep='NA')
def read8(s):
 info=json.loads((A/'hd_bin8_comparison_v1'/s/'input_summary.json').read_text())
 path=H/info['matrix_path']
 with h5py.File(path) as f:
  g=f['matrix'];x=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:])).T.tocsr();ids=g['barcodes'][:].astype(str);genes=g['features/id'][:].astype(str);symbols=g['features/name'][:].astype(str)
 d=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet').set_index('bin_id').loc[ids].reset_index();return x,d,genes,symbols,path

def panel(x,d,sy):
 for k,v in PANELS.items():
  z=x[:,np.isin(sy,v)];d[k+'_detected']=np.asarray((z>0).sum(1)).ravel();d[k+'_UMI']=np.asarray(z.sum(1)).ravel()
 for g in sorted(set(sum(PANELS.values(),[])+TOP+['TIGIT'])):
  d['raw_'+g]=np.asarray(x[:,sy==g].sum(1)).ravel()
 d['MP10_detected']=np.asarray((x[:,np.isin(sy,TOP)]>0).sum(1)).ravel();return d

def gate():
 from datetime import datetime
 from zoneinfo import ZoneInfo
 import time
 while 7<=datetime.now(ZoneInfo('Asia/Shanghai')).hour<23:time.sleep(30)
