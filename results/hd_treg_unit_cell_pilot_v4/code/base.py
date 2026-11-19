import os
for name in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS','ARROW_NUM_THREADS']: os.environ[name]='1'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:4])
from pathlib import Path
import json, yaml, numpy as np, pandas as pd, h5py
from scipy import sparse
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
R=Path(__file__).resolve().parents[1]; A=R.parent; V=A/'hd_treg_mp10_spatial_all_v1'; Z=A/'hd_treg_mp10_enrichment_v3'; S=A/'segmentation_assignment_v2'; N=A/'segmentation_restart_v1'; O=A/'hd_mp10_treg_endo_precommunication_v1'; F=A/'hd_treg_mp10_spatial_followup_v2'
sections=['Cervilla','D1_Johnson']
T='CD3D CD3E CD3G TRAC CD2'.split(); core='FOXP3 IL2RA CTLA4 IKZF2'.split(); epi='EPCAM KRT7 KRT19 UPK2 KRT5 KRT14 KRT17 KRT6A KRT6B COL17A1'.split(); fib='COL1A1 COL1A2 DCN LUM'.split(); mye='LST1 LYZ TYROBP'.split()
exclude='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1 TIGIT NECTIN2 CD96 NECTIN1 PVR CD226'.split()
def rd(p):return pd.read_csv(p,sep='\t')
def put(d,p):d.to_csv(R/p,sep='\t',index=False,na_rep='NA')
def js(p):return json.loads(Path(p).read_text())
def roi_for(x,y,roi):
 out=np.full(len(x),'',dtype=object)
 for q in roi.itertuples():
  on=(x>=q.x0)&(x<q.x1)&(y>=q.y0)&(y<q.y1);out[on]=[f'{v}|{q.roi}'.strip('|') for v in out[on]]
 return out
