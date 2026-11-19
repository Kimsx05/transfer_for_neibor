import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','RAYON_NUM_THREADS','ARROW_NUM_THREADS']:os.environ[k]='1'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:16])
from pathlib import Path
import json,sys
import numpy as np,pandas as pd,h5py,joblib
from scipy import sparse,ndimage
from scipy.spatial import cKDTree
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
R=Path(__file__).resolve().parents[1];A=R.parent;U=A/'hd_treg_unit_cell_pilot_v4';V=A/'hd_treg_mp10_spatial_all_v1';S=A/'segmentation_assignment_v2';N=A/'segmentation_restart_v1';O=A/'hd_mp10_treg_endo_precommunication_v1';Z=A/'hd_treg_mp10_enrichment_v3'
sections=['Cervilla','D1_Johnson']
T='CD3D CD3E CD3G TRAC CD2'.split();core='FOXP3 IL2RA CTLA4 IKZF2'.split();reg=core+'RTKN2 CCR8 TNFRSF18 TNFRSF4 ICOS LAYN MAGEH1'.split()
epi='EPCAM KRT7 KRT19 UPK2 KRT5 KRT14 KRT17 KRT6A KRT6B COL17A1'.split();fib='COL1A1 COL1A2 DCN LUM'.split();mye='LST1 LYZ TYROBP'.split()
top=(A.parents[1]/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split()
def rd(p):return pd.read_csv(p,sep='\t')
def put(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
def js(p):return json.loads(p.read_text())
def qc(s):
 d=pd.read_parquet(S/s/'object_qc.parquet',columns=['cell_id','original_nucleus_id','nucleus_x_um','nucleus_y_um','nuclear_area_um2','image_confirmed_invalid','geometry_review']).reset_index().rename(columns={'original_nucleus_id':'id','nucleus_x_um':'x','nucleus_y_um':'y'}).sort_values('id').reset_index(drop=True)
 assert np.array_equal(d.id,np.arange(1,len(d)+1));return d
def inputdata(s):
 audit=js(N/s/'metadata/input_audit.json')
 cols=['source_bin_index','array_row','array_col','pxl_col_in_fullres','pxl_row_in_fullres','raw_counts','nucleus_label','cell_label','in_tissue','status']
 d=pd.read_parquet(S/s/'bin_to_cell.parquet',columns=cols).sort_values('source_bin_index').reset_index(drop=True)
 assert np.array_equal(d.source_bin_index,np.arange(len(d)))
 d['x']=d.pxl_col_in_fullres*audit['native_mpp'];d['y']=d.pxl_row_in_fullres*audit['native_mpp'];d['usable']=d.in_tissue.eq(1)&~d.status.str.startswith('excluded')
 return d,audit
def genedata(audit):
 with h5py.File(audit['input_matrix'])as f:
  g=f['matrix/features'];return g['id'][:].astype(str),g['name'][:].astype(str)
def slabs(audit,chunk=100000):
 with h5py.File(audit['input_matrix'])as f:
  g=f['matrix'];ptr=g['indptr'][:];ng=len(g['features/id']);n=len(ptr)-1
  for a in range(0,n,chunk):
   b=min(n,a+chunk);lo,hi=ptr[a],ptr[b];yield a,b,sparse.csr_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr[a:b+1]-lo),shape=(b-a,ng))
def classify(x,sy):
 m=joblib.load(U/'model.joblib');ref=pd.Index(m['normalization_symbols']);idx=ref.get_indexer(sy);i=np.flatnonzero(idx>=0);mapping=sparse.csr_matrix((np.ones(len(i)),(i,idx[i])),shape=(len(sy),len(ref)))
 assert not set(m['features'])&set(top+'TIGIT NECTIN2 CD96 NECTIN1 PVR CD226'.split())
 feat=ref.get_indexer(m['features']);model=m['model'];records=[]
 for a in range(0,x.shape[0],4096):
  z=(x[a:a+4096]@mapping).tocsr();lib=np.asarray(z.sum(axis=1)).ravel();z=z[:,feat].multiply(10000/np.maximum(1,lib)[:,None]).tocsr();z.data=np.log1p(z.data)
  p=model.predict_proba(z);order=np.argsort(p,axis=1);n=np.arange(len(lib));it=list(model.classes_).index('Treg')
  rec=pd.DataFrame(dict(pred=model.classes_[order[:,-1]],runner=model.classes_[order[:,-2]],support=p[n,order[:,-1]],margin=p[n,order[:,-1]]-p[n,order[:,-2]],pTreg=p[:,it]))
  rec.loc[lib==0,['pred','runner']]=None;rec.loc[lib==0,['support','margin','pTreg']]=np.nan;records.append(rec)
 return pd.concat(records,ignore_index=True)
def ist(l):return l.fillna('').str.match('^(CD4|CD8|Treg|Tfh)')
