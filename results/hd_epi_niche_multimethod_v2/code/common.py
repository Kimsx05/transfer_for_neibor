from pathlib import Path
import os,json,time
os.environ['OMP_NUM_THREADS']='4';os.environ['NUMBA_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['TF_NUM_INTRAOP_THREADS']='4';os.environ['TF_NUM_INTEROP_THREADS']='1'
import numpy as np,pandas as pd,yaml
import pyarrow as pa
pa.set_cpu_count(4);pa.set_io_thread_count(4)
R=Path(__file__).resolve().parents[1];V=R.parent/'hd_epi_niche_dual_methods_v1';E=R.parent/'hd_bin8_banksy_extension_v1';HD=R.parents[1]
CFG=yaml.safe_load((R/'config.yaml').read_text());SECTIONS=CFG['sections'];TYPES=json.loads((R/'type_mapping.json').read_text())
def anchors():return pd.read_parquet(R/'tables/frozen_anchors.parquet')
def train_index(a):
 t=pd.read_csv(R/'tables/T_v2.tsv',sep='\t');ids=pd.Index(a.section+'::'+a.bin_id).get_indexer(t.section+'::'+t.bin_id);assert (ids>=0).all();return ids

def composition(ref,r):
 d=pd.concat([pd.read_parquet(V/s/f'r{r}.parquet') for s in SECTIONS],ignore_index=True);a=anchors();assert np.array_equal(d.section+'::'+d.bin_id,a.key);return d,d[[ref+'_P_'+x for x in TYPES[ref]]].to_numpy(),(d[ref+'_n_known']>=10).to_numpy()
def gate():
 from datetime import datetime
 from zoneinfo import ZoneInfo
 while 7<=datetime.now(ZoneInfo('Asia/Shanghai')).hour<23:
  print('NIGHT_WINDOW_WAIT',datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),flush=True);time.sleep(30)
 os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[int(os.environ.get('NICHE_CORE_OFFSET','0')):int(os.environ.get('NICHE_CORE_OFFSET','0'))+4]))
def source(s):
 base=R.parent/'hd_bin8_comparison_v1' if s in SECTIONS[:2] else E
 info=json.loads((base/s/'input_summary.json').read_text());path=HD/info['matrix_path'] if s in SECTIONS[:2] else (base/s/info['counts_relative']).resolve()
 return base/s,path

def counts(s,genes=None,ids=None):
 import h5py
 from scipy.sparse import csc_matrix
 _,p=source(s)
 with h5py.File(p) as f:
  g=f['matrix'];b=g['barcodes'][:].astype(str);gn=g['features/id'][:].astype(str);sy=g['features/name'][:].astype(str)
  # sparse genes x supplied bins; discard original view after selected rows/columns.
  x=csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]))
 if ids is None:
  n=pd.read_parquet(V/s/'nodes.parquet',columns=['bin_id','analysis_qc']);ids=n.loc[n.analysis_qc,'bin_id'].to_numpy()
 ix=pd.Index(b).get_indexer(ids);assert (ix>=0).all();x=x[:,ix]
 if genes is not None:
  gi=pd.Index(gn).get_indexer(genes);assert (gi>=0).all();x=x[gi];gn=gn[gi];sy=sy[gi]
 return x.T.tocsr(),np.asarray(ids),gn,sy

def buffer_split(a,T,r):
 from scipy.spatial import cKDTree
 d=a.iloc[T];h=((d.x_um//400).astype(int)*73856093+(d.y_um//400).astype(int)*19349663)%7
 test=(h==0).to_numpy();keep=~test;nearest=np.full(len(T),np.inf)
 for s in SECTIONS[:4]:
  q=np.flatnonzero((d.section==s).to_numpy());tt=q[test[q]];tr=q[~test[q]]
  if len(tt) and len(tr):nearest[tr]=cKDTree(d.iloc[tt][['x_um','y_um']]).query(d.iloc[tr][['x_um','y_um']],workers=4)[0];keep[tr]=nearest[tr]>2*r+1e-7
 return keep,test,nearest
