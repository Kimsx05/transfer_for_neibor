from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1];HD=ROOT.parents[1];V1=ROOT.parent/'cell_annotation_freeze_v1';SEED=123
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8',TF_NUM_INTRAOP_THREADS='8',TF_NUM_INTEROP_THREADS='2',PYTHONHASHSEED='123')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
def section(name):
 key='Cervilla' if name=='Cervilla' else 'D1';meta=json.loads((HD/'input_recovery_followup'/f'{key}_recovered_image.json').read_text());base=HD/('data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um' if name=='Cervilla' else 'data/GSE300692/D1_Johnson/outs/binned_outputs/square_002um');out=ROOT/name
 for d in ['metadata','figures','inputs']: (out/d).mkdir(exist_ok=True,parents=True)
 return meta,base,out

def setup_gpu_env():
 import sys
 overlay=ROOT/'.gpu_packages';sys.path.insert(0,str(overlay));os.environ['KERAS_HOME']=str(ROOT/'models/keras_cache');os.environ['CELLTYPIST_FOLDER']=str(ROOT/'models/celltypist_cache');os.environ['CUDA_VISIBLE_DEVICES']='0'

def jwrite(path,x):path.write_text(json.dumps(x,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)))
