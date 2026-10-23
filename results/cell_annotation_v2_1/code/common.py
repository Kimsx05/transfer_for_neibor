from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1];V2=ROOT.parent/'cell_annotation_v2';V1=ROOT.parent/'cell_annotation_freeze_v1';SEED=123
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8',PYTHONHASHSEED='123')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
def js(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda z:z.item() if hasattr(z,'item') else str(z)))
def read_section(name):
 import anndata as ad
 return ad.read_h5ad(V2/name/'cells.h5ad')
