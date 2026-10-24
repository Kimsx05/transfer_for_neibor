from pathlib import Path
import os,json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent;SEED=123
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
def js(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda a:a.item() if hasattr(a,'item') else str(a)))
def load():
 import numpy as np,pandas as pd
 o=pd.read_parquet(BASE/'joint_umap_major_review/working_annotations.parquet');z=np.load(BASE/'cell_annotation_v2_1/metadata/existing_coordinates.npz');idx=pd.Index(z['cell_id']).get_indexer(o.cell_id);assert (idx>=0).all() and o.cell_id.is_unique
 return o,z['X_scVI'][idx],z['spatial_um'][idx]
