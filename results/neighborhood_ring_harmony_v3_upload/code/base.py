import os,sys,json,hashlib,time
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism'); OLD=ROOT/'neighborhood_ring_comparison_v2'; OUT=ROOT/'neighborhood_ring_harmony_v3'
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(OLD/'runtime'),str(ROOT/'spatial_ecology_top20_v2/runtime'),str(ROOT/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
sys.path.append(str(ROOT/'neighborhood_v1/.venv/lib/python3.12/site-packages'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(OUT/'logs/mplconfig')
import numpy as np,pandas as pd
from scipy import sparse
import joblib

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n')
def tab(p,x):pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def log(x):print(time.strftime('%F %T'),x,flush=True)
COLS=['spot_uid','dataset_id','section_id','barcode','patient_uid','source_patient_identity_known','array_row','array_col','x','y','epithelial_fraction','ring1_valid_n','ring2_valid_n','complete_rings','reference','roi_id','copykat_prediction_original','cnv_score_original']
def inputs():
 m=pd.read_parquet(OLD/'03_features/common_anchors.parquet',columns=COLS);m['sample_uid']=m.dataset_id+'::'+m.section_id
 a=np.load(OLD/'03_features/features.npz');f={k:a[k] for k in a.files};return m,f,f['reference']
TYPES=json.loads((OLD/'03_features/feature_columns.json').read_text())['q1']
