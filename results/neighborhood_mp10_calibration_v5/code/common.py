import os,sys,json,hashlib,time
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism'); OLD=ROOT/'neighborhood_ring_comparison_v2'; V3=ROOT/'neighborhood_ring_harmony_v3'; V4=ROOT/'neighborhood_mp10_source_v4'; OUT=ROOT/'neighborhood_mp10_calibration_v5'
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

from scipy.stats import rankdata,spearmanr
TOP=pd.read_csv(OLD/'00_preflight/MP10_top20.tsv',sep='\t').gene.tolist()
PAIRS=[(1,2),(1,6),(7,6)]
def safe(out,p):
 p=Path(p).resolve();assert p.is_relative_to(Path(out).resolve());p.parent.mkdir(parents=True,exist_ok=True);return p
def write_tsv(out,p,x):pd.DataFrame(x).to_csv(safe(out,p),sep='\t',index=False)
def write_json(out,p,x):js(safe(out,p),x)
def auc_parts(x,genes):
 x=np.asarray(x,float);ix=pd.Index(genes).get_indexer(TOP);assert (ix>=0).all();k=int(np.ceil(.05*len(genes)));den=sum(k-i for i in range(1,21));r=rankdata(-x,method='average',axis=1)[:,ix];c=np.maximum(0,k-r)/den;s=c.sum(1);s[x.sum(1)==0]=np.nan;return s,c,r

def params(section):
 import anndata as ad,torch
 p=ROOT/'cell2location_full_run/v2_30type/02_sections'/section;a=ad.read_h5ad(p/'spatial_posterior.h5ad');v=torch.load(p/'model/model.pt',weights_only=False,map_location='cpu')['attr_dict'];C=v['cell_state_df_'];assert list(C.index)==list(a.var_names);assert list(C.columns)==list(a.uns['mod']['factor_names']);s=a.uns['mod']['post_sample_means'];return a,C,s

def allocate(Y,C,w,mg,ambient):
 mu=w@np.asarray(C).T*np.asarray(mg)+ambient;num=w[:,8,None]*np.asarray(C)[None,:,8]*np.asarray(mg);assert not ((mu<=0)&(Y>0)).any();E=np.divide(Y*num,mu,out=np.zeros_like(Y,dtype=float),where=mu>0);assert np.isfinite(E).all() and (E>=0).all() and (E<=Y+1e-6).all();return E

def epi_cache():
 p=OUT/'01_score_fidelity/counts_cache';g=p.joinpath('genes.txt').read_text().splitlines();cells=p.joinpath('cells.txt').read_text().splitlines();X=sparse.csr_matrix((np.fromfile(p/'data.bin',dtype='<f8'),np.fromfile(p/'indices.bin',dtype='<i4'),np.fromfile(p/'indptr.bin',dtype='<i4')),shape=(len(cells),len(g)));md=pd.read_csv(p/'metadata.tsv.gz',sep='\t',low_memory=False);assert md.cell_id.tolist()==cells;return X,g,md
