import os,sys,json,hashlib,time,platform
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');OUT=Path(__file__).resolve().parents[1]
V2=ROOT/'neighborhood_ring_comparison_v2';V3=ROOT/'neighborhood_ring_harmony_v3';V4=ROOT/'neighborhood_mp10_source_v4';V5=ROOT/'neighborhood_mp10_calibration_v5';V6=ROOT/'neighborhood_niche_celltype_screen_v6'
sys.path=[p for p in sys.path if 'dist-packages' not in p];sys.path[:0]=[str(V2/'runtime'),str(ROOT/'spatial_ecology_top20_v2/runtime'),str(ROOT/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')];sys.path.append(str(ROOT/'neighborhood_v1/.venv/lib/python3.12/site-packages'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(OUT/'logs/mplconfig')
import numpy as np,pandas as pd,scipy
from scipy import sparse
from scipy.stats import rankdata,wilcoxon
TYPES=json.loads((V2/'03_features/feature_columns.json').read_text())['q1'];EPI=TYPES.index('Epithelial');NON=[j for j in range(30) if j!=EPI];RINGS=['R1','R2','R3'];LAYERS=['spot_MP10_model_universe','epi_MP10_model_universe'];SCOPES=['ALL_DATASETS']+sorted(pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t').dataset_id.unique());TIERS=['POSITIVE_CORE','POSITIVE_BROAD','EXPLORATORY','DEPLETION_CONTROL'];SETS=['MAIN','COMPLETE_RINGS','TME_MATCHED_RAW','TME_CONDITIONAL']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n');tmp.replace(p)
def tab(p,x):
 p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def log(s):print(time.strftime('%F %T'),s,flush=True)
def bh(p,n):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ii=np.flatnonzero(np.isfinite(p));order=ii[np.argsort(p[ii])];q[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,len(order)+1))[::-1])[::-1]);return q

def signed(x):
 x=np.asarray(x,float);x=x[np.isfinite(x)];nz=x[x!=0];n=len(x);d=dict(n_units=n,n_nonzero=len(nz),n_zero=int((x==0).sum()),n_positive=int((x>0).sum()),n_negative=int((x<0).sum()),mean_effect=x.mean() if n else np.nan,median_effect=np.median(x) if n else np.nan,IQR=np.ptp(np.quantile(x,[.25,.75])) if n else np.nan,p=np.nan,W=np.nan,actual_method='NOT_TESTED',status='INSUFFICIENT_UNITS')
 if n>=5:
  if len(nz)==0:d.update(p=1.,W=0.,status='ALL_ZERO_DIFFERENCES')
  elif len(nz)<5:d['status']='INSUFFICIENT_NONZERO_PAIRS'
  else:
   method='exact' if n<50 and n==len(nz) and len(np.unique(abs(x)))==n else 'approx';r=wilcoxon(x,alternative='two-sided',zero_method='wilcox',method=method,correction=True);d.update(p=r.pvalue,W=r.statistic,actual_method=method+'; two-sided signed-rank; Wilcox zeros; continuity=True',status='CALCULABLE')
 return d

def stars(q):return 'N/A' if not np.isfinite(q) else '****' if q<.0001 else '***' if q<.001 else '**' if q<.01 else '*' if q<.05 else 'ns'
def add_bh(d,keys,n,qname):
 d=d.copy();d[qname]=d.groupby(keys+['scope'],dropna=False).p.transform(lambda a:bh(a,n));d[qname+'_8scopes']=d.groupby(keys,dropna=False).p.transform(lambda a:bh(a,n*8)) if keys else bh(d.p,n*8)
 ix=d.scope.ne('ALL_DATASETS');d.loc[ix,qname+'_7datasets']=d[ix].groupby(keys,dropna=False).p.transform(lambda a:bh(a,n*7)) if keys else bh(d.loc[ix,'p'],n*7)
 return d

def align(uid,target,x):
 ix=pd.Index(uid);assert ix.is_unique and pd.Index(target).is_unique,'DUPLICATE_UID';take=ix.get_indexer(target);assert (take>=0).all(),'MISSING_UID';return x[take]
def input_record(p,purpose):
 p=Path(p);return dict(path=str(p.resolve()),size=p.stat().st_size,mtime=pd.Timestamp(p.stat().st_mtime,unit='s',tz='UTC').isoformat(),sha256=sha(p),purpose=purpose)
def checkpoint(name):js('logs/'+name+'.json',dict(status=name,input_hash=sha(OUT/'00_audit/input_manifest.tsv'),spec_hash=sha(OUT/'00_audit/spec_frozen.json'),time_UTC=pd.Timestamp.now(tz='UTC').isoformat()))
def done(name):
 p=OUT/'logs'/f'{name}.json'
 if not p.exists():return False
 d=json.loads(p.read_text());assert d['input_hash']==sha(OUT/'00_audit/input_manifest.tsv') and d['spec_hash']==sha(OUT/'00_audit/spec_frozen.json'),'CHECKPOINT_INPUT_CONFLICT';return True

def pairs():return pd.read_csv(OUT/'00_approval/approved_pairs_v7.tsv',sep='\t')
def load():
 m=pd.read_parquet(OUT/'01_ring123/anchor_metadata.parquet');a=np.load(OUT/'01_ring123/ring123_arrays.npz');assert np.array_equal(a['spot_uid'],m.spot_uid);return m,a['q'],a['center_mean']
def masks(m,q):
 base=m.main_eligible.to_numpy();tme=base&(np.sum(q[:,:,NON],axis=2)>=.05).all(1)
 return dict(MAIN=base,COMPLETE_RINGS=base&m.complete123.to_numpy(),TME_MATCHED_RAW=tme,TME_CONDITIONAL=tme)
def queries():
 p=pairs();return [dict(query_id=f'all_c{j:02d}',range='ALL_ELIGIBLE_EPI',accepted_label=-99,celltype=TYPES[j],cell_index=j) for j in NON]+[dict(query_id=r.pair_id,range='WITHIN_SELECTED_NICHE',accepted_label=r.accepted_label,celltype=r.celltype,cell_index=r.cell_index) for r in p.itertuples()]
