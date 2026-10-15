import os,sys,json,hashlib,time,platform
from pathlib import Path
ROOT=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');OUT=Path(__file__).resolve().parents[1]
V2=ROOT/'neighborhood_ring_comparison_v2';V3=ROOT/'neighborhood_ring_harmony_v3';V4=ROOT/'neighborhood_mp10_source_v4';V5=ROOT/'neighborhood_mp10_calibration_v5';V6=ROOT/'neighborhood_niche_celltype_screen_v6';V7=ROOT/'neighborhood_ring123_broad_v7'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(OUT/'logs/mplconfig');os.environ['NUMBA_CACHE_DIR']=str(OUT/'logs/numba_cache')
import numpy as np,pandas as pd,scipy
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon,rankdata
TYPES=json.loads((V2/'03_features/feature_columns.json').read_text())['q1'];EPI=TYPES.index('Epithelial');NON=[j for j in range(30) if j!=EPI];LAYERS=['spot_MP10_model_universe','epi_MP10_model_universe'];TRACKS=['SOURCE_P05','LOCAL_BH05'];KS=[10,6,18];COV=pd.read_csv(V2/'00_preflight/coverage63.tsv',sep='\t');SCOPES=['ALL_DATASETS']+sorted(COV.dataset_id.unique());VARIABLES=[('MP10_SPOT',LAYERS[0],'SOURCE'),('MP10_EPI',LAYERS[1],'SOURCE')]+[(f'CELL_{j:02d}',TYPES[j],'CELL') for j in NON]
FAMILY_N={'NICHE_PRIMARY':36,'CELL_PRIMARY_DISTANCE':58,'NICHE_OUTSIDE_SECONDARY':18,'CELL_OVERLAP_SECONDARY':58,'CELL_OUTSIDE_SECONDARY':58}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):
 p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n');tmp.replace(p)
def tab(p,x):
 p=OUT/p;p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(x).to_csv(p,sep='\t',index=False)
def log(s):print(time.strftime('%F %T'),s,flush=True)
def bh(p,n=None):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ii=np.flatnonzero(np.isfinite(p));order=ii[np.argsort(p[ii])];n=len(p) if n is None else n;q[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,len(order)+1))[::-1])[::-1]);return q

def signed(x):
 x=np.asarray(x,float);x=x[np.isfinite(x)];near=(abs(x)<=1e-12)&(x!=0);x=x.copy();x[abs(x)<=1e-12]=0;n=len(x);nz=x[x!=0];d=dict(n_units=n,n_nonzero=len(nz),n_zero=int((x==0).sum()),n_rounded_nearzero=int(near.sum()),n_positive=int((x>0).sum()),n_negative=int((x<0).sum()),mean_effect=x.mean() if n else np.nan,median_effect=np.median(x) if n else np.nan,IQR=np.ptp(np.quantile(x,[.25,.75])) if n else np.nan,p=np.nan,W=np.nan,actual_method='NOT_TESTED',status='INSUFFICIENT_UNITS')
 if n>=5:
  if len(nz)==0:d.update(p=1.,W=0.,status='ALL_ZERO_DIFFERENCES')
  elif len(nz)<5:d['status']='INSUFFICIENT_NONZERO_PAIRS'
  else:
   method='exact' if n<50 and n==len(nz) and len(np.unique(abs(x)))==n else 'approx';r=wilcoxon(x,alternative='two-sided',zero_method='wilcox',method=method,correction=True);d.update(p=r.pvalue,W=r.statistic,actual_method=method+';two-sided;Wilcox zeros;continuity=True',status='CALCULABLE')
 return d

def stars(q):return 'N/A' if not np.isfinite(q) else '****' if q<.0001 else '***' if q<.001 else '**' if q<.01 else '*' if q<.05 else 'ns'
def input_record(p,why):
 p=Path(p);return dict(path=str(p.resolve()),size=p.stat().st_size,sha256=sha(p),purpose=why)
def fingerprint():return hashlib.sha256((''.join(sha(OUT/p) for p in ['00_audit/input_manifest.tsv','00_audit/spec_frozen.json','00_audit/source_versions.json'])).encode()).hexdigest()
def checkpoint(name):js('logs/'+name+'.json',dict(status=name,fingerprint=fingerprint(),time_UTC=pd.Timestamp.now(tz='UTC').isoformat()))
def done(name):
 p=OUT/'logs'/f'{name}.json'
 if not p.exists():return False
 assert json.loads(p.read_text())['fingerprint']==fingerprint(),'CHECKPOINT_CHANGED_INPUT';return True

def load_section(sec):
 p=OUT/'01_geometry'/sec;a=np.load(p/'values.npz');m=pd.read_parquet(p/'nodes.parquet');assert np.array_equal(a['spot_uid'],m.spot_uid);g=sparse.load_npz(p/'R1.npz');return m,a,g

def nearest(xy,components,source,target):
 d=np.full(len(source),np.nan)
 for comp in np.unique(components[source]):
  si=np.flatnonzero(components[source]==comp);ti=target[components[target]==comp]
  if len(ti):d[si]=cKDTree(xy[ti]).query(xy[source[si]],k=1,workers=1)[0]
 return d

def induced_components(g,ids,minsize=5):
 labels=np.full(g.shape[0],-1,int);formal=np.zeros(g.shape[0],bool);records=[]
 if len(ids):
  n,lab=connected_components(g[ids][:,ids],directed=False)
  for c in range(n):
   members=ids[lab==c];labels[members]=c;formal[members]=len(members)>=minsize;records.append(dict(component_id=c,size=len(members),formal=len(members)>=minsize))
 return labels,formal,records
