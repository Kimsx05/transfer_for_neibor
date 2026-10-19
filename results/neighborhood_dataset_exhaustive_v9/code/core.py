import os,sys,json,hashlib,time,ast,subprocess,shutil,resource
from pathlib import Path
OUT=Path(__file__).resolve().parents[1]; ROOT=OUT.parent
V2=ROOT/'neighborhood_ring_comparison_v2';V3=ROOT/'neighborhood_ring_harmony_v3';V7=ROOT/'neighborhood_ring123_broad_v7';V8=ROOT/'neighborhood_spottedpy_v8'
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.update(MPLBACKEND='Agg',MPLCONFIGDIR=str(OUT/'tmp/mplconfig'),NUMBA_CACHE_DIR=str(OUT/'cache/numba'),TMPDIR=str(OUT/'tmp'),PYTHONDONTWRITEBYTECODE='1')
import numpy as np,pandas as pd
DATASETS=['CanDis__UNVERIFIED','GSE171351','GSE238145','GSE285715','GSE319536','Inhousedata__UNVERIFIED','NPJ__UNVERIFIED'];INHOUSE=DATASETS[5]
TYPES=json.loads((V2/'03_features/feature_columns.json').read_text())['q1'];NON=[i for i,t in enumerate(TYPES) if t!='Epithelial'];EPI=TYPES.index('Epithelial');LAYERS=['spot_MP10_model_universe','epi_MP10_model_universe'];TRACKS=['SOURCE_P05','LOCAL_BH05']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,default=str).encode()).hexdigest()
def arrhash(x):return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
def js(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.pending');tmp.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n');tmp.replace(p)
def table(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);x=pd.DataFrame(x);tmp=p.with_name(p.name+'.pending')
 if p.suffix=='.parquet':
  x.to_parquet(tmp,index=False,compression='zstd',compression_level=5);import pyarrow.parquet as pq;assert pq.read_metadata(tmp).num_rows==len(x)
 else:x.to_csv(tmp,sep='\t',index=False)
 tmp.replace(p)
def npz(p,**kw):
 p=Path(p);tmp=p.with_name(p.name+'.pending')
 with open(tmp,'wb') as f:np.savez_compressed(f,**kw)
 with np.load(tmp) as a:assert set(a.files)==set(kw)
 tmp.replace(p)
def extract(path,names,ns):
 tree=ast.parse(Path(path).read_text());nodes=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in names];assert {n.name for n in nodes}==set(names)
 exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns);return ns
def log(s):print(time.strftime('%F %T'),s,flush=True)
def bytes_used(p=OUT):return sum(f.stat().st_size for f in p.rglob('*') if f.is_file() and not f.is_symlink())
def quota():
 try:
  raw=subprocess.check_output(['quota','-w'],text=True,stderr=subprocess.STDOUT,timeout=15)
  for line in raw.splitlines():
   s=line.split()
   if s and s[0]=='/dev/sda':return dict(status='KNOWN',remaining_bytes=max(0,(min(int(s[2]),int(s[3]))-int(s[1].rstrip('*')))*1024),raw=raw)
  return dict(status='UNKNOWN',raw=raw)
 except Exception as e:return dict(status='UNKNOWN',error=str(e))
class ResourcePaused(RuntimeError):pass
def guard(stage,peak=0):
 cfg=json.loads((OUT/'config.json').read_text());q=quota();free=shutil.disk_usage(OUT).free;available=min(free,q.get('remaining_bytes',free));used=bytes_used();cache=bytes_used(OUT/'cache');tmp=bytes_used(OUT/'tmp');figures=bytes_used(OUT/'review');r=dict(time=time.time(),stage=stage,used_bytes=used,available_bytes=available,quota_status=q['status'],predicted_peak_write=peak,rss_max_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,cache_bytes=cache,temp_bytes=tmp)
 p=OUT/'resource_usage.tsv';pd.DataFrame([r]).to_csv(p,sep='\t',index=False,mode='a',header=not p.exists())
 if used+peak>cfg['budget_bytes'] or available-peak<cfg['reserve_bytes'] or cache>cfg['cache_budget_bytes'] or tmp+peak>cfg['temp_budget_bytes'] or figures>cfg['figure_budget_bytes']:raise ResourcePaused(json.dumps(r))
 return r
def temp_path(combo,name):
 p=OUT/'tmp'/combo/name;p.parent.mkdir(parents=True,exist_ok=True)
 js(p.parent/'.active.json',dict(pid=os.getpid(),run_id='neighborhood_dataset_exhaustive_v9_20261018'))
 marker=json.loads((OUT/'tmp/.ownership.json').read_text());assert marker['run_id']==json.loads((OUT/'config.json').read_text())['run_id']
 reg=OUT/'tmp_registry.tsv';pd.DataFrame([dict(path=str(p.relative_to(OUT)),combo=combo,reason='regenerable representation exchange')]).to_csv(reg,sep='\t',index=False,mode='a',header=not reg.exists());return p
def clean(combo=None):
 reg=OUT/'tmp_registry.tsv'
 if not reg.exists():return
 marker=json.loads((OUT/'tmp/.ownership.json').read_text());assert marker['run_id']==json.loads((OUT/'config.json').read_text())['run_id']
 for r in pd.read_csv(reg,sep='\t').drop_duplicates('path').itertuples():
  if combo and r.combo!=combo:continue
  p=OUT/r.path
  if not p.exists():continue
  active=p.parent/'.active.json'
  if active.exists():
   pid=json.loads(active.read_text())['pid']
   if pid!=os.getpid():
    try:os.kill(pid,0)
    except ProcessLookupError:pass
    else:raise RuntimeError('TEMP_STILL_REFERENCED_BY_PID '+str(pid))
  assert not p.is_symlink() and p.resolve().is_relative_to((OUT/'tmp').resolve()) and p.is_file()
  z=dict(time=time.time(),path=r.path,bytes=p.stat().st_size,stage=r.combo,reason=r.reason)
  dest=OUT/'cleanup_manifest.tsv';pd.DataFrame([z]).to_csv(dest,sep='\t',index=False,mode='a',header=not dest.exists());p.unlink()
def config_hash():return sha(OUT/'config.json')
def mark(folder,stage,files,extra=None):
 d=dict(stage=stage,time=time.time(),config_hash=config_hash(),input_hash=sha(OUT/'audit/input_manifest.tsv'),files={str(Path(p).relative_to(OUT)):sha(p) for p in files},extra=extra or {})
 js(folder/(stage+'.json'),d);js(folder/'state.json',dict(stage=stage,time=d['time'],extra=d['extra']))
def valid(folder,stage):
 p=folder/(stage+'.json')
 if not p.exists():return False
 d=json.loads(p.read_text());assert d['config_hash']==config_hash() and d['input_hash']==sha(OUT/'audit/input_manifest.tsv'),'CHECKPOINT_CONFIG_OR_INPUT_CHANGED'
 assert all((OUT/f).exists() and sha(OUT/f)==h for f,h in d['files'].items()),'CHECKPOINT_ARTIFACT_MISSING_OR_CHANGED'
 if stage=='COMPLETE':
  for child in ['NICHE_FROZEN','RING_DONE','SPOTTED_MAIN_DONE','FROZEN_DONE']:assert valid(folder,child)
 return True
def correction(p,by=False):
 p=np.asarray(p,float);q=np.full(len(p),np.nan);ok=np.isfinite(p);a=np.where(ok,p,1.);order=np.argsort(a,kind='stable');n=len(a)
 if n:q[order]=np.minimum(1,np.minimum.accumulate((a[order]*n/np.arange(1,n+1)*(np.sum(1/np.arange(1,n+1)) if by else 1))[::-1])[::-1]);q[~ok]=np.nan
 return q
