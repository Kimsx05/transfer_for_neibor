from pathlib import Path
import os,json,hashlib,datetime,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['CUDA_VISIBLE_DEVICES']='0'
import numpy as np,pandas as pd
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');V=P/'mp10_k32_reference_identifiability_v5';S=V/'scope_frozen_resume_v1';R6=P/'mp10_k32_reference_diagnostic_v6';R=P/'mp10_k32_lower80_spatial_v7'
K='Epi_K32_anchor';L='Epi_MP10_lower80';H='withheld_MP10high_nonK32'
def say(*x):print(datetime.datetime.now(datetime.timezone.utc).isoformat(),*x,flush=True)
def js(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist() if isinstance(x,np.ndarray) else str(x),ensure_ascii=False)+'\n');t.replace(p)
def wt(x,p):x.to_csv(p,sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def key(x):return hashlib.sha256(('\n'.join(map(str,x))+'\n').encode()).hexdigest()
def seed(purpose):return int.from_bytes(hashlib.sha256(('123|'+purpose).encode()).digest()[:4],'little')%(2**31-1)
def config():return json.loads((R/'CONFIG.yaml').read_text())
def outmap(name):return R/'03_dev_mapping'/name if name.startswith('DEV') else R/'04_spatial_mapping'/name.split('_')[0]/name
