import os, json, hashlib, sys, time, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse as sp
D=Path(__file__).resolve().parents[1];os.chdir(D)
V=D.parent/'mp10_srebf_nectin_virtual_perturbation_v1';P=D.parents[1]
prior=P/'resources/mp10_regulatory_priors';C=json.loads(Path('config/main.json').read_text())
K=C['source'];E=C['endpoints'];cg=Path('/sys/fs/cgroup/user.slice')/f'user-{os.getuid()}.slice'/f'user@{os.getuid()}.service'/'mp10vptwo.slice'
def rd(f,**kw):return pd.read_csv(f,sep='\t',**kw)
def save(x,f):x.to_csv(f,sep='\t',index=False,na_rep='')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(2**22),b''):h.update(block)
 return h.hexdigest()
def stamp(files):
 return hashlib.sha256(json.dumps([(str(p),sha(p)) for p in files]).encode()).hexdigest()
def guard():
 if (cg/'memory.current').exists():
  v=int((cg/'memory.current').read_text())
  if v>=C['soft_stop_bytes']:raise MemoryError('task reached 190 GB')
  while v>=C['soft_pause_bytes']:
   time.sleep(2);v=int((cg/'memory.current').read_text())
   if v>=C['soft_stop_bytes']:raise MemoryError('task reached 190 GB')

