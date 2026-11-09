from pathlib import Path
import json,sys,hashlib
r=Path(__file__).resolve().parents[1]
for s in ['Cervilla','D1_Johnson']:
 o=r/s
 f=o/'cache_identity.json'
 if f.exists():
  z=json.loads(f.read_text());src=r/z['source_relative'];st=src.stat();assert st.st_size==z['source_size'] and st.st_mtime_ns==z['source_mtime_ns'],(s,'source RDS changed; do not combine stale cached branches with updated input')
  assert hashlib.sha256((o/'frozen_HVG.tsv').read_bytes()).hexdigest()==z['HVG_list_sha256'],(s,'frozen HVG list changed')
 for tag,lam in [('l00',0),('l02',.2),('l08',.8)]:
  f=o/f'{tag}_complete.json'
  if not f.exists():continue
  z=json.loads(f.read_text());expected=dict(section=s,lambda_=lam,npcs=30,HVG=3000,seed=123,resolution=.5,use_agf=False,split_scale=False,k_param=20)
  expected.pop('seed') # seed is in UMAP and global input record
  for k,v in expected.items():
   key='lambda' if k=='lambda_' else k
   assert z[key]==v,(s,tag,key,'cache parameters differ')
  assert z['UMAP']['seed']==123 and z['UMAP']['epochs']==200
  assert (o/f'{tag}_embeddings.rds').exists() and (o/f'{tag}_clusters.tsv.gz').exists()
 print(s,'completed branch cache parameters checked')
