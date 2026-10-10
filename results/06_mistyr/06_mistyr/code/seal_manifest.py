"""Seal output hashes after finalize.log is closed; do not hash a live log."""
import bootstrap
from bootstrap import ROOT
import hashlib,csv
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
rows=[]
for f in sorted(ROOT.rglob('*')):
 rel=f.relative_to(ROOT)
 if not f.is_file() or rel.parts[0] in ['runtime','sources'] or '__pycache__' in rel.parts or f.name=='output_manifest.tsv':continue
 role='model' if f.suffix=='.rds' else 'figure' if f.suffix in ['.png','.pdf'] else 'code' if rel.parts[0]=='code' else 'table' if f.suffix=='.tsv' else 'documentation_or_log'
 rows.append(dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f),role=role))
with open(ROOT/'output_manifest.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256','role'],delimiter='\t');w.writeheader();w.writerows(rows)
print('SEALED',len(rows),'outputs; runtime/source package files identified by source_registry.json')
