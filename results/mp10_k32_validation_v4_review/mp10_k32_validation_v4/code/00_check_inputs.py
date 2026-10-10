from pathlib import Path
import csv,hashlib
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4')
f=R/'00_audit/INPUT_MANIFEST.tsv'
if f.exists():
 for r in csv.DictReader(f.open(),delimiter='\t'):
  p=Path(r['path']);s=p.stat();assert s.st_size==int(r['bytes']) and s.st_mtime_ns==int(r['mtime_ns']),str(p)+' changed; use a new run directory'
 print('Recorded input identity unchanged')
else:print('First run: identity will be frozen in the input audit')
