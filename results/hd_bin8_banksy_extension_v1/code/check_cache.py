from pathlib import Path
import sys,json,hashlib
r=Path(__file__).resolve().parents[1];s=sys.argv[1];o=r/s
v=dict(input=json.loads((o/'input_summary.json').read_text())['source_signature'],parameters_sha256=hashlib.sha256((r/'parameters.json').read_bytes()).hexdigest(),HVG_sha256=hashlib.sha256((o/'frozen_HVG.tsv').read_bytes()).hexdigest())
p=o/'cache_identity.json'
if p.exists():assert json.loads(p.read_text())==v,'Input/parameters/HVG changed; refusing stale results'
else:p.write_text(json.dumps(v,indent=2))
print('CACHE_IDENTITY_OK',s)
