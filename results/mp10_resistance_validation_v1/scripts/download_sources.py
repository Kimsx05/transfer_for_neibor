from pathlib import Path
import json,hashlib,urllib.request,zipfile
B=Path(__file__).resolve().parents[1]
for d in json.loads((B/'config/download_manifest.json').read_text()):
 p=Path(d['path']);expected=d['sha256']
 if p.exists():
  if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise RuntimeError('Existing source changed: '+str(p))
  continue
 data=urllib.request.urlopen(d['url'],timeout=120).read()
 if hashlib.sha256(data).hexdigest()!=expected:raise RuntimeError('Upstream source changed; explicit freeze amendment required: '+d['url'])
 p.write_bytes(data)
z=zipfile.ZipFile(B/'00_audit/sources/AMC_supp.zip')
for n in ['mmc9.xlsx','mmc4.xlsx']:
 p=B/'00_audit/sources'/n
 if not p.exists():p.write_bytes(z.read(n))
print('All frozen sources verified')
