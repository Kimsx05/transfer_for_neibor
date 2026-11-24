from fetch import *
import argparse
args=sys.argv[1:]
for d,n,u in items:
 if not any(tok in n or tok==d for tok in args):continue
 dst=C/d/n
 if dst.exists():continue
 dst.parent.mkdir(parents=True,exist_ok=True)
 with tempfile.NamedTemporaryFile(delete=False,dir=dst.parent,prefix='.'+n+'.repair.') as h: p=Path(h.name)
 urls=[u]
 if d=='DepMap_24Q4_v1' and n!='manifest.json':
  meta=next(a for a in json.loads((dst.parent/'manifest.json').read_bytes())['files'] if a['name']==n);urls+=[meta['downloadUrl' if 'downloadUrl' in meta else 'download_url'],u.replace('ndownloader.figshare.com/files','figshare.com/ndownloader/files')]
 if d=='HMGCS1_2026': urls+=['https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13396504/supplementaryFiles']
 if u.startswith('https://ftp.ncbi') and '/suppl/' in u:urls+=['https://www.ncbi.nlm.nih.gov/geo/download/?acc='+d+'&format=file&file='+n]
 for j in range(12 if d=="DepMap_24Q4_v1" else 6):
  url=urls[j%len(urls)];cmd=['/usr/bin/curl','-fsSL','--connect-timeout','8','--max-time','600','--speed-limit','2000','--speed-time','30','-o',str(p),url]
  if j%2==1:cmd.insert(1,'-C');cmd.insert(2,'-')
  r=sp.run(cmd,capture_output=True)
  try:
   if r.returncode:raise RuntimeError(r.stderr.decode()[-180:])
   if url.endswith('supplementaryFiles'):
    with zipfile.ZipFile(p) as zz: binary=zz.read(next(x for x in zz.namelist() if x.endswith(n)))
    p.write_bytes(binary)
   validate(p,n,d);os.replace(p,dst);print('verified',n,flush=True);break
  except Exception as e:print(n,j,repr(e),flush=True)
 if p.exists():p.unlink()

