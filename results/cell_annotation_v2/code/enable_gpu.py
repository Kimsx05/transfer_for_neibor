"""Reuse v1 Python packages; add only same-version GPU wheels in a v2 overlay."""
from pathlib import Path
import requests,json,re,concurrent.futures,zipfile,hashlib,os,threading,time
root=Path(__file__).resolve().parents[1];dest=root/'.gpu_packages';dest.mkdir(exist_ok=True);temp=Path('/tmp/visiumhd_v2_gpu_wheels');temp.mkdir(exist_ok=True)
def session():
 s=requests.Session();s.trust_env=False;s.proxies={'https':'http://10.10.10.215:7897'};return s
def metadata(n,v):
 cache=temp/(n+'-'+v+'.metadata.json')
 if cache.exists():return json.loads(cache.read_text())
 for attempt in range(8):
  try:
   ss=session()
   if attempt%2:ss.proxies={}
   x=ss.get(f'https://pypi.org/pypi/{n}/{v}/json',timeout=(15,45));x.raise_for_status();j=x.json();cache.write_text(json.dumps(j));return j
  except Exception:
   if attempt==7:raise
   time.sleep(2)
meta={} 
for n,v in [('tensorflow','2.20.0'),('torch','2.8.0')]:
 meta[n]=metadata(n,v)
requirements={'tensorflow':'2.20.0','torch':'2.8.0','triton':'3.4.0','nvidia-cuda-nvcc-cu12':'12.8.93'}
for req in meta['torch']['info']['requires_dist']:
 m=re.match(r'(nvidia-[\w-]+)==([\d.]+)',req)
 if m:requirements[m[1]]=m[2]
def install(kv):
 n,v=kv;mark=dest/(n+'.installed.json')
 if mark.exists():return json.loads(mark.read_text())
 s=session();j=meta.get(n)
 if not j:
  j=metadata(n,v)
 files=[z for z in j['urls'] if z['filename'].endswith('.whl') and 'x86_64' in z['filename'] and ('cp312-cp312' in z['filename'] or 'py3-none' in z['filename'])];assert len(files)==1,(n,[z['filename'] for z in files]);w=files[0];f=temp/w['filename'];mark=dest/(n+'.installed.json')
 if mark.exists():return json.loads(mark.read_text())
 print('FETCH',n,w['size'],flush=True)
 ledger=f.with_suffix('.ranges.json');lock=threading.Lock()
 if ledger.exists():state=json.loads(ledger.read_text())
 else:
  offset=f.stat().st_size if f.exists() else 0;state={'prefix':offset,'done':[]};ledger.write_text(json.dumps(state))
  with f.open('ab') as out:out.truncate(w['size'])
 ranges=[(a,min(a+8*1024*1024,w['size'])) for a in range(state['prefix'],w['size'],8*1024*1024)];fd=os.open(f,os.O_RDWR)
 def part(ab):
  a,b=ab
  if a in state['done']:return
  ss=session()
  for attempt in range(5):
   try:
    ss.proxies={'https':'http://10.10.10.215:7897'} if attempt%2==0 else {}
    x=ss.get(w['url'],headers={'Range':f'bytes={a}-{b-1}'},timeout=(15,60));assert x.status_code==206 and x.headers.get('Content-Range')==f'bytes {a}-{b-1}/{w["size"]}';assert len(x.content)==b-a;os.pwrite(fd,x.content,a)
    with lock:state['done'].append(a);ledger.write_text(json.dumps(state))
    return
   except Exception:
    if attempt==4:raise
    time.sleep(1)
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:list(e.map(part,ranges))
 finally:os.close(fd)
 h=hashlib.sha256()
 with f.open('rb') as prior:
  for b in iter(lambda:prior.read(1024*1024),b''):h.update(b)
 assert f.stat().st_size==w['size'] and h.hexdigest()==w['digests']['sha256']
 ledger.unlink()
 with zipfile.ZipFile(f) as z:z.extractall(dest)
 info=dict(package=n,version=v,wheel=w['filename'],url=w['url'],bytes=w['size'],pypi_digest_verified=True);mark.write_text(json.dumps(info));f.unlink();print('READY',n,flush=True);return info
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:results=list(e.map(install,requirements.items()))
(root/'metadata/gpu_overlay.json').write_text(json.dumps(results,indent=2));print('GPU_OVERLAY_READY',flush=True)
