# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
from pathlib import Path
import requests,hashlib,json,time
root=Path(__file__).resolve().parents[1];dest=root/'.pip_cache/torch-2.8.0+cpu-cp312-cp312-manylinux_2_28_x86_64.whl'
url='https://download.pytorch.org/whl/cpu/torch-2.8.0%2Bcpu-cp312-cp312-manylinux_2_28_x86_64.whl'
expected='cb9a8ba8137ab24e36bf1742cb79a1294bd374db570f09fc15a5e1318160db4e'
s=requests.Session();s.trust_env=False;s.proxies={'https':'http://10.10.10.215:7897'}
r=s.get(url,headers={'Range':'bytes=0-0'},timeout=40);assert r.status_code==206,r.status_code;total=int(r.headers['Content-Range'].split('/')[-1]);print('TOTAL',total,flush=True)
part=dest.with_suffix('.part');offset=part.stat().st_size if part.exists() else 0
with part.open('ab') as f:
 while offset<total:
  end=min(total,offset+(4<<20))-1
  for attempt in range(3):
   try:
    r=s.get(url+f'?chunk={offset}',headers={'Range':f'bytes={offset}-{end}'},timeout=(15,50));assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {offset}-{end}/{total}';assert len(r.content)==end-offset+1;break
   except Exception:
    if attempt==2:raise
  f.write(r.content);f.flush();offset=end+1;print(offset,total,flush=True)
h=hashlib.sha256(part.read_bytes()).hexdigest();assert h==expected,(h,expected);part.replace(dest)
(root/'metadata/torch_wheel_verification.json').write_text(json.dumps(dict(url=url,bytes=total,expected_sha256=expected,actual_sha256=h,verified=True),indent=2))
