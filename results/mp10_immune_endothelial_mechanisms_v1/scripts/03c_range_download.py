"""Four asynchronous IO streams; no worker processes or compute thread pool."""
import asyncio,aiohttp,pathlib,os,json,hashlib,math
P=pathlib.Path('resources/mp10_regulatory_priors');name='hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather';part=P/(name+'.part');out=P/name;state=P/'scenic_range_download_state.json';total=1327758248;url='https://zenodo.org/records/7562143/files/'+name+'?download=1'
async def main():
 if out.exists():print('Already downloaded');return
 if state.exists():st=json.loads(state.read_text())
 else:
  have=part.stat().st_size if part.exists() else 0
  step=math.ceil((total-have)/4);st={'url':url,'total':total,'original_contiguous_bytes':have,'ranges':[{'next':a,'end':min(a+step,total)-1} for a in range(have,total,step)]};state.write_text(json.dumps(st,indent=2))
 fd=os.open(part,os.O_RDWR|os.O_CREAT,0o664)
 timeout=aiohttp.ClientTimeout(total=None,sock_connect=30,sock_read=90)
 async with aiohttp.ClientSession(timeout=timeout,trust_env=True) as s:
  async def run(i,r):
   for attempt in range(4):
    try:
     if r['next']>r['end']:return
     start=r['next']
     async with s.get(url,headers={'Range':f"bytes={start}-{r['end']}"}) as res:
      assert res.status==206,(res.status,dict(res.headers));assert res.headers.get('Content-Range','').startswith(f'bytes {start}-'),res.headers
      async for chunk in res.content.iter_chunked(1024*1024):
       os.pwrite(fd,chunk,r['next']);r['next']+=len(chunk);state.write_text(json.dumps(st,indent=2))
     assert r['next']==r['end']+1;print('Completed range',i,flush=True);return
    except Exception as e:
     print('Retry range',i,attempt,str(e),flush=True)
     if attempt==3:raise
     await asyncio.sleep(5)
  await asyncio.gather(*(run(i,r) for i,r in enumerate(st['ranges'])))
 os.close(fd)
 h=hashlib.md5();sh=hashlib.sha256()
 with part.open('rb') as f:
  for z in iter(lambda:f.read(1024*1024),b''):h.update(z);sh.update(z)
 assert h.hexdigest()=='081f2b37de0c4281679ce2aaeca6432f',h.hexdigest();part.rename(out)
 log=json.loads((P/'scenic_mirror_manifest.json').read_text());log=[x for x in log if x['file']!=name];log.append(dict(file=name,url=url,md5=h.hexdigest(),sha256=sh.hexdigest(),bytes=total,status='completed',source='Zenodo author archive of cisTarget hg38 refseq-r80 mc9nr'));(P/'scenic_mirror_manifest.json').write_text(json.dumps(log,indent=2));print('COMPLETE',name,flush=True)
asyncio.run(main())
