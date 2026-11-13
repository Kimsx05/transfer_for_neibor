import pathlib,requests,hashlib,json,time
p=pathlib.Path('resources/mp10_regulatory_priors');log=[]
files=[('motifs-v9-nr.hgnc-m0.001-o0.0.tbl','634bf3c4328d300549b58d391c24a6e7'),('hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather','081f2b37de0c4281679ce2aaeca6432f')]
for name,expected in files:
 url='https://zenodo.org/records/7562143/files/'+name+'?download=1';dest=p/name
 try:
  if not dest.exists():
   r=requests.get(url,timeout=(30,120),stream=True);r.raise_for_status()
   with open(str(dest)+'.part','wb') as f:
    for z in r.iter_content(1024*1024):f.write(z)
   pathlib.Path(str(dest)+'.part').rename(dest)
  h=hashlib.md5();sha=hashlib.sha256()
  with dest.open('rb') as f:
   for z in iter(lambda:f.read(1024*1024),b''):h.update(z);sha.update(z)
  assert h.hexdigest()==expected,(h.hexdigest(),expected)
  log.append(dict(file=name,url=url,md5=h.hexdigest(),sha256=sha.hexdigest(),bytes=dest.stat().st_size,status='completed',source='Zenodo author archive of cisTarget hg38 refseq-r80 mc9nr'))
 except Exception as e:log.append(dict(file=name,url=url,status='technically_blocked',reason=str(e)))
 (p/'scenic_mirror_manifest.json').write_text(json.dumps(log,indent=2));print(log[-1],flush=True)
