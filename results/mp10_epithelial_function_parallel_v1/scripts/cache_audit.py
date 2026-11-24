from fetch import items,C,O,digest
import csv,hashlib,lzma,json,gzip,zipfile
from pathlib import Path
rows=[]
extra=[('GSE149933','GSE149933_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149933/soft/GSE149933_family.soft.gz'),('Lawson2020','mValidationAnalysis.Rmd','https://raw.githubusercontent.com/NMikolajewicz/Lawson2020/82662ff8183307ec09439dc001834537ec00bda3/Catherine/mValidationAnalysis.Rmd')]
for d,n,u in items+extra:
 p=C/d/n; packed=False
 if not p.exists() and (p.parent/(n+'.xz')).exists():p=p.parent/(n+'.xz');packed=True
 if not p.exists():
  assert d=='GSE149933' and 'series_matrix' in n
  rows.append([d,n,u,'由完整 family SOFT 替代','','','','']);continue
 ch=digest(p,'sha256');md=hashlib.md5()
 if d=='DepMap_24Q4_v1' and n!='manifest.json':
  m=next(i for i in json.loads((p.parent/'manifest.json').read_text())['files'] if i['name']==n)
  size=0
  with (lzma.open(p,'rb') if packed else p.open('rb')) as f:
   for block in iter(lambda:f.read(1<<22),b''):md.update(block);size+=len(block)
  assert m['computed_md5']==md.hexdigest() and m['size']==size
  msg='官方 MD5 和原始长度通过'
 else:
  if n.endswith(('.zip','.xlsx')):
   with zipfile.ZipFile(p) as z:assert z.testzip() is None
  elif n.endswith('.gz'):
   with gzip.open(p) as z:
    while z.read(1<<20):pass
  msg='已校验可读并记录 SHA256'
 rows.append([d,n,u,msg,str(p),p.stat().st_size,ch,md.hexdigest() if msg.startswith('官方') else ''])
with (O/'tables/resources_final.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['resource','requested_filename','source','final_status','cache_path','stored_bytes','stored_sha256','uncompressed_official_md5']);w.writerows(rows)

