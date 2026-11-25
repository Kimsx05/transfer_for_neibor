"""Reuse checked-in task inputs; retrieve individual official processed exports if missing."""
from fetch import *
from geo import run
import shutil
from urllib.request import urlopen
urls=D/'config/resource_manifest.tsv'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def prepare():
 for rec in csv.DictReader(urls.open(),delimiter='\t'):
  p=Path(rec['path'])
  if rec['role']=='TCGA cache' and not p.is_file():p=I/p.name
  if not p.is_file() and rec['url'].startswith('https://'):
   if p.name=='GPL6947_extracted.soft':
    transfer=I/'GPL6947_extracted.soft.downloading'
    with urlopen(rec['url'],timeout=60) as src,gzip.GzipFile(fileobj=src) as g,transfer.open('wb') as out:
     for ln in g:
      out.write(ln)
      if ln.startswith(b'!platform_table_end'):break
    if sha(transfer)!=rec['sha256']:raise ValueError('GPL6947 annotation hash mismatch')
    transfer.replace(p)
   else:p=get(rec['url'],p.name)
  if not p.is_file() or sha(p)!=rec['sha256']:raise ValueError('Missing or changed input '+str(p))
 run()
if __name__=='__main__':prepare()

