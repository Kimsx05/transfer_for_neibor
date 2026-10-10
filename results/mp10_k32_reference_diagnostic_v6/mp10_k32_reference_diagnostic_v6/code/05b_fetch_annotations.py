from common import *
import urllib.request,re,concurrent.futures
O=R/'01_feature_audit/public_annotations';O.mkdir(exist_ok=True);jobs=[]
for gsm,sample in zip(range(8707234,8707240),['T1','T2','T3','N1','N2','N3']):
 name=f'GSM{gsm}_{sample}_features.tsv.gz';jobs.append((f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8707nnn/GSM{gsm}/suppl/{name}',name))
# Small public annotation resources only; never download expression matrices.
def fetch(job):
 url,name=job;p=O/name
 try:
  if not p.exists():
   content=urllib.request.urlopen(url,timeout=35).read();assert len(content)<15000000;p.write_bytes(content)
  return dict(url=url,file=str(p),bytes=p.stat().st_size,sha256=sha(p),status='DOWNLOADED',date='2026-10-10')
 except Exception as e:return dict(url=url,file=str(p),status='UNAVAILABLE',reason=str(e))
rows=list(concurrent.futures.ThreadPoolExecutor(max_workers=6).map(fetch,jobs))
url='https://www.10xgenomics.com/support/software/space-ranger/downloads'
try:
 html=urllib.request.urlopen(url,timeout=35).read().decode();(O/'10x_download_page.html').write_text(html);links=sorted(set(re.findall(r'https[^\s<>"\\]+(?:Probe_Set|probe_set)[^\s<>"\\]+?\.csv',html)));js(O/'10x_probe_links.json',links)
 for link in links:
  if ('v1.0' in link or 'v2.0' in link) and ('Human' in link or 'human' in link):rows.append(fetch((link,link.split('/')[-1])))
except Exception as e:rows.append(dict(url=url,status='UNAVAILABLE',reason=str(e)))
wt(pd.DataFrame(rows),O/'download_manifest.tsv');say('FETCH_DONE',rows)
