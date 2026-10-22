from pathlib import Path
import requests,json
from concurrent.futures import ThreadPoolExecutor
S=requests.Session();S.trust_env=False;S.proxies={'https':'http://127.0.0.1:17890','http':'http://127.0.0.1:17890'}
urls={
'mpbc_paper.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13232474/fullTextXML',
'cervilla_github_tree.json':'https://api.github.com/repos/scervilla/SpatialBenchmarking/git/trees/main?recursive=1',
'cervilla_dataverse.json':'https://dataverse.csuc.cat/api/datasets/:persistentId/?persistentId=doi:10.34810/DATA3063',
'geo_299573_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9040nnn/GSM9040999/suppl/',
'geo_318916_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE318nnn/GSE318916/suppl/',
'geo_300692_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9066nnn/GSM9066380/suppl/'}
def f(kv):
 k,u=kv
 try:
  r=S.get(u,timeout=45);r.raise_for_status();Path('metadata',k).write_bytes(r.content);print(k,r.status_code,len(r.content),flush=True)
 except Exception as e: print(k,str(e),flush=True)
with ThreadPoolExecutor(max_workers=3) as p:list(p.map(f,urls.items()))
