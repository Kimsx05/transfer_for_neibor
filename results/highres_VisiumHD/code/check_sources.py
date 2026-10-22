"""Fetch small official metadata with explicit proxy selection; no shared env changes."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests,json,time
urls={'zenodo_record':'https://zenodo.org/api/records/17999961',**{a:'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+a+'&targ=self&form=text&view=full' for a in ['GSE300692','GSM9066380','GSE299573','GSM9040999','GSE318916','GSM9505917','GSM9505918']}}
def fetch(kv):
 k,u=kv
 for proxy in ['http://127.0.0.1:17890',None]:
  try:
   s=requests.Session();s.trust_env=False
   if proxy:s.proxies={'https':proxy,'http':proxy}
   r=s.get(u,timeout=(10,40));r.raise_for_status()
   Path('metadata/'+k+'.txt').write_text(r.text)
   print(k,r.status_code,len(r.content),proxy,flush=True)
   with Path('logs/network.jsonl').open('a') as f:f.write(json.dumps({'url':u,'status':r.status_code,'bytes':len(r.content),'proxy':proxy})+'\n')
   return
  except Exception as e:
   print(k,proxy,str(e),flush=True)
   with Path('logs/network.jsonl').open('a') as f:f.write(json.dumps({'url':u,'error':str(e),'proxy':proxy})+'\n')
with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(fetch,urls.items()))
