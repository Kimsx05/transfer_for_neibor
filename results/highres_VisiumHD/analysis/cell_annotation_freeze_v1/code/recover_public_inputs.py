# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import requests,json,datetime,hashlib
root=Path(__file__).resolve().parents[1]
urls={
'cervilla_zenodo':'https://zenodo.org/api/records/17999961',
'cervilla_github':'https://api.github.com/repos/scervilla/SpatialBenchmarking/git/trees/main?recursive=1',
'cervilla_loader':'https://raw.githubusercontent.com/scervilla/SpatialBenchmarking/main/preprocessing/load_visiumhd.R',
'cervilla_mirror':'https://dataverse.csuc.cat/api/datasets/:persistentId/?persistentId=doi:10.34810/DATA3063',
'cervilla_paper':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12888464/fullTextXML',
'uc_scc_sample':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9040999&targ=self&form=text&view=full',
'uc_scc_listing':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9040nnn/GSM9040999/suppl/',
'd1_sample':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9066380&targ=self&form=text&view=full',
'd1_listing':'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9066nnn/GSM9066380/suppl/',
'd1_paper':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12325965/fullTextXML',
'd1_author_tree':'https://api.github.com/repos/biopzhang/thor_scripts_ncomm/git/trees/main?recursive=1'}
def fetch(item):
 k,u=item;rpt={'key':k,'url':u,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 s=requests.Session();s.trust_env=False;s.proxies={'http':'http://127.0.0.1:17890','https':'http://127.0.0.1:17890'}
 try:
  r=s.get(u,timeout=(10,40));rpt.update(status=r.status_code,content_type=r.headers.get('Content-Type'),bytes=len(r.content),final_url=r.url)
  p=root/'metadata'/('public_'+k+'.txt');p.write_bytes(r.content);rpt['path']=str(p);rpt['sha256']=hashlib.sha256(r.content).hexdigest()
 except Exception as e:rpt['error']=str(e)
 print(json.dumps(rpt),flush=True);return rpt
if __name__=='__main__':
 with ThreadPoolExecutor(max_workers=4) as pool:res=list(pool.map(fetch,urls.items()))
 (root/'metadata/public_recovery_manifest.json').write_text(json.dumps(res,indent=2))
