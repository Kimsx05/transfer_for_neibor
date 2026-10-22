from pathlib import Path
import requests,json,xml.etree.ElementTree as E
from concurrent.futures import ThreadPoolExecutor
S=requests.Session();S.trust_env=False;S.proxies={'https':'http://127.0.0.1:17890','http':'http://127.0.0.1:17890'}
urls={
 'thor_paper.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12325965/fullTextXML',
 'thor_scripts_tree.json':'https://api.github.com/repos/biopzhang/thor_scripts_ncomm/git/trees/main?recursive=1',
 'cervilla_github_tree.json':'https://api.github.com/repos/scervilla/SpatialBenchmarking/git/trees/master?recursive=1',
 'thor_scripts_readme.txt':'https://raw.githubusercontent.com/biopzhang/thor_scripts_ncomm/main/README.md',
 'GSE299573_ena_runs.txt':'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRX29144058&result=read_run&fields=run_accession,study_accession,experiment_accession,fastq_ftp,fastq_bytes&format=tsv',
}
def f(kv):
 k,u=kv
 try:
  r=S.get(u,timeout=45);r.raise_for_status();Path('metadata',k).write_bytes(r.content);print(k,r.status_code,len(r.content),flush=True)
 except Exception as e:print(k,str(e),flush=True)
with ThreadPoolExecutor(max_workers=3) as p:list(p.map(f,urls.items()))
root=E.parse('metadata/mpbc_paper.xml').getroot()
text=' '.join(root.itertext())
for term in ['Data availability','Availability of data','GSE318916']:
 pos=text.find(term);print(term,pos,text[max(0,pos-250):pos+900])
