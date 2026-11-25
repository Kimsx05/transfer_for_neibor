from common import *
from concurrent.futures import ThreadPoolExecutor
def run(pair):
 j,s=pair;guard()
 directory=Path('cache/oracle')/j/s; marker=directory/'complete.json'
 if marker.is_file():
  previous=json.loads(marker.read_text())
  inputs=['code/oracle.py','cache/base_edges.tsv.gz','tables/control_matching.tsv','config/main.json','config/input_fingerprints.tsv','tables/panel_members.tsv']
  if previous.get('signature')==stamp([Path(f) for f in inputs]) and all((directory/f).is_file() and sha(directory/f)==h for f,h in previous['outputs'].items()):return
 log=Path('logs')/f'oracle_{j}_{s}.log'
 with log.open('w') as f:subprocess.run(['/usr/bin/time','-v','-o','logs/oracle_'+j+'_'+s+'.usage','cache/co310/bin/python','code/oracle.py',j,s],stdout=f,stderr=subprocess.STDOUT,check=True)
 print(j,s,'ok',flush=True)
if __name__=='__main__':
 m=rd('cache/metadata.tsv.gz');run(('pilot','C0'))
 js=['main','seed456']+['drop_'+s for s in sorted(m.dataset.unique())]
 with ThreadPoolExecutor(max_workers=min(8,C['threads_total']-1)) as pool:
  tasks=[pool.submit(run,(j,s)) for j in js for s in sorted(m.state.unique())]
  for task in tasks:task.result()

