import os,sys,subprocess,json,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[3];O=ROOT/'results/mp10_pathology_regions_v3'
sys.path.insert(0,str(ROOT/'spatial_ecology_top20_v2/runtime'))
import pandas as pd
m=pd.read_csv(ROOT/'mp10_pathology_regions_v2/tables/section_manifest.tsv',sep='\t');c=json.loads((O/'config.json').read_text());secs=c['pilot_sections']+[s for s in m.section_id if s not in c['pilot_sections']]
def run(sec):
 with open(O/'logs'/f'{sec}.run.log','w') as log:
  cmd=['/usr/bin/time','-v','/usr/bin/python3.12',str(O/'scripts/run_v3.py'),sec]
  if '--force' in sys.argv:cmd.append('--force')
  result=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
 print(time.strftime('%H:%M:%S'),sec,'exit',result.returncode,flush=True)
with ThreadPoolExecutor(max_workers=c['parallel_sections']) as pool:list(pool.map(run,secs))
