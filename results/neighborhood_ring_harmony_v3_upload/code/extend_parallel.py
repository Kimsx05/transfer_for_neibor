from base import *
from concurrent.futures import ThreadPoolExecutor
import subprocess,time
p=OUT/'02_R12_pilot/selected_harmony_parameters.json'
while not p.exists():time.sleep(5)
s=json.loads(p.read_text());tasks=[(ring,mode,'NONE' if mode=='NONE' else s[mode]['strength']) for ring in ['R1','R2'] for mode in ['NONE','DATASET','SAMPLE']]
def run(t):
 ring,mode,strength=t
 with open(OUT/'logs'/f'extend_{mode}_{ring}.log','w') as f:subprocess.run(['/usr/bin/python3.12',str(OUT/'code/cluster.py'),ring,mode,strength],stdout=f,stderr=subprocess.STDOUT,check=True)
 return t
with ThreadPoolExecutor(max_workers=3) as pool:
 for t in pool.map(run,tasks):log('Extension completed '+str(t))
