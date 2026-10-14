from common import *
from concurrent.futures import ThreadPoolExecutor
import subprocess

def run(scope):
 with (OUT/'logs'/('plots_'+scope+'.log')).open('a') as f:
  subprocess.run(['/usr/bin/python3.12',str(OUT/'code/plot_results.py'),'--scope',scope],stdout=f,stderr=subprocess.STDOUT,check=True)
 return scope
if done('M03_COMPOSITION_BOXES_DONE'):
 intact=True
 for scope in SCOPES:
  f=OUT/'06_figures'/scope/'display_name_map.tsv'
  if not f.exists():intact=False;continue
  d=pd.read_csv(f,sep='\t');intact=intact and len(d)==307 and all((OUT/x).exists() for x in d.png) and all((OUT/x).exists() for x in d.pdf) and (OUT/'06_figures'/scope/'overview.png').exists() and (OUT/'06_figures'/scope/'structural_Epi_TME.png').exists()
 if intact:log('All figure checkpoints intact; skipped');sys.exit(0)
with ThreadPoolExecutor(max_workers=4) as pool:
 for scope in pool.map(run,SCOPES):log('finished figures '+scope)
checkpoint('M03_COMPOSITION_BOXES_DONE');log('ALL SCOPE FIGURES COMPLETE')
