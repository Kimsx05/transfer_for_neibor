from common import *
import subprocess, concurrent.futures
def launch(job):
 log=Path('logs')/('fit_'+job+'.log')
 with log.open('w') as f:subprocess.run(['taskset','-c','0-9','bash','code/fit.sh',job],stdout=f,stderr=subprocess.STDOUT,check=True)
if __name__=='__main__':
 import sys
 if '--pilot-done' not in sys.argv:launch('pilot')
 meta=rd('cache/metadata.tsv.gz');m=meta[meta.main]
 jobs=['main','alt']+['drop_'+x for x in sorted(m.dataset.unique())]+['state_'+x for x in sorted(m.state.unique())]+['sub%02d'%i for i in range(1,21)]
 usage=Path('logs/usage_pilot.txt').read_text()
 import re
 rss=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',usage)[1])
 stage=rd('logs/time_pilot.tsv');t=stage.wall_s.sum();ng=rd('qc/jobs/pilot_coverage.tsv.gz').in_network.sum()
 assert rss<32*1024**2 and ng>0
 fac=(2000/ng)**3
 save(pd.DataFrame([dict(pilot_cells=2500,pilot_genes=ng,pilot_stage_wall_s=t,pilot_RSS_KiB=rss,formal_jobs=len(jobs),planned_concurrent_jobs=3,estimated_wall_s=t*fac*len(jobs)/3,estimate_basis='pilot stages x cubic gene ratio, conservative for sparse states; excludes patient export',working_memory_cap_per_job_GB=38,maximum_parallel_address_space_GB=114,additional_disk_budget_GB=5)]),'tables/resource_plan.tsv')
 def safe(job):
  try:launch(job);print(job,'complete',flush=True)
  except Exception as e:print(job,repr(e),flush=True);raise
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
  futures=[pool.submit(safe,j) for j in jobs]
  for f in futures:f.result()

