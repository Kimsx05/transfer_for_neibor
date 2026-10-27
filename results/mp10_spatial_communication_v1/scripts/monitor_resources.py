from pathlib import Path
import subprocess,time,csv,json,os
O=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_spatial_communication_v1');rootpid=1993775;peak=0;psspeak=0
with open(O/'logs/commot_Inhouse_process_tree.tsv','w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['unix_time','n_processes','sum_rss_GiB_conservative','sum_pss_GiB','pids'])
 while Path(f'/proc/{rootpid}').exists():
  lines=subprocess.check_output(['ps','-u','t070721','-o','pid=,ppid=,rss='],text=True).splitlines();rows=[tuple(map(int,l.split())) for l in lines];pids={rootpid}
  for _ in range(4):pids.update(pid for pid,ppid,rss in rows if ppid in pids)
  chosen=[(pid,rss) for pid,ppid,rss in rows if pid in pids];rss=sum(r for _,r in chosen)/1024**2;pss=0
  for pid,_ in chosen:
   try:
    for line in Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines():
     if line.startswith('Pss:'):pss+=int(line.split()[1])/1024**2
   except (OSError,PermissionError):pass
  peak=max(peak,rss);psspeak=max(psspeak,pss);w.writerow([time.time(),len(chosen),rss,pss,';'.join(str(p) for p,_ in chosen)]);f.flush();time.sleep(10)
(O/'commot/Inhouse_resource_tree_summary.json').write_text(json.dumps({'sampling_seconds':10,'peak_sampled_sum_rss_GiB_conservative':peak,'peak_sampled_sum_pss_GiB':psspeak,'worker_processes':4,'BLAS_threads_each':1},indent=2))
