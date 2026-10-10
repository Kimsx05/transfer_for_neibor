from pathlib import Path
import subprocess,time,csv,datetime
r=Path(__file__).resolve().parents[1]
p=r/'logs/GPU_resource_history.csv'
with p.open('a',buffering=1) as f:
 w=csv.writer(f)
 if p.stat().st_size==0:w.writerow(['timestamp','gpu_util_percent','used_MiB','total_MiB','GPU_processes'])
 while True:
  q=subprocess.run(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,memory.total','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True).stdout.strip().split(', ')
  procs=subprocess.run(['nvidia-smi','--query-compute-apps=pid,used_memory','--format=csv,noheader,nounits'],capture_output=True,text=True,check=True).stdout.strip().replace('\n','; ')
  w.writerow([datetime.datetime.now().isoformat(),*q,procs])
  if all((r/f'05_deconvolution/{b}_TEST/MAPPING_DONE.json').exists() for b in ['MAIN','NO_P07T']):break
  time.sleep(30)
