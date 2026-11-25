import sys,time,csv,psutil
p=psutil.Process(int(sys.argv[1]))
with open('logs/simultaneous_memory.tsv','w',buffering=1) as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['epoch','live_workers','total_RSS_bytes'])
 while p.is_running() and p.status()!=psutil.STATUS_ZOMBIE:
  try:
   workers=[x for x in p.children(recursive=True) if x.name()=='R']
   rss=p.memory_info().rss+sum(x.memory_info().rss for x in workers)
   w.writerow([time.time(),len(workers),rss])
  except psutil.Error:pass
  time.sleep(3)

