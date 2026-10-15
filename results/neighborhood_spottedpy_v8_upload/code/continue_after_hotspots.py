from common import *
import psutil,subprocess
pid=int(sys.argv[1]);log('Waiting for active hotspot execution '+str(pid))
while psutil.pid_exists(pid):
 try:
  if psutil.Process(pid).status()==psutil.STATUS_ZOMBIE:break
 except psutil.NoSuchProcess:break
 time.sleep(10)
assert done('M03_ALL_HOTSPOTS_COMPLETE'),'Hotspot process ended without complete checkpoint; inspect logs/hotspots.log'
log('All hotspots complete; continuing relations, figures, full validation and packages')
subprocess.run([str(OUT/'code/run_spottedpy_v8.sh'),'--resume'],check=True)
