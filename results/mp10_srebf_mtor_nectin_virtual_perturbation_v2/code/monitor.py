from common import *
import signal, datetime
f=Path('logs/cgroup_usage.tsv')
with f.open('a',buffering=1) as o:
 if o.tell()==0:o.write('epoch\tmemory_bytes\tpeak_bytes\tcpu_usec\tdisk_bytes_periodic\n')
 i=0
 while not Path(os.environ.get('MP10_MONITOR_STOP','tmp/monitor.done')).exists():
  v=int((cg/'memory.current').read_text());peak=int((cg/'memory.peak').read_text());cpu=int((cg/'cpu.stat').read_text().split()[1]);disk=int(subprocess.run(['du','-sb','.'],stdout=subprocess.PIPE,stderr=subprocess.DEVNULL).stdout.split()[0]) if i%15==0 else ''
  o.write(f'{time.time()}\t{v}\t{peak}\t{cpu}\t{disk}\n')
  if v>=C['soft_stop_bytes']:
   for sc in cg.glob('run-*.scope'):subprocess.run(['systemctl','--user','stop',sc.name],check=False)
   break
  i+=1;time.sleep(4)

