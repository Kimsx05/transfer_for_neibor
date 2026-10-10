from common import *
import subprocess,time,psutil
jobs=[('S31','FIXED_CAPTURE'),('S30','FULL_DEV'),('S30','TARGET_UMI'),('S30','FIXED_CAPTURE')];children={}
def active():
 out=[]
 for p in psutil.process_iter(['pid','cmdline']):
  try:
   cmd=' '.join(p.info['cmdline'] or [])
   if R.name+'/code/02_merged_reference.py' in cmd or R.name+'/code/04_mapping.py' in cmd:
    if p.info['cmdline'] and 'python3.12' in p.info['cmdline'][0]:out.append(p.pid)
  except (psutil.NoSuchProcess,psutil.AccessDenied):pass
 return out
while jobs or children:
 for k,p in list(children.items()):
  if p.poll() is not None:
   say('CHILD_EXIT',k,p.returncode);children.pop(k)
   if p.returncode:js(R/'logs'/('ERROR_'+k+'.json'),{'exit_code':p.returncode,'time':datetime.datetime.now().isoformat()})
 try:used=int(subprocess.check_output(['nvidia-smi','--query-gpu=memory.used','--format=csv,noheader,nounits'],text=True).strip().splitlines()[0])
 except Exception:used=23040
 for job in list(jobs):
  ls,col=job;tag=ls+'_'+col
  if (R/f'05_deconvolution/{tag}/MAPPING_DONE.json').exists():jobs.remove(job);continue
  running=False
  for pid in active():
   try:
    args=psutil.Process(pid).cmdline()
    if args[-2:]==[ls,col]:running=True
   except psutil.NoSuchProcess:pass
  if running:continue
  if ls=='S30' and not (R/'03_merged_reference/REFERENCE_DONE.json').exists():continue
  if len(active())>=3 or used>14000:break
  log=open(R/'logs'/('04_'+tag+'.log'),'a');p=subprocess.Popen(['bash',str(R/'code/python.sh'),str(R/'code/04_mapping.py'),ls,col],stdout=log,stderr=subprocess.STDOUT);log.close();children[tag]=p;jobs.remove(job);say('LAUNCHED',tag,p.pid,'GPU_used_before',used);time.sleep(15);break
 time.sleep(10)
say('SCHEDULE_FINISHED')
