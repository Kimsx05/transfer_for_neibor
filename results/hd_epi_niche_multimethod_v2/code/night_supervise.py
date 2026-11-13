from pathlib import Path
import os,subprocess,time,signal,json,psutil,threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from zoneinfo import ZoneInfo
R=Path(__file__).resolve().parents[1];py=str(R/'.venv/bin/python');active={};lock=threading.Lock()
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def allowed():return not 7<=datetime.now(ZoneInfo('Asia/Shanghai')).hour<23

def run(stage,offset=0,heavy=True):
 while heavy and not allowed():
  (R/'logs'/f'{stage}_status.json').write_text(json.dumps(dict(stage=stage,state='waiting_night_window',time=now())));time.sleep(30)
 with open(R/'logs'/f'{stage}.log','a') as log:
  p=subprocess.Popen([py,str(R/'code'/f'{stage}.py')],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env={**os.environ,'NICHE_CORE_OFFSET':str(offset)});paused=False
  with lock:active[stage]=p.pid
  while p.poll() is None:
   if heavy and not allowed() and not paused:os.killpg(p.pid,signal.SIGSTOP);paused=True
   if paused and allowed():os.killpg(p.pid,signal.SIGCONT);paused=False
   rss=0
   with lock:pids=list(active.values())
   for pid in pids:
    try:pp=psutil.Process(pid);rss+=sum(v.memory_info().rss for v in [pp]+pp.children(recursive=True) if v.is_running())
    except psutil.Error:pass
   state='paused_outside_window' if paused else 'running';(R/'logs'/f'{stage}_status.json').write_text(json.dumps(dict(stage=stage,state=state,pid=p.pid,combined_RSS_GB=rss/1e9,time=now())))
   if rss>195e9:
    for pid in pids:os.killpg(pid,signal.SIGTERM)
    raise RuntimeError('Combined task memory limit approached')
   time.sleep(15)
  with lock:active.pop(stage,None)
  (R/'logs'/f'{stage}_status.json').write_text(json.dumps(dict(stage=stage,state='complete' if p.returncode==0 else 'failed',returncode=p.returncode,time=now())))
  if p.returncode:raise RuntimeError(f'{stage} failed {p.returncode}; inspect log')
def composition():
 for stage in ['native_diagnostic','rings']:run(stage,0)
 with ThreadPoolExecutor(max_workers=4) as pool:
  ff=[pool.submit(run,f'fit_combo_{ref}_{rad}',offset) for (ref,rad,offset) in [('blca',50,0),('blca',200,4),('hpca',50,8),('hpca',200,12)]]
  for f in ff:f.result()
def molecular():
 for stage in ['expression','cellcharter_models']:run(stage,16)
if __name__=='__main__':
 errors=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(composition),pool.submit(molecular)]
  for f in futures:
   try:f.result()
   except Exception as e:errors.append(str(e))
 (R/'logs/runtime.json').write_text(json.dumps(dict(state='heavy_failed' if errors else 'heavy_complete',errors=errors,time=now())))
 for stage in ['mender_ring_sensitivity','summarize','figures','validate','report']:
  try:run(stage,0,stage in ['validate','mender_ring_sensitivity'])
  except Exception as e:errors.append(str(e))
 (R/'logs/runtime.json').write_text(json.dumps(dict(state='completed' if not errors else 'incomplete',errors=errors,time=now())))
 if errors:raise SystemExit(1)
