#!/usr/bin/python3
"""Single-owner resumable all-section queue. No model choices changed by scheduling."""
import os,sys,csv,json,time,subprocess,signal,fcntl,traceback
from pathlib import Path
sys.path.insert(0,'/usr/local/lib/python3.12/site-packages')
try: import psutil
except ImportError:
 sys.path.insert(0,'/home/data/t070721/codex_workspace/Bladder Metabolism/cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages');import psutil
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');O=P/'mp10_spatial_communication_v2_all_sections';V=P/'mp10_spatial_communication_v1';L=O/'logs';S=O/'scripts';L.mkdir(exist_ok=True)
lock=open(L/'scheduler.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);lock.write(str(os.getpid()));lock.flush()
def now():return time.strftime('%Y-%m-%dT%H:%M:%S%z')
def atomic_json(path,value):
 t=Path(str(path)+'.tmp');t.write_text(json.dumps(value,ensure_ascii=False,indent=2));t.replace(path)
def output(method,sec):return O/('commot' if method=='COMMOT' else 'spatialcellchat')/sec/('primary2dNN' if method=='COMMOT' else 'main')
def marker(t):return Path(t['output_dir'])/('COMPLETED.json' if t['method']=='COMMOT' else 'COMPLETED.txt')
fields='task_id method section_id status attempts reason cache_dir output_dir pid started_at ended_at elapsed_seconds peak_tree_rss_GiB peak_tree_pss_GiB cpu_budget postprocess_status postprocess_pid'.split()
old={x['task_id']:x for x in json.loads((L/'task_status.json').read_text())} if (L/'task_status.json').exists() else {}
rows=list(csv.DictReader(open(V/'input_review/candidate_sections.tsv'),delimiter='\t')); tasks=[]
for row in rows:
 for method in ['COMMOT','SpatialCellChat']:
  sec=row['section_id'];tid=method+'::'+sec;out=output(method,sec)
  t=old.get(tid,dict.fromkeys(fields,''));t.update(task_id=tid,method=method,section_id=sec,output_dir=str(out));t.setdefault('n_spots',int(row['n_spots']));t['n_spots']=int(row['n_spots']);t['attempts']=int(t.get('attempts') or 0)
  if not t.get('status'):t['status']='queued'
  if '--retry-failed' in sys.argv and t['status'] in ['failed','blocked']:
   t.update(status='queued',attempts=0,reason='explicit retry retaining native checkpoints')
  if '--retry-failed' in sys.argv and t.get('postprocess_status')=='failed':t.update(postprocess_status='pending',postprocess_attempts=0)
  if sec in ['GSE171351_A1','Inhousedata_T4469479']:
   source=V/out.relative_to(O);expected=source/('COMPLETED.json' if method=='COMMOT' else 'COMPLETED.txt')
   if not expected.exists():raise RuntimeError('Missing v1 completion: '+str(expected))
   out.mkdir(parents=True,exist_ok=True);atomic_json(out/'CACHE_REFERENCE.json',dict(source_dir=str(source),read_only=True,reason='same frozen inputs, parameters and aggregation; no reinference',v1_snapshot='8813dc280bcba59867a7afd49bf412805350922b'));t.update(status='reused',cache_dir=str(source),reason='verified v1 primary result reuse')
  elif marker(t).exists():t.update(status='completed',cache_dir=str(out))
  elif t['status']=='running':
   pid=int(t.get('pid') or 0)
   if pid and psutil.pid_exists(pid):raise RuntimeError(f'Previous worker still exists {pid}; inspect before restart')
   t.update(status='queued',reason='resume after interrupted worker; native checkpoint retained')
  if (out/'POSTPROCESS_COMPLETED.json').exists():t['postprocess_status']='completed'
  elif t.get('postprocess_status')=='running':t['postprocess_status']='pending'
  tasks.append(t)
# The first two newly inferred jobs deliberately cover a small and a larger section.
priority={('COMMOT','CanDis_2026_Sample11'):0,('SpatialCellChat','GSE238145_B10'):1}
tasks.sort(key=lambda t:(priority.get((t['method'],t['section_id']),2),t['n_spots'],t['section_id'],t['method']))
def publish():
 atomic_json(L/'task_status.json',tasks)
 tmp=L/'task_status.tsv.tmp'
 with open(tmp,'w') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(tasks)
 tmp.replace(L/'task_status.tsv')
env=os.environ.copy();env.update({k:'1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS']});env.update(SCC_WORKERS='6',PYTHONHASHSEED='0',MPLBACKEND='Agg',MPLCONFIGDIR=str(L/'mplconfig'))
running={};post=None;launched=0;started=time.time();resource_file=open(L/'process_tree_resources.jsonl','a',buffering=1)
def launch(t):
 global launched
 out=Path(t['output_dir']);out.mkdir(parents=True,exist_ok=True);attempt=t['attempts']+1
 cmd=([sys.executable,str(S/'commot_run.py'),t['section_id']] if t['method']=='COMMOT' else ['/usr/local/bin/Rscript',str(S/'scc_run.R'),t['section_id']])
 budget=5 if t['method']=='COMMOT' else 7;te=env.copy()
 if attempt>1:
  if t['method']=='COMMOT':cmd+=['--sequential-components'];budget=1
  else:te['SCC_WORKERS']='2';budget=3
 log=open(L/(t['method']+'_'+t['section_id']+f'_attempt{attempt}.log'),'a');p=subprocess.Popen(cmd,cwd=P,env=te,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 t.update(status='running',attempts=attempt,pid=p.pid,started_at=now(),ended_at='',reason='',cpu_budget=budget,cache_dir=str(out),peak_tree_rss_GiB=0.,peak_tree_pss_GiB=0.)
 running[t['task_id']]=(p,log,time.time());launched+=1;print(now(),'START',t['task_id'],p.pid,flush=True);publish()
def measure(p):
 rss=pss=cpu=0.;n=0
 try:
  procs=[psutil.Process(p.pid)];procs+=procs[0].children(recursive=True)
  for q in procs:
   try:
    m=q.memory_full_info();rss+=m.rss;pss+=getattr(m,'pss',m.rss);ct=q.cpu_times();cpu+=ct.user+ct.system;n+=1
   except (psutil.NoSuchProcess,psutil.AccessDenied):pass
 except psutil.NoSuchProcess:pass
 return dict(rss_GiB=rss/2**30,pss_GiB=pss/2**30,cpu_seconds=cpu,n_processes=n)
def startpost(t):
 log=open(L/(t['method']+'_'+t['section_id']+'_postprocess.log'),'a');p=subprocess.Popen([sys.executable,str(S/'postprocess.py'),t['method'],t['section_id']],cwd=P,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);t.update(postprocess_status='running',postprocess_pid=p.pid,postprocess_attempts=t.get('postprocess_attempts',0)+1);return(t,p,log)
publish();print(now(),'SCHEDULER START 126 tasks, initial 2 jobs; cap 24 CPU /128 GiB',flush=True)
while True:
 used=0;psstotal=0
 for t in tasks:
  if t['task_id'] not in running:continue
  p,log,ts=running[t['task_id']];r=measure(p);used+=int(t['cpu_budget']);psstotal+=r['pss_GiB'];t['elapsed_seconds']=round(time.time()-ts,1)
  for key,rkey in [('peak_tree_rss_GiB','rss_GiB'),('peak_tree_pss_GiB','pss_GiB')]:t[key]=max(float(t.get(key) or 0),r[rkey])
  resource_file.write(json.dumps(dict(timestamp=now(),task_id=t['task_id'],pid=p.pid,**r))+'\n')
  rc=p.poll()
  if rc is not None:
   log.close();del running[t['task_id']];t['ended_at']=now();t['pid']=''
   if rc==0 and marker(t).exists():t.update(status='completed',reason='',postprocess_status='pending')
   elif t['attempts']<2:t.update(status='queued',reason=f'Attempt {t["attempts"]} exit {rc}; retry with fewer processes and existing checkpoint')
   else:t.update(status='failed',reason=f'Both attempts failed; last exit {rc}; see per-attempt logs; recover with scheduler --retry-failed')
   print(now(),'EXIT',t['task_id'],rc,t['status'],flush=True)
 # Corrective resource stop affects only this queue's process groups; checkpoints remain reusable.
 if psstotal>120:
  largest=max((t for t in tasks if t['task_id'] in running),key=lambda t:float(t['peak_tree_pss_GiB']))
  p=running[largest['task_id']][0];print(now(),'MEMORY GUARD terminate own worker',largest['task_id'],flush=True);os.killpg(p.pid,signal.SIGTERM)
 if post:
  t,p,log=post;rc=p.poll()
  if rc is not None:
   log.close();t['postprocess_status']='completed' if rc==0 and (Path(t['output_dir'])/'POSTPROCESS_COMPLETED.json').exists() else ('pending' if t['postprocess_attempts']<2 else 'failed');t['postprocess_pid']='';post=None
 if not post:
  for t in tasks:
   if t['status'] in ['reused','completed'] and t.get('postprocess_status') not in ['completed','failed'] and (O/'input_review/sections'/t['section_id']/'READY.json').exists():post=startpost(t);break
 # Initial two workers are measured for at least 10 minutes. Increase automatically only with headroom.
 concurrency=2
 if time.time()-started>=600 and psstotal<60 and psutil.virtual_memory().available>160*2**30:concurrency=3
 used=sum(int(t['cpu_budget']) for t in tasks if t['task_id'] in running)
 for t in tasks:
  if len(running)>=concurrency:break
  if t['status']!='queued':continue
  inp=O/'input_review/sections'/t['section_id']
  if not (inp/'READY.json').exists():continue
  need=5 if t['method']=='COMMOT' else 7
  if used+need+3>24 or psstotal>85 or psutil.virtual_memory().available<100*2**30:continue
  launch(t);used+=need
 # Input completion sentinel distinguishes still-preparing from truly blocked.
 if (O/'input_review/INPUT_PREPARATION_COMPLETED.json').exists():
  for t in tasks:
   if t['status']=='queued' and not (O/'input_review/sections'/t['section_id']/'READY.json').exists():t.update(status='blocked',reason='Input preparation finished without READY; see section FAILED.json',ended_at=now())
 publish()
 if all(t['status'] in ['completed','reused','failed','blocked'] for t in tasks) and not running and not post and all(t.get('postprocess_status') in ['completed','failed'] for t in tasks if t['status'] in ['completed','reused']):break
 time.sleep(10)
atomic_json(L/'INFERENCE_QUEUE_FINISHED.json',dict(time=now(),statuses={k:sum(t['status']==k for t in tasks) for k in ['completed','reused','failed','blocked']},postprocess_failed=sum(t.get('postprocess_status')=='failed' for t in tasks)))
for script in ['common_lr.py','summarize.py','report.py']:
 log=open(L/(script+'.log'),'a');rc=subprocess.call([sys.executable,str(S/script)],cwd=P,env=env,stdout=log,stderr=subprocess.STDOUT);log.close()
 if rc:atomic_json(L/'FINALIZATION_FAILED.json',dict(script=script,returncode=rc,time=now()));sys.exit(rc)
atomic_json(L/'ALL_DELIVERABLES_FINISHED.json',dict(time=now(),tasks=126));print(now(),'ALL DELIVERABLES FINISHED',flush=True)
