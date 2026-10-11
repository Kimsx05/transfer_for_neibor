from common import *
import subprocess,sys,traceback
cfg=config();names=['GSE171351_B1','DEV','GSE238145_B17'];active={};errors={};runid='v7_'+str(int(cfg['t0_epoch']));js(R/'00_manifest/RUN_ID.json',{'run_id':runid,'launcher_pid':os.getpid(),'started':time.time()});log=R/'logs'
def start(name,mode='block'):
 f=open(log/(name+'.log'),'a');cmd=['bash',str(R/'code/python.sh'),str(R/'code/04_mapping.py'),name,mode];p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True);active[name]=(p,mode,time.time());say('LAUNCH',name,mode,p.pid);wt(pd.DataFrame([{'name':n,'pid':v[0].pid,'mode':v[1],'start_epoch':v[2]} for n,v in active.items()]),R/'00_manifest/active_pids.tsv')
def epoch(n):
 p=outmap(n)/'convergence.json';return json.loads(p.read_text()).get('epochs',0) if p.exists() else 0
def stable(n):
 p=outmap(n)/'MAPPING_DONE.json';return p.exists() and json.loads(p.read_text()).get('numerical_status')=='STABLE'
def gpu():
 try:
  s=subprocess.check_output(['nvidia-smi','--query-gpu=memory.used,memory.total,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip().split(',');return list(map(int,s))
 except:return [99999,0,100]
def refresh():
 if (R/'code/07_report.py').exists():
  with open(log/'report_updates.log','a') as f:subprocess.run(['bash',str(R/'code/python.sh'),str(R/'code/07_report.py'),'incremental'],stdout=f,stderr=subprocess.STDOUT)
while not (R/'01_reference/REFERENCE_DONE.json').exists():
 if time.time()>cfg['training_deadline']:raise RuntimeError('No reference at deadline')
 time.sleep(10)
js(R/'00_manifest/REFERENCE_FROZEN_FOR_MAPPING.json',{'signature_sha256':sha(R/'01_reference/signatures.tsv.gz'),'frozen_epoch':time.time(),'NB_status':json.loads((R/'01_reference/REFERENCE_DONE.json').read_text())})
while not (R/'02_pseudospots/SIMULATION_FROZEN.json').exists() or not (R/'00_manifest/SPATIAL_INPUTS_READY.json').exists():time.sleep(5)
while time.time()<cfg['training_deadline']:
 for n,(p,mode,st) in list(active.items()):
  if p.poll() is not None:
   if p.returncode:errors[n]=p.returncode;say('FAILED',n,p.returncode)
   else:say('FINISHED_BLOCK',n,epoch(n),time.time()-st)
   del active[n];refresh()
 pending=[n for n in names if n not in active and n not in errors and not stable(n) and epoch(n)<30000]
 if not pending and not active:break
 usage=gpu();maxjobs=2 if usage[1]-usage[0]>8500 and usage[2]<95 else 1
 if pending and len(active)<maxjobs:
  pending.sort(key=lambda n:(epoch(n),names.index(n)));start(pending[0])
 js(R/'STATUS.json',{'execution_status':'RUNNING','numerical_status':'MODEL_SPECIFIC','biological_evidence':'EVALUATION_PENDING','active':{n:{'pid':v[0].pid,'mode':v[1],'elapsed_s':time.time()-v[2]} for n,v in active.items()},'epochs':{n:epoch(n) for n in names},'errors':errors,'gpu':usage,'training_deadline':cfg['training_deadline'],'updated':time.time()});time.sleep(10)
# Each running callback exits at deadline, preserving current actual parameters.
while active:
 for n,(p,mode,st) in list(active.items()):
  if p.poll() is not None:
   if p.returncode:errors[n]=p.returncode
   del active[n];refresh()
 time.sleep(5)
for n in names:
 if n not in errors and (outmap(n)/'model/model.pt').exists() and not (outmap(n)/'MAPPING_DONE.json').exists():
  start(n,'final');p=active[n][0];p.wait();errors.update({n:p.returncode} if p.returncode else {});del active[n];refresh()
js(R/'00_manifest/DEFAULT_MAPPING_FINISHED.json',{'names':names,'errors':errors,'epochs':{n:epoch(n) for n in names},'finished':time.time(),'all_stable':all(stable(n) for n in names)});refresh();say('DEFAULT MAPPING FINISHED',errors)
