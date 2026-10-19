"""Bounded recovery for the already-authorized serial run, never bypass resource pauses."""
import os,time,json,subprocess,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
python=ROOT/'neighborhood_spottedpy_v8/env/bin/python';initial=int(sys.argv[1]) if len(sys.argv)>1 else 0
if initial:
 while Path(f'/proc/{initial}').exists():time.sleep(15)
for attempt in range(2):
 resume=json.loads((OUT/'RESUME.json').read_text()) if (OUT/'RESUME.json').exists() else {}
 if resume.get('status')=='EXECUTION_COMPLETE':break
 states=[]
 for p in (OUT/'combinations').glob('*/state.json'):
  try:states.append(json.loads(p.read_text()).get('stage'))
  except Exception:pass
 if 'RESOURCE_PAUSED' in states:
  print('RESOURCE_PAUSED: preserved checkpoint; no automatic retry or budget expansion',flush=True);break
 print('Bounded automatic recovery attempt',attempt+1,flush=True)
 with open(OUT/'logs/controller_recovery.log','ab') as log:
  proc=subprocess.run([str(python),str(OUT/'code/run.py'),'--resume'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
 print('recovery exit',proc.returncode,flush=True)
else:print('Finite retry allowance exhausted; inspect SEARCH_INCOMPLETE and combination error files',flush=True)
resume=json.loads((OUT/'RESUME.json').read_text()) if (OUT/'RESUME.json').exists() else {}
if resume.get('status')=='EXECUTION_COMPLETE':
 with open(OUT/'logs/final_validation.log','ab') as log:
  subprocess.run([str(python),str(OUT/'final_verify.py')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
 resume=json.loads((OUT/'RESUME.json').read_text())
print('Final saved status:',resume.get('status','NO_CHECKPOINT'),flush=True)
