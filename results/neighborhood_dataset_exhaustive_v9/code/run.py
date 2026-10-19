from core import *
import argparse,traceback,fcntl

def main():
 p=argparse.ArgumentParser();p.add_argument('--resume',action='store_true');p.add_argument('--combo');p.add_argument('--report-only',action='store_true');p.add_argument('--clean-temp',action='store_true');p.add_argument('--review-only',action='store_true');args=p.parse_args()
 lock=open(OUT/'logs/controller.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 from preflight import prepare
 prepare()
 from report import summarize
 if args.report_only:summarize();return
 if args.review_only:
  assert summarize();from review import run;run();summarize();return
 if args.clean_temp:clean();return
 from validate import validate
 validate()
 rows=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str}).sort_values('execution_order')
 if args.combo:rows=rows[rows.combo_id.eq(args.combo)];assert len(rows)==1
 codehash={str(p.relative_to(OUT)):sha(p) for p in (OUT/'code').glob('*') if p.is_file()}
 if (OUT/'audit/code_freeze.json').exists():assert json.loads((OUT/'audit/code_freeze.json').read_text())==codehash,'CODE_CHANGED_AFTER_LAUNCH: audit before resuming'
 else:js(OUT/'audit/code_freeze.json',codehash)
 from frozen import run as frozen
 from rings import run as ring
 from spotted import run as spotted
 for row in rows.itertuples():
  f=OUT/'combinations'/row.combo_id;f.mkdir(exist_ok=True)
  if valid(f,'COMPLETE'):continue
  start=time.time()
  try:
   guard(row.combo_id+':start',0);log('START '+row.combo_id+' '+row.datasets)
   if not valid(f,'NICHE_FROZEN'):
    with open(OUT/'logs'/f'{row.combo_id}_refit.log','a') as lf:proc=subprocess.run(['/usr/bin/python3',str(OUT/'code/refit.py'),row.combo_id],stdout=lf,stderr=subprocess.STDOUT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    if proc.returncode==75:raise ResourcePaused('representation guard requested pause')
    proc.check_returncode()
   ring(row.combo_id);spotted(row.combo_id);frozen(row.combo_id)
   for st in ['NICHE_FROZEN','RING_DONE','SPOTTED_MAIN_DONE','FROZEN_DONE']:assert valid(f,st)
   mark(f,'COMPLETE',[f/(st+'.json') for st in ['NICHE_FROZEN','RING_DONE','SPOTTED_MAIN_DONE','FROZEN_DONE']],dict(seconds_controller=time.time()-start));clean(row.combo_id);log('COMPLETE '+row.combo_id)
   if row.execution_order==0:
    used=bytes_used(f);js(OUT/'audit/pilot_resource_estimate.json',dict(actual_combo_bytes=used,upper_estimated_99_bytes=used*99,seconds=time.time()-start,rss_max_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,budget_bytes=json.loads((OUT/'config.json').read_text())['budget_bytes'],interpretation='pilot measured; smaller combinations have fewer sections; guards enforce actual cumulative usage'))
   # No provisional significance affects ordering or continuation.
  except ResourcePaused as e:
   js(f/'state.json',dict(stage='RESOURCE_PAUSED',error=str(e)));summarize();return
  except Exception as e:
   if isinstance(e,OSError) and e.errno in [28,122]:
    js(f/'state.json',dict(stage='RESOURCE_PAUSED',error=repr(e)));summarize();return
   error=traceback.format_exc();js(f/'state.json',dict(stage='TECHNICAL_FAILED',error=error));log('TECHNICAL_FAILED '+row.combo_id+' '+repr(e));(f/'error.txt').write_text(error)
  if row.execution_order%10==0:summarize()
 all_done=summarize()
 if all_done:
  from review import run
  run();summarize()
if __name__=='__main__':main()
