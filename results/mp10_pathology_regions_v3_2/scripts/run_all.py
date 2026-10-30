from concurrent.futures import ProcessPoolExecutor,as_completed
from run_v3_2 import *
def one(sec):
 try:process(sec);return sec,True
 except Exception:
  (O/'logs'/f'{sec}.error').write_text(traceback.format_exc());return sec,False
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as ex:
  futures=[ex.submit(one,s) for s in M.index]
  failures=[]
  for f in as_completed(futures):
   sec,ok=f.result()
   if not ok:failures.append(sec);print('FAILED',sec,flush=True)
 print('Failures:',failures,flush=True)
 if failures:raise SystemExit(1)
