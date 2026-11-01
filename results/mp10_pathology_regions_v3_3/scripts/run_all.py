from concurrent.futures import ProcessPoolExecutor,as_completed
from run_v3_3 import *
def one(sec):
 try:process(sec);return sec,True
 except Exception:
  (O/'logs'/f'{sec}.error').write_text(traceback.format_exc());return sec,False
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=2) as ex:
  fs=[ex.submit(one,s) for s in M.index];failed=[]
  for f in as_completed(fs):
   sec,ok=f.result()
   if not ok:failed.append(sec);print('FAILED',sec,flush=True)
 print('Failures',failed,flush=True)
 if failed:raise SystemExit(1)
