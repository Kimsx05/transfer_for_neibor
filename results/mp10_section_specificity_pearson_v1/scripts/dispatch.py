from core import *
import time,traceback
if __name__=='__main__':
 pending=set(M.section_id);active={}
 with ProcessPoolExecutor(max_workers=4) as ex:
  while pending or active:
   for sec in M.section_id:
    if sec not in pending:continue
    if (O/'tables'/f'{sec}_status.tsv').exists():pending.remove(sec);continue
    if len(active)>=4:break
    if (O/'scores'/f'{sec}_gene_coverage.tsv').exists():active[ex.submit(run,sec)]=sec;pending.remove(sec)
   for fut in list(active):
    if fut.done():
     sec=active.pop(fut)
     try:print(fut.result(),flush=True)
     except Exception:traceback.print_exc();raise
   if pending or active:time.sleep(5)
 print('ALL 63 COMPLETE',flush=True)
