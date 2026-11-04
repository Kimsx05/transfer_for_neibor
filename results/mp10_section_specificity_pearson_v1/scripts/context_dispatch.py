from context import *
import time
if __name__=='__main__':
 pending=set(M.section_id);active={}
 with ProcessPoolExecutor(max_workers=3) as ex:
  while pending or active:
   for sec in M.section_id:
    if sec not in pending:continue
    if (O/'tables'/f'{sec}_context_done.tsv').exists():pending.remove(sec);continue
    if len(active)>=3:break
    if (O/'spots'/f'{sec}.parquet').exists():active[ex.submit(run_context,sec)]=sec;pending.remove(sec)
   for f in list(active):
    if f.done():print(f.result(),flush=True);active.pop(f)
   if pending or active:time.sleep(5)
 print('ALL CONTEXT COMPLETE',flush=True)
