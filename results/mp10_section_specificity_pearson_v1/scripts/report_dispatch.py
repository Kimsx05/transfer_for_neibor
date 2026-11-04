from report import *
import time
if __name__=='__main__':
 pending=set(M.section_id);active={}
 with ProcessPoolExecutor(max_workers=2) as ex:
  while pending or active:
   for sec in M.section_id:
    if sec not in pending:continue
    if (O/'reports'/f'{sec}_figure_index.tsv').exists():pending.remove(sec);continue
    if len(active)>=2:break
    if (O/'tables'/f'{sec}_context_done.tsv').exists() and (O/'tables'/f'{sec}_status.tsv').exists():active[ex.submit(report,sec)]=sec;pending.remove(sec)
   for f in list(active):
    if f.done():print(f.result(),flush=True);active.pop(f)
   if pending or active:time.sleep(5)
