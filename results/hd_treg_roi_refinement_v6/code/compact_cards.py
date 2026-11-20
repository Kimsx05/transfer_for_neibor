"""Repack complete review shards using WebP without changing card dimensions or review records."""
import sys,sqlite3,io,os,multiprocessing
from pathlib import Path
from PIL import Image
import pandas as pd
root=Path(__file__).resolve().parents[1]/'proseg_observed'
def convert(args):
 f,keys=args
 with sqlite3.connect(f,timeout=40) as db:
  present={r[0] for r in db.execute('select id from cards')}
  if present!=set(keys):return
  if db.execute('select substr(jpeg,1,4) from cards limit 1').fetchone()[0]==b'RIFF':return
  tmp=f.with_suffix('.repack');tmp.unlink(missing_ok=True)
  with sqlite3.connect(tmp) as dest:
   dest.execute('create table cards (id text primary key,jpeg blob not null)')
   for key,pic in db.execute('select id,jpeg from cards'):
    with Image.open(io.BytesIO(pic)) as im:
     out=io.BytesIO();im.save(out,'webp',quality=75,method=0);dest.execute('insert into cards values (?,?)',(key,out.getvalue()))
  with sqlite3.connect(tmp) as dest:
   assert set(z[0] for z in dest.execute('select id from cards'))==present
  os.replace(tmp,f)
if __name__=='__main__':
 tasks=[]
 for s in ['Cervilla','D1_Johnson']:
  o=root/s;q=pd.read_parquet(o/'review_queue.parquet',columns=['review_id']).review_id.tolist()
  for f in sorted((o/'card_archive').glob('*.sqlite')):
   i=int(f.stem)
   if s=='D1_Johnson' and i>=125:continue
   tasks.append((f,q[i*1024:(i+1)*1024]))
 with multiprocessing.Pool(8) as pool:list(pool.imap_unordered(convert,tasks))
