"""Render the entire frozen queue in resumable batches. Rendering never sets review state."""
from cards import *
import sqlite3,io,multiprocessing,functools
def part(b):
 keys=queue[b*1024:(b+1)*1024];db=folder/f'{b:04d}.sqlite';cn=sqlite3.connect(db);cn.execute('CREATE TABLE IF NOT EXISTS cards (id TEXT PRIMARY KEY, jpeg BLOB NOT NULL)');done=set(z[0] for z in cn.execute('SELECT id FROM cards'))
 for k in keys:
  if k in done:continue
  img=c.card(k);buf=io.BytesIO();img.save(buf,format='WEBP',quality=75,method=0);cn.execute('INSERT INTO cards VALUES (?,?)',(k,buf.getvalue()))
 cn.commit();assert cn.execute('SELECT count(*) FROM cards').fetchone()[0]==len(keys);cn.close();return b
if __name__=='__main__':
 import cards as module
 module.font=functools.lru_cache(maxsize=8)(module.font)
 s=sys.argv[1];folder=R/s/'card_archive';folder.mkdir(exist_ok=True);c=Cards(s);queue=list(c.q.index)
 with multiprocessing.get_context('fork').Pool(int(sys.argv[2])) as pool:
  for b in pool.imap_unordered(part,range(int(sys.argv[3]) if len(sys.argv)>3 else 0,(len(queue)+1023)//1024)):
   print('completed shard',b,flush=True)
