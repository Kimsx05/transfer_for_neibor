"""Read final review archives, align each ID, and keep rendering distinct from viewing."""
from base6 import *
import sqlite3,io
from PIL import Image
from hashlib import sha256
out=[]
for s in sections:
 o=R/s;q=pd.read_parquet(o/'review_index.parquet');q['archive_shard']=np.arange(len(q))//1024
 for shard,z in q.groupby('archive_shard',sort=False):
  f=o/'card_archive'/f'{shard:04d}.sqlite'
  with sqlite3.connect('file:'+str(f)+'?mode=ro',uri=True) as c:
   stored={u[0] for u in c.execute('select id from cards')};ids=set(z.review_id)
   assert stored==ids,(s,shard,len(stored),len(ids))
   k,b=c.execute('select id,jpeg from cards order by id limit 1').fetchone()
   with Image.open(io.BytesIO(b)) as img:img.load();w,h=img.size;fmt=img.format
   assert w>=700 and h>=470
   out.append((s,shard,len(ids),k,fmt,w,h))
 q['individual_jpeg']=q.rendered;q['rendered']=True;q.to_parquet(o/'review_index.parquet',index=False)
tab(pd.DataFrame(out,columns=['section','shard','unique_cards','sample_id','encoding','width','height']),'archive_readback.tsv')
