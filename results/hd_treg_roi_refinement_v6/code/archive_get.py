from pathlib import Path
import sqlite3,sys
import pandas as pd
section,key,out=sys.argv[1:];root=Path(__file__).resolve().parents[1]/'proseg_observed'/section
idx=pd.read_parquet(root/'review_queue.parquet',columns=['review_id']).set_index('review_id',verify_integrity=True).index.get_loc(key)
shard=root/'card_archive'/f'{idx//1024:04d}.sqlite'
with sqlite3.connect(f'file:{shard}?mode=ro',uri=True) as c:
 row=c.execute('select jpeg from cards where id=?',(key,)).fetchone()
 if row is None:raise ValueError('Card not committed in this shard')
 from PIL import Image
 import io
 with Image.open(io.BytesIO(row[0])) as image:image.save(out)
