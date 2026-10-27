from runtime import *
import pandas as pd,numpy as np,json
from scipy import sparse
rows=[]
for out in (O/'commot').glob('*/*'):
 if not (out/'COMPLETED.json').exists():continue
 info=json.loads((out/'COMPLETED.json').read_text());sec=info['section'];xy=pd.read_csv(O/'input_review/sections'/sec/'coords.tsv',sep='\t')[['x_um','y_um']].to_numpy();idx=pd.read_csv(out/'position_lr_index.tsv',sep='\t');lr=pd.read_csv(out/'included_lr.tsv',sep='\t')
 for i,r in idx.iterrows():
  T=sparse.load_npz(out/r['file']).tocoo();pos=T.data>0;d=np.linalg.norm(xy[T.row[pos]]-xy[T.col[pos]],axis=1);lim=float(lr.iloc[i].distance_um);maxd=float(d.max()) if len(d) else 0.;assert maxd<=lim+1e-8
  rows.append(dict(section_id=sec,variant=info['variant'],lr_index=i,category=r.category,max_positive_link_um=maxd,declared_range_um=lim,n_positive_links=int(pos.sum()),status='PASS'))
pd.DataFrame(rows).to_csv(O/'commot/position_link_distance_checks.tsv',sep='\t',index=False);print('Distance checks PASS',len(rows))
