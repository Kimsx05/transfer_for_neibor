from base6 import *
import hashlib
orig=ROOT/'proseg_pre_edge';n=rd(orig/'tables/visual_notes.tsv');p=R/'tables/visual_notes.tsv';cur=rd(p);rows=[];add=[]
for rec in n.to_dict('records'):
 s=rec['section'];name='cell_'+str(rec['id'])+'.jpg';a=orig/s/'cards'/name;b=R/s/'cards'/name
 if not (a.exists() and b.exists()):continue
 digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
 x=digest(a);y=digest(b);rows.append((s,rec['id'],str(a),str(b),x,y,x==y))
 if x==y and not ((cur.section==s)&(cur.id==rec['id'])).any():add.append(rec)
if add:put(pd.concat([cur,pd.DataFrame(add)],ignore_index=True),p)
tab(pd.DataFrame(rows,columns=['section','id','previous','final','sha_before','sha_final','identical']),'same_image_review.tsv')
