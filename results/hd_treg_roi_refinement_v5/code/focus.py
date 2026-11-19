from pictures import *
v=View('D1_Johnson')
d=pd.read_parquet(v.out/'cells_scored.parquet')
q=d[d.group.eq('RNA_candidate')].sort_values('id');w=285;h=265;cv=Image.new('RGB',(5*w,2*h),'white');dt=ImageDraw.Draw(cv)
for i,z in enumerate(q.itertuples()):
 ox=(i%5)*w;oy=(i//5)*h
 cv.paste(Image.fromarray(v.sample(z.x,z.y,55,220,'edge')),(ox,oy))
 label(dt,(ox,oy+221),f'{z.id}  {z.pred}',12)
 label(dt,(ox,oy+238),f'T {z.T_UMI}  FOXP3 {z.FOXP3}  core {z.core_width}  pTreg {z.pTreg:.3f}',11)
cv.save(R/'figures/D1_Johnson_RNA_candidates.jpg',quality=93)
