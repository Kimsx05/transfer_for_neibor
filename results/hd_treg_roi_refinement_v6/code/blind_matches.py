from base6 import *
out=[]
for s in sys.argv[1:] or sections:
 o=R/s;d=pd.read_parquet(o/'identity_frozen.parquet');notes=R/'tables/visual_notes.tsv'
 if notes.exists():
  f=rd(notes);n=f.loc[f.section.eq(s)&f.action.isin(['unusable','assignment_uncertain']),'id'];d.loc[d.id.isin(n),['extended','control','epi']]=False
 native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['x_um','y_um','component','eligible_pool']);d['scorable_epi']=d.epi&d.nFeature.ge(50)
 centres=d[d.extended|d.control]
 for r in [50,200]:
  rows=[]
  for k,z in centres.groupby('component'):
   if k<0:continue
   ep=d[d.scorable_epi&d.component.eq(k)];pool=native[native.component.eq(k)&native.eligible_pool];xy=z[['x','y']].to_numpy()
   z=z.copy();z['op']=cKDTree(ep[['x','y']]).query_ball_point(xy,r,return_length=True);z['cv']=np.minimum(1,cKDTree(pool[['x_um','y_um']]).query_ball_point(xy,r,return_length=True)*64/(np.pi*r*r));rows.append(z[(z.op>0)&(z.cv>=.5)])
  if not rows:continue
  q=pd.concat(rows);f=np.array([np.log1p(q.op),np.log1p(q.nCount),q.cv]).T
  if q.empty:continue
  sd=f.std(axis=0);sd[sd==0]=1;q[['f1','f2','f3']]=f/sd
  for z in q[q.extended].itertuples():
   for cx in ['T','CD4']:
    cs=q[q.control&q.component.eq(z.component)&(q.nonTreg_type.str.startswith('CD4') if cx=='CD4' else True)].copy();cs['delta_um']=np.hypot(cs.x-z.x,cs.y-z.y);cs=cs[cs.delta_um<=1000];cs['diff']=np.linalg.norm(cs[['f1','f2','f3']]-np.array([z.f1,z.f2,z.f3]),axis=1)
    for i,ctrl in enumerate(cs.sort_values(['diff','delta_um','id']).head(3).itertuples()):out.append((s,r,z.id,ctrl.id,i+1,cx))
put(pd.DataFrame(out,columns=['section','radius','case','control','rank','pool']),R/'tables/blind_comparators.tsv')
