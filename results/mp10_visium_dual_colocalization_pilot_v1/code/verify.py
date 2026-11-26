from shared import *
from PIL import Image
import re
fig=read(O/'figure_index.tsv')
cross=read(T/'four_way.tsv')
mp=read(T/'functional_correlations.tsv')
q=read(O/'frozen_query.tsv.gz')
status=[];rates=[]
for s in S:
 d=read(O/f'spots/{s}.tsv.gz');up=read(P/f'spots/{s}.tsv.gz')
 z=join(d[KEY+['query',MP]],up[KEY+['query',MP]])
 assert (z.query_x==z.query_y).all() and np.allclose(z[MP+'_x'],z[MP+'_y'],atol=1e-12,equal_nan=True)
 e=d[d['query']];assert len(e)==S[s] and set(e.barcode)==set(q.loc[q.section_id==s,'barcode'])
 assert np.max(abs(d[H]+d[N]-1))<1e-7
 cells=list(read(O/'annotation_mapping.tsv').query('subtype!="Epithelial"').subtype)
 assert np.max(abs(d[[scol(t) for t in cells]].sum(axis=1)-1))<1e-7
 assert np.isclose(d[[ccol(t) for t in cells]].sum(axis=1),d.c2l_mean_TME29).all()
 for t in cells:assert np.allclose(d[lcol(t)],np.log1p(d[ccol(t)]))
 for model,n in [('A',31),('B',30)]:
  g=read(T/f'{model}_pairwise.tsv.gz').query('section_id==@s');assert len(g)==n*(n-1)//2 and g.n.eq(len(e)).all()
 g=cross[(cross.section_id==s)&(cross.level=='subtype')]
 assert len(g)==4*29 and g.groupby('cell_type').n.nunique().eq(1).all()
 for t in cells:
  for combo,x,y in [('A',H,scol(t)),('B',MP,lcol(t)),('MP10_transfer',MP,scol(t)),('H05_c2l',H,lcol(t))]:
   r=g[(g.cell_type==t)&(g.combination==combo)].iloc[0]
   assert np.isclose(r.pearson,corr(e[x],e[y]),atol=1e-8,equal_nan=True)
 for x in [H,MP]:
  for mod in ['Synthesis','Homeostasis','Fatty_acid','Efflux']:
   for v in ['full','minus20']:
    c=mod+'__'+v;b=e[c]
    assert b.between(0,1).all()
    r=mp[(mp.section_id==s)&(mp.exposure==x)&(mp.module_score==c)].iloc[0]
    assert r.n==len(e) and np.isclose(r.pearson,corr(e[x],b),equal_nan=True)
 nfig=0
 for r in fig[fig.section_id==s].itertuples():
  with Image.open(O/r.png) as img:
   assert min(img.size)>400
   assert abs(img.info.get('dpi',[0])[0]-300)<1
   img.verify()
  assert (O/r.pdf).read_bytes().startswith(b'%PDF')
  nfig+=1
 assert nfig>=55
 status.append(dict(section_id=s,full=len(d),frozen_query=len(e),assigned_niche=int(e.niche.notna().sum()),TME_subtypes=29,completed_models=2,completed_function_variants=8,standalone_figures=nfig,status='verified'))
overview=fig[fig.section_id=='four_sections']
assert len(overview)==4
for row in overview.itertuples():
 with Image.open(O/row.png) as img:
  assert abs(img.info.get('dpi',[0])[0]-300)<1
  img.verify()
 assert (O/row.pdf).read_bytes().startswith(b'%PDF')
assert not fig.duplicated(['section_id','figure']).any()
assert len(read(O/'reference_manifest.tsv.gz'))==21718
for f in ['README.md','METHODS.md','REPORT.md','REPORT.html','figures.html']:
 assert (O/f).stat().st_size>1000
for f in (O/'logs').glob('*.resources'):
 t=f.read_text(); rss=re.search(r'Maximum resident set size \(kbytes\): (\d+)',t);elapsed=re.search(r'Elapsed \(wall clock\) time .*?: (.+)',t)
 if rss:rates.append(dict(stage=f.stem,peak_RSS_KiB=int(rss[1]),elapsed=elapsed[1],exit_status=re.search(r'Exit status: (\d+)',t)[1] if 'Exit status:' in t else 'unknown'))
save(status,'status.tsv');save(rates,'resources.tsv')
print(pd.DataFrame(status).to_string(index=False))

