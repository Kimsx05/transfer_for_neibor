from common import *
from PIL import Image
import re
figs=read(O/'figure_index.tsv'); sm=read(O/'section_metrics.tsv')
audit=read(O/'barcode_audit.tsv').set_index('section_id')
status=[]
for s,expected in S.items():
 d=read(O/f'spots/{s}.tsv.gz');assert set(d.section_id)=={s} and d.barcode.is_unique
 assert len(d)==audit.loc[s,'n_expression'] and d['query'].sum()==expected
 for branch,p in [('main',d),('noTop20',read(O/f'spots/{s}.noTop20.tsv.gz'))]:
  assert p.barcode.is_unique and set(p.barcode)==set(d.barcode)
  scores=p.filter(regex='^support_')
  assert scores.shape[1]==10 and np.allclose(scores.sum(axis=1),1)
  assert np.array_equal(p.H05_conditional_support.isna(),p.Epi_support==0)
  assert ((p.Epi_support<.1)==p.low_denominator).all()
  feat=read(O/f'features/{s}_{branch}.tsv')
  assert feat.in_anchors.all() and feat.in_PCA.all()
  assert set(feat.gene)&set(TOP)==(set(TOP) if branch=='main' else set())
 f=figs[figs.section_id==s];assert len(f)==11
 for row in f.itertuples():
  with Image.open(O/row.png) as im:
   im.verify()
  assert (O/row.pdf).read_bytes().startswith(b'%PDF')
  if row.plot not in ['tissue','scatter_H05_support','scatter_H05_conditional_support']:
   for suffix in ['_legend_only.pdf','_panel_only.pdf']:assert (O/(row.pdf[:-4]+suffix)).stat().st_size>1000
 assert len(sm[sm.section_id==s])==4
 status.append(dict(unit=s,status='complete',detail=f'{len(d)} mapped; {expected} frozen query; 2 branches; 11 standalone figures'))
pa=read(O/'pseudo_spot_predictions.tsv.gz')
assert len(pa)==480 and pa.pseudospot_id.is_unique
assert set(read(O/'features/pseudo_medium_main.tsv').gene)==set(read(O/'features/pseudo_low_main.tsv').gene)
status.append(dict(unit='pseudo_spot_validation',status='complete',detail='reused TEST, 240 mixtures x 2 matched depths; 3 specimens, 2 held-out source groups'))
for f in ['REPORT.md','REPORT.html','METHODS.md','README.md','input_manifest.tsv','reference_sampling_ids.tsv.gz','composition_depth_diagnostics.tsv','feature_ablation.tsv','pseudo_spot_validation.tsv']:
 assert (O/f).stat().st_size>100
resources=[]
for f in sorted((O/'logs').glob('*.resources')):
 t=f.read_text();m=re.search(r'Maximum resident set size \(kbytes\): (\d+)',t)
 if m:
  resources.append(dict(stage=f.stem,peak_RSS_KiB=int(m[1]),elapsed=re.search(r'Elapsed \(wall clock\) time .*?: (.+)',t)[1]))
save(resources,'resources.tsv');save(status,'status.tsv')
print(pd.DataFrame(status).to_string(index=False))

