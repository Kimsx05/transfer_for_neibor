from pathlib import Path
import pandas as pd,json
r=Path(__file__).resolve().parents[1];hd=r.parents[1];rows=[]
for s in ['UC_SCC','MPBC1','MPBC2']:
 j=json.loads((r/s/'input_summary.json').read_text());rows.append(dict(section=s,role='original raw counts',path_relative_to_HD=j['source_signature']['path'],input_resolution_um=j['input_resolution_um'],tissue_mask_known=j['in_tissue_known'],software=j['software_version']))
 rows.append(dict(section=s,role='analysis counts loader source',path_relative_to_HD=str((r/s/j['counts_relative']).resolve().relative_to(hd)),input_resolution_um=8,tissue_mask_known=j['in_tissue_known'],software='native8um' if s=='UC_SCC' else 'fixed original grid raw sum'))
for s in ['Cervilla','D1_Johnson']:rows.append(dict(section=s,role='unchanged main bin object + BANKSY increments',path_relative_to_HD=f'analysis/hd_bin8_comparison_v1/{s}; analysis/hd_bin8_banksy_v1/{s}',input_resolution_um=8,tissue_mask_known=True,software='reuse; no model rerun'))
pd.DataFrame(rows).to_csv(r/'tables/input_manifest.tsv',sep='\t',index=False)
rows=[]
for s in ['UC_SCC','MPBC1','MPBC2']:
 for route in ['blca','hpca']:rows.append(dict(section=s,reference=route,**json.loads((r/s/(route+'_prediction_summary.json')).read_text())))
pd.DataFrame(rows).to_csv(r/'tables/new_section_reference_evaluability.tsv',sep='\t',index=False)
# Reuse two existing pipeline summaries; absent new-section metrics stay missing.
p=hd/'data/Cervilla/visiumhd_Bladder/web_summary.html';txt=p.read_text();d=json.JSONDecoder().raw_decode(txt.split('const data =',1)[1].lstrip())[0]
def walk(x):
 if isinstance(x,dict):
  for v in x.values():yield from walk(v)
 elif isinstance(x,list):
  if len(x)==2 and isinstance(x[0],str) and isinstance(x[1],str):yield tuple(x)
  for v in x:yield from walk(v)
c=dict(walk(d));p2=hd/'data/GSE300692/D1_Johnson/outs/metrics_summary.csv';dj=pd.read_csv(p2).iloc[0];rows=[]
for k in ['Number of Reads','Sequencing Saturation','Reads Mapped Confidently to Probe Set','Reads Mapped Confidently to the Filtered Probe Set']:
 rows.extend([dict(section='Cervilla',metric=k,value_as_reported=c[k],source=str(p.relative_to(hd)),scope='original pipeline library, not re-QC analysis bins'),dict(section='D1_Johnson',metric=k,value_as_reported=str(dj[k]),source=str(p2.relative_to(hd)),scope='original pipeline library; fractions decimal in source')])
for s in ['UC_SCC','MPBC1','MPBC2']:rows.append(dict(section=s,metric='pipeline sequencing/alignment metrics',value_as_reported='not supplied in located local inputs',source='targeted dataset directory search',scope='not inferred from expression counts'))
pd.DataFrame(rows).to_csv(r/'tables/available_pipeline_metrics.tsv',sep='\t',index=False)
