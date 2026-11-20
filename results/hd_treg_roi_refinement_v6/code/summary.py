from base6 import *
from importlib.metadata import version
# Summaries derived only from stored data.
rows=[];funnel=[];geo=[];unrel=[];coverage=[]
for s in sections:
 o=R/s;d=pd.read_parquet(o/'cells_scored.parquet');k=rd(o/'priority_dominant_successor.tsv');src=pd.read_parquet(o/'sources.parquet');cv=pd.read_parquet(o/'candidate_grid_coverage.parquet');reviews=pd.read_parquet(o/'review_index.parquet')
 for v,col in [('source','sources'),('roi','roi')]:
  if v=='source':
   t=src[col].str.split('|',regex=False).explode().value_counts().rename_axis('source').reset_index(name='n');t['section']=s;rows.append(t)
 for z in ['in_D','geometry_assignable','UMI_ok','direct_ok','major_T_ok','A','B','main','strict','extended','control','epi','score_ok']:
  funnel.append((s,z,int((d[z]&d.in_D).sum())))
 for n,zz in d[d.in_D].groupby('mixed_v6'):geo.append((s,n,len(zz)))
 for z in reviews.outcome.dropna().unique():
  f=reviews[reviews.outcome.eq(z)];unrel.append((s,z,len(f),f.visual_note.notna().sum()))
 for r in [50,200]:
  a=cv[f'usable_{r}']/cv[f'theoretical_{r}'];coverage.append((s,r,len(a),int((a>=.5).sum()),np.median(a),a.min()))
tab(pd.concat(rows,ignore_index=True),'screening_sources.tsv');tab(pd.DataFrame(funnel,columns=['section','condition_in_D','n']),'identity_conditions.tsv');tab(pd.DataFrame(geo,columns=['section','mixed_v6','n']),'mixed_states.tsv');tab(pd.DataFrame(unrel,columns=['section','outcome_before_visual','n','opened']),'review_status.tsv');tab(pd.DataFrame(coverage,columns=['section','radius','sources','input_support_at_least_half','median','min']),'coverage.tsv')
