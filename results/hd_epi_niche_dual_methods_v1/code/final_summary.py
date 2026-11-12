from pathlib import Path
import pandas as pd,numpy as np,json
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
R=Path(__file__).resolve().parents[1];rows=[];cov=[];pair=[]
for ref in ['blca','hpca']:
 for radius in [50,200]:
  o=R/f'{ref}_r{radius}';j=json.loads((o/'complete.json').read_text());m=pd.read_csv(o/'candidate_metrics.tsv',sep='\t');d=pd.read_parquet(o/'labels.parquet')
  for method in ['A','B']:
   z=m[(m.method==method)&(m.parameter==(j['K'] if method=='A' else j['resolution']))].iloc[0];counts=d[method+'_assignment_source'].value_counts();rows.append(dict(reference=ref,radius=radius,method=method,selected_parameter=z.parameter,n_clusters=int(z.n_clusters),seed_ARI=z.seed_ARI,block_ARI=z.block_ARI,silhouette=z.silhouette,shared_anchors=len(d),evaluable=int(d.evaluable.sum()),native_train=int(counts.get('native_train',0)),external=int(counts.get('model_predict',0)+counts.get('external_projection',0)),unresolved=int(counts.get('unresolved',0)),complete=True))
   for s,g in d.groupby('section'):cov.append(dict(section=s,reference=ref,radius=radius,method=method,anchors=len(g),evaluable=int(g.evaluable.sum()),training=int((g[method+'_assignment_source']=='native_train').sum()),external=int(g[method+'_assignment_source'].isin(['model_predict','external_projection']).sum()),unresolved=int((g[method+'_assignment_source']=='unresolved').sum())))
for kind in ['scale','reference']:
 for val in (['blca','hpca'] if kind=='scale' else [50,200]):
  a=pd.read_parquet(R/(f'{val}_r50' if kind=='scale' else f'blca_r{val}')/'labels.parquet');b=pd.read_parquet(R/(f'{val}_r200' if kind=='scale' else f'hpca_r{val}')/'labels.parquet')
  for method in ['A','B']:
   for scope in ['common_evaluable_including_unresolved','common_assigned_both']:
    ok=a.evaluable&b.evaluable
    if scope=='common_assigned_both':ok&=(a[method+'_niche']!='unresolved')&(b[method+'_niche']!='unresolved')
    x=a.loc[ok,method+'_niche'];y=b.loc[ok,method+'_niche'];pair.append(dict(comparison=kind,reference_or_radius=val,method=method,scope=scope,n=int(ok.sum()),ARI=adjusted_rand_score(x,y),NMI=normalized_mutual_info_score(x,y)))
pd.DataFrame(rows).to_csv(R/'tables/main_combinations.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(R/'tables/coverage_wide.tsv',sep='\t',index=False);pd.DataFrame(pair).to_csv(R/'tables/paired_comparisons_explicit_denominators.tsv',sep='\t',index=False)
