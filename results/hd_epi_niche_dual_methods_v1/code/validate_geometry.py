from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1];types=json.loads((R/'type_mapping.json').read_text());rows=[]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 d=pd.read_parquet(R/s/'nodes.parquet');a=pd.read_parquet(R/s/'r50.parquet').set_index('bin_id');b=pd.read_parquet(R/s/'r200.parquet').set_index('bin_id');assert a.index.equals(b.index)
 for ref,cats in types.items():
  for r,p in [(50,a),(200,b)]:
   c=p[[ref+'_C_'+t for t in cats]];known=p[ref+'_n_known'];assert np.all(c.sum(1)==known);v=p[[ref+'_P_'+t for t in cats]];assert np.allclose(v[known>0].sum(1),1);assert v[known==0].isna().all().all()
  assert (a[ref+'_n_known']<=b[ref+'_n_known']).all()
 # Include lowest support/edge and random centers. Explicit non-chunked Euclidean enumeration in local halo.
 ids=list(a.n_analysis_neighbors.nsmallest(4).index)+a.sample(8,random_state=123).index.tolist()
 for bid in ids:
  z=d.set_index('bin_id').loc[bid];delta=d[['x_um','y_um']].to_numpy()-z[['x_um','y_um']].to_numpy(float);ds=np.einsum('ij,ij->i',delta,delta)
  for r,p in [(50,a),(200,b)]:
   inside=(ds<=(r+1e-7)**2)&(d.bin_id!=bid)&(d.component==z.component)&d.analysis_qc
   for ref,cats in types.items():
    for typ in cats:assert int(((d['type_'+ref]==typ)&inside).sum())==int(p.loc[bid,ref+'_C_'+typ])
 rows.append(dict(section=s,n_anchors=len(a),direct_edge_and_random_checks=len(ids),counts_P_and_nested_radii=True,center_excluded=True,components_separate=True))
pd.DataFrame(rows).to_csv(R/'tables/geometry_validation.tsv',sep='\t',index=False)
(R/'tables/api_validation.json').write_text(json.dumps({'toy_radius_center_unknown_all_unknown':'passed','real_ROI_full_vs_query_row_composition':'passed','same_seed_native_ARI':1.0,'source':'logs/verify_api.log; code/verify_api.py'},indent=2))
print('GEOMETRY_VALIDATED')
