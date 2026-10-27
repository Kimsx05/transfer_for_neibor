import sys,json,gzip
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'mp10_cnv_regions_interface_v1';sys.path.insert(0,str(ROOT/'spatial_ecology_top20_v2/runtime'))
import pandas as pd,numpy as np
from shapely.geometry import shape,Point,LineString
from shapely.ops import unary_union
p=pd.read_csv(O/'regions/position_regions.tsv.gz',sep='\t',low_memory=False);d=pd.read_csv(O/'interfaces/position_interface_distances.tsv.gz',sep='\t',low_memory=False);b=pd.read_csv(O/'interfaces/boundary_summary.tsv.gz',sep='\t',low_memory=False);f=pd.read_csv(O/'input_review/position_features.tsv.gz',sep='\t',low_memory=False);src=pd.read_parquet(ROOT/'spatial_ecology_top20_v2/00_inputs/frozen_spot_inputs.parquet');keys=['section_id','position_id'];ss=src.rename(columns={'barcode':'position_id'}).set_index(keys);pp=p.set_index(keys).loc[ss.index];ff=f.set_index(keys).loc[ss.index]
assert np.allclose(pp[['x','y']],ss[['x','y']],rtol=0,atol=1e-10)
assert np.allclose(ff.EpiMP10_Top20_AUCell,ss.EpiMP10_Top20_AUCell,rtol=0,atol=1e-14,equal_nan=True)
assert np.array_equal(pp.final_mask,ss.epi_fraction_gate.ge(.25))
assert b.loc[b.eligible_interface,'boundary_type'].eq('Epi-rich_region_interface').all()
results=[]
for sec,gg in d[d.signed_distance_available].groupby('section_id'):
 with gzip.open(O/'regions/masks'/f'{sec}.geojson.gz','rt') as z:geo=json.load(z)
 polys={x['properties']['region_id']:shape(x['geometry']) for x in geo['features']};bd=b[(b.section_id==sec)&b.eligible_interface];maxerr=0;tested=0
 for _,r in gg.sample(min(12,len(gg)),random_state=2026).iterrows():
  point=Point(r.x,r.y);group=bd[bd.region_id==r.region_id];lines=unary_union([LineString([(x.x1,x.y1),(x.x2,x.y2)]) for x in group.itertuples()]);got=point.distance(lines);maxerr=max(maxerr,abs(got-r.distance_to_boundary));assert abs(got-r.distance_to_boundary)<1e-7
  assert polys[r.region_id].covers(point)==(r.signed_distance<0)
  tested+=1
 results.append(dict(section_id=sec,n_random_positions=tested,max_independent_distance_error_pixel=maxerr,geometry_sign_consistency=True))
pd.DataFrame(results).to_csv(O/'logs/independent_geometry_checks.tsv',sep='\t',index=False)
(O/'logs/source_parity_validation.json').write_text(json.dumps(dict(n_positions=len(p),source_xy_unchanged=True,source_raw_MP10_unchanged=True,frozen_epi_gate_unchanged=True,independent_shapely_distance_check_positions=sum(x['n_random_positions'] for x in results),max_distance_error=max(x['max_independent_distance_error_pixel'] for x in results)),indent=2))
print('PASS',len(results),'sections; source parity and independent point-to-line checks')
