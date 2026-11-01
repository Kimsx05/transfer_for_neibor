from run_v3_3 import *
from html.parser import HTMLParser
from urllib.parse import unquote
rows=[]
for sec in M.index:
 p=O/'by_dataset'/M.loc[sec,'dataset_id']/sec;d=pd.read_parquet(p/'spot_regions_v3_3.parquet');old=pd.read_parquet(OLD/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_1.parquet');ed=pd.read_csv(p/'lesion_boundary_edges.tsv.gz',sep='\t');co=json.loads((p/'coordinate_scale.json').read_text())
 assert len(d)==len(old) and not d.barcode.duplicated().any()
 original=[c for c in old if c.startswith(('gene_','marker_','detected_','c2l_','SpaCET_')) or c in ['MP10_original','malignant_fraction','CNV_support','baseline_confidence','patient_id','barcode','x','y','array_row','array_col']]
 pd.testing.assert_frame_equal(d[original],old[original],check_exact=True)
 excluded=old.lesion_location.isin(['unlocalized','quality_excluded']);assert d.loc[excluded,'lesion_location'].equals(old.loc[excluded,'lesion_location'])
 adj,lookup=graph(d);tc,tcs=components(adj,np.ones(len(d),bool));rl,rs=components(adj,d.lesion_location.eq('inside').to_numpy())
 assert len(rs)==d.loc[d.lesion_location.eq('inside'),'lesion_roi_id'].nunique()
 for ix in rs:assert d.lesion_roi_id.iloc[ix].nunique()==1
 for e in ed.itertuples():
  assert d.lesion_location.iloc[e.inside_index]=='inside' and d.lesion_roi_id.iloc[e.inside_index]==e.lesion_roi_id
  if e.outside_index>=0:assert e.outside_index in adj[e.inside_index] and tc[e.inside_index]==tc[e.outside_index]
  assert np.isclose(np.hypot(e.x2_um-e.x1_um,e.y2_um-e.y1_um),100/np.sqrt(3))
 assert set(d.reference_boundary_segment_id)-{''}<=set(ed.segment_id)
 assert d.signed_distance_um.notna().equals(d.exact_distance_status.eq('estimable'))
 assert d.loc[d.lesion_location.eq('inside'),'signed_distance_um'].dropna().lt(0).all();assert d.loc[d.lesion_location.eq('outside'),'signed_distance_um'].dropna().gt(0).all()
 no=d.exact_distance_status.isin(['outer_edge_not_observed','unlocalized','quality_excluded'])
 for w in [100,200,300]:
  near=d[f'distance_band_{w}um'].isin(['Inner margin','Outer margin']);far=d[f'distance_band_{w}um'].isin(['Interior (>w confirmed)','Outside (>w confirmed)']);assert d.loc[near,'abs_distance_upper_bound_um'].le(w).all();assert d.loc[far,'abs_distance_lower_bound_um'].gt(w).all();assert not (near[no]|far[no]).any()
 # Every changed baseline position needs one explicit allowed local or author rule.
 changed=d.lesion_location.ne(d.v3_1_lesion_location);assert d.loc[changed,'change_class'].isin(['weak_borrowed_roi_to_surrounding_outside','author_coarse_region_corrected']).all()
 for i in np.flatnonzero(d.change_class.eq('author_coarse_region_corrected')):
  assert d.author_coarse_reference.iloc[i]==('author_coarse_inside' if d.lesion_location.iloc[i]=='inside' else 'author_coarse_outside')
 # Recompute selected nearest references from NEW segments, to catch stale v3.1/v3.2 distances.
 if sec in ['GSE238145_B16','GSE319536_BC13','GSE285715_T3']:
  ii=np.flatnonzero(d.distance_to_observed_outer_candidate_um.notna())[::37]
  for i in ii:
   ee=ed[ed.boundary_type.eq('candidate_observed_outer')]
   ee=ee[ee.lesion_roi_id.eq(d.lesion_roi_id.iloc[i])] if d.lesion_location.iloc[i]=='inside' else ee[ee.tissue_component.eq(tc[i])]
   g=nearest(d[['x_um','y_um']].iloc[[i]].to_numpy(),segments(ee))[0][0]
   assert np.isclose(g,d.distance_to_observed_outer_candidate_um.iloc[i])
 # Geometric support includes candidate outside even when v3.2 called it unresolved.
 good=ed[ed.boundary_type.eq('candidate_observed_outer')]
 if len(good):assert d.distance_to_observed_outer_candidate_um.notna().any()
 for f in p.glob('*.png'):
  with Image.open(f) as im:assert min(im.info.get('dpi',[0,0]))>299
 rows.append({'section_id':sec,'n_spots':len(d),'candidate_edges':len(good),'identities_coordinates_original_values':'pass','new_geometry_same_tissue_component':'pass','unlocalized_QC_preserved':'pass','position_change_basis':'pass','distance_band_consistency':'pass'})
assert sum(r['candidate_edges'] for r in rows)>0,'Whole cohort candidate edges vanished: position/confidence/geometry must be re-examined.'
# Focused independent sign/scale and unknown-region examples.
g=np.array([[[0.,-1000.],[0.,1000.]]]);dist=nearest(np.array([[-250.,0.],[50.,0.],[0.,0.]]),g)[0];assert np.allclose(dist,[250,50,0]);assert min(500,300)>200 and not (300>500+1e-6);assert not (400>400+1e-6)
sec='GSE238145_B16';d=pd.read_parquet(O/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_3.parquet');z=d[d.barcode.eq('AGATAACACGGTATAA-1_4')].iloc[0];assert z.lesion_location==z.v3_1_lesion_location=='outside' and z.v3_2_lesion_location=='inside';assert z.lesion_roi_id=='';assert abs(z.domain_local_evidence_score-.463318)<1e-5
identity=pd.read_csv(O/'tables/preserved_identity_crosswalk.tsv',sep='\t').set_index('section_id');assert identity.loc['GSE238145_B5','analysis_unit_uid_v3_1']==identity.loc['GSE238145_B18','analysis_unit_uid_v3_1']=='GSE238145::P1';assert identity.loc[['GSE285715_T1','GSE285715_T2','GSE285715_T3'],'analysis_unit_uid_v3_1'].nunique()==3
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  for k,v in attrs:
   if k in ['href','src']:self.links.append(v)
for p in O.rglob('*.html'):
 parser=Links();parser.feed(p.read_text())
 for link in parser.links:
  if '://' not in link and not link.startswith('#'):assert (p.parent/unquote(link.split('#')[0])).exists(),(p,link)
pd.DataFrame(rows).to_csv(O/'tables/checks.tsv',sep='\t',index=False)
print('PASS',len(rows),'sections',sum(r['n_spots'] for r in rows),'spots; B16 correction, original identities/values, geometry and all HTML links')
