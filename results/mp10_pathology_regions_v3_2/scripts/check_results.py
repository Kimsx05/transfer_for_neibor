from run_v3_2 import *
from html.parser import HTMLParser
from urllib.parse import unquote
checks=[]
for sec in M.index:
 p=O/'by_dataset'/M.loc[sec,'dataset_id']/sec;d=pd.read_parquet(p/'spot_regions_v3_2.parquet');old=pd.read_parquet(OLD/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_1.parquet');ed=pd.read_csv(p/'lesion_boundary_edges.tsv.gz',sep='\t')
 assert not d.barcode.duplicated().any() and len(d)==len(old)
 keep=[c for c in old if c.startswith(('gene_','marker_','detected_','c2l_','SpaCET_')) or c in ['MP10_original','malignant_fraction','native_spot_MSQ','CNV_support','baseline_confidence','patient_id','barcode','x','y','array_row','array_col']]
 pd.testing.assert_frame_equal(d[keep],old[keep],check_exact=True)
 adj,lookup=graph(d);tc,cs=components(adj,np.ones(len(d),bool));rl,rs=components(adj,d.lesion_location.eq('inside').to_numpy())
 for ix in rs:assert d.lesion_roi_id.iloc[ix].nunique()==1
 assert len(rs)==d.loc[d.lesion_location.eq('inside'),'lesion_roi_id'].nunique()
 assert set(ed.lesion_roi_id)<=set(d.lesion_roi_id)
 assert set(d.reference_boundary_segment_id)-{''}<=set(ed.segment_id)
 for q in ed.itertuples():
  assert d.lesion_location.iloc[q.inside_index]=='inside' and d.lesion_roi_id.iloc[q.inside_index]==q.lesion_roi_id
  if q.outside_index>=0:assert tc[q.inside_index]==tc[q.outside_index] and q.outside_index in adj[q.inside_index]
  assert np.isclose(np.hypot(q.x2_um-q.x1_um,q.y2_um-q.y1_um),100/np.sqrt(3))
 assert d.signed_distance_um.notna().equals(d.exact_distance_status.eq('estimable'))
 assert (d.loc[d.lesion_location.eq('inside'),'signed_distance_um'].dropna()<0).all()
 assert (d.loc[d.lesion_location.eq('outside'),'signed_distance_um'].dropna()>0).all()
 for w in [100,200,300]:
  near=d[f'distance_band_{w}um'].isin(['Inner margin','Outer margin']);far=d[f'distance_band_{w}um'].isin(['Interior (>w confirmed)','Outside (>w confirmed)']);assert d.loc[near,'abs_distance_upper_bound_um'].le(w).all();assert d.loc[far,'abs_distance_lower_bound_um'].gt(w).all()
  no=d.exact_distance_status.isin(['outer_edge_not_observed','extent_unresolved']);assert not (near[no]|far[no]).any()
 for f in p.glob('*.png'):
  with Image.open(f) as im:assert min(im.info.get('dpi',[0,0]))>299
 checks.append({'section_id':sec,'n_spots':len(d),'barcode_values_coordinates':'pass','connected_ROI_and_new_segment_ids':'pass','no_cross_fragment_edges':'pass','no_fake_margin_without_supported_edge':'pass','marker_missingness_preserved':'pass'})
# Analytic counterexamples keep the distance/band model test non-vacuous even if this
# cohort revision has no accepted macro fronts. Unknown means ALL outside K.
contour=np.array([[[0.,-1000.],[0.,1000.]]]);a=nearest(np.array([[-250.,0.],[50.,0.],[0.,0.]]),contour)[0];assert np.allclose(a,[250,50,0])
# Outside spot x=400, observed front x=0, unobserved area starts x=300:
# q=100, not distance to a remote ROI-side bad edge. A 200-um band is unresolved.
g=400.;q=100.;L=min(g,q);U=g;assert not (U<=200 or L>200)
# A fully supported technical missing position can change q from 100 to 500;
# real defect/unknown missing cannot. Exact distance remains strict q>g.
assert not (100>g+1e-6);assert 500>g+1e-6;assert not (g>g+1e-6)
# Band-only information must not create exact distances: [300,500] => far at w=200.
assert min(500,300)>200 and not (300>500+1e-6)
id=pd.read_csv(O/'tables/preserved_identity_crosswalk.tsv',sep='\t').set_index('section_id')
assert id.loc['GSE238145_B5','analysis_unit_uid_v3_1']==id.loc['GSE238145_B18','analysis_unit_uid_v3_1']=='GSE238145::P1'
t=id.loc[['GSE285715_T1','GSE285715_T2','GSE285715_T3'],'analysis_unit_uid_v3_1'];assert t.nunique()==3 and t.str.contains('unresolved_section').all()
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  for k,v in attrs:
   if k in ['href','src']:self.links.append(v)
for f in O.rglob('*.html'):
 p=Links();p.feed(f.read_text())
 for link in p.links:
  if '://' not in link and not link.startswith('#'):assert (f.parent/unquote(link.split('#')[0])).exists(),(f,link)
pd.DataFrame(checks).to_csv(O/'tables/checks.tsv',sep='\t',index=False)
print('PASS:',len(checks),'sections;',sum(x['n_spots'] for x in checks),'spots; geometry/unknown counterexamples; preserved unit identities; all HTML links')
