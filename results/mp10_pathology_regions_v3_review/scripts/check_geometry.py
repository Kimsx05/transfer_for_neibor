from run_v3 import *
from shapely.geometry import Point,shape
square=np.array([[[0,0],[1000,0]],[[1000,0],[1000,1000]],[[1000,1000],[0,1000]],[[0,1000],[0,0]]],float)
pts=np.array([[500,500],[1100,500],[0,500],[1000,1000]])
d,_,_=nearest(pts,square);inside=np.array([Polygon(square[:,0]).contains(Point(x)) for x in pts]);signed=np.where(inside,-d,d)
assert np.allclose(signed,[-500,100,0,0])
# Adjacent 100-um lattice centers must have a shared contour 50 um from each center.
hex0=np.c_[VERT[:,0]/6,VERT[:,1]/(2*np.sqrt(3))]*100
hex1=hex0+[100,0];common=np.array(sorted(set(map(tuple,hex0))&set(map(tuple,hex1))))
assert len(common)==2 and np.allclose(nearest(np.array([[0,0],[100,0]]),common[None])[0],50)
assert np.isclose(Polygon(hex0).area,10000*np.sqrt(3)/2)
checks=[]
for p in (O/'by_dataset').glob('*/*/complete.json'):
 z=pd.read_parquet(p.parent/'spot_regions_v3.parquet');v=z.signed_distance_um
 assert not z.barcode.duplicated().any() and not z.main_region.isin(['Mixed','Other']).any()
 assert (v[z.lesion_location=='inside'].dropna()<0).all() and (v[z.lesion_location=='outside'].dropna()>0).all()
 assert (v.notna()==z.distance_status.eq('estimable')).all()
 assert z.MP10_original.notna().all()
 for f in json.loads((p.parent/'lesion_rois_um.geojson').read_text())['features']:assert shape(f['geometry']).is_valid
 checks.append({'section_id':z.section_id.iloc[0],'n_spots':len(z),'geometry_and_labels':'pass'})
pd.DataFrame(checks).to_csv(O/'tables/geometry_checks.tsv',sep='\t',index=False)
print('PASS square sign/contour/um, 100um pitch -> 50um edge distance; cohort checks',len(checks))
