from run_v3_2 import *
cols=['marker_epithelial']+['gene_'+g for g in GENES]+['malignant_fraction','native_spot_MSQ'];vals={x:[] for x in cols};cov=[]
for sec in M.index:
 d=pd.read_parquet(OLD/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_1.parquet',columns=cols)
 for f in cols:vals[f].append(d[f].dropna().to_numpy())
 cov.append({'section_id':sec,'epi_available_genes':';'.join(g for g in GENES if d['gene_'+g].notna().any()),**{g:('available' if d['gene_'+g].notna().any() else 'unavailable') for g in GENES}})
sc={f:{'limits':[0,1 if f=='malignant_fraction' else float(np.quantile(np.concatenate(v),.99))],'upper_limit_basis':'cohort finite values 99th percentile, identical across all sections','missing':'gray/unavailable, never zero'} for f,v in vals.items()}
(O/'color_scales.json').write_text(json.dumps(sc,indent=2));pd.DataFrame(cov).to_csv(O/'tables/marker_availability.tsv',sep='\t',index=False)
# Known geometry sign / units test, independent of real labels.
seg=np.array([[[0.,-1000.],[0.,1000.]]]);dd,_,_=nearest(np.array([[-250.,0.],[50.,0.]]),seg);assert np.allclose(dd,[250,50]);assert np.array_equal(np.where(np.array([True,False]),-dd,dd),[-250,50]);assert dd[0]>200 and dd[1]<=100
(O/'tables/geometry_check.json').write_text(json.dumps({'straight_contour_test':'passed','inside_signed_um':-250,'outside_signed_um':50,'zero':'contour, not spot center','seed':20261030},indent=2))
