from base_data import *
rows=[]
for s in HD:
 a=dataset(s);old=ad.read_h5ad(H/s/'cells.h5ad',backed='r')
 idx=pd.Index(a.obs.old_id[a.obs.source=='HE_primary']).get_indexer(old.obs_names)
 # Primary rows are located by the exact old ID, independent of index position.
 lut=pd.Series(np.arange(len(a)),index=a.obs_names);byold=pd.Series(np.arange(len(a))[a.obs.source.eq('HE_primary')],index=a.obs.loc[a.obs.source.eq('HE_primary'),'old_id']);ix=byold.loc[old.obs_names].values
 assert old.var_names.equals(pd.Index(a.var.feature_id))
 changes=0
 for i in range(0,len(ix),4096):
  delta=old.X[i:i+4096]-a.X[ix[i:i+4096]];delta.eliminate_zeros();changes+=delta.nnz
 dx=np.linalg.norm(old.obs[['x_um','y_um']].values-a.obs.iloc[ix][['x_um','y_um']].values,axis=1);assert changes==0
 rows.append((s,len(ix),changes,*np.quantile(dx,[.5,.99,1])))
write(pd.DataFrame(rows,columns=['sample','old_retained_matched','changed_count_entries','median_centroid_delta_um','p99_centroid_delta_um','max_centroid_delta_um']),'primary_exact_comparison.tsv')

