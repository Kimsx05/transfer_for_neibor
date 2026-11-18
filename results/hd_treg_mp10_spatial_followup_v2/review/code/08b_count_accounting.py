from fu_common import *
import anndata as ad
x,d,ge,sy,_=v1.load8('Cervilla');di=pd.Series(np.arange(len(d)),index=d.bin_id)
a=ad.read_h5ad(A/'segmentation_assignment_v2/Cervilla/cells.h5ad',backed='r');cm=pd.read_parquet(R/'tables/Cervilla_cell_region_source_cells.parquet');ids=cm.cell_id.unique();ii=a.obs_names.get_indexer(ids);oo=np.argsort(ii);cx=a.X[ii[oo]].tocsr();ci=pd.Series(np.arange(len(ids)),index=ids[oo]);a.file.close()
m=pd.read_parquet(R/'tables/common_support_reconstruction_members.parquet');units=pd.read_csv(R/'tables/common_support_reconstructed_units.tsv',sep='\t').set_index('unit_id');bp=pd.read_parquet(V/'Cervilla/region_pairs.parquet');bp=bp[bp.T_grid_um.eq(32)&bp.radius_um.eq(200)&~bp.shared_original_bins];cp=pd.read_parquet(R/'tables/Cervilla_cell_region_pairs.parquet');cp=cp[cp.radius_um.eq(200)&~cp.shared_original_bins];em=pd.read_parquet(V/'Cervilla/Epi64_members.parquet');out=[]
for tid in pd.read_parquet(R/'tables/Cervilla_common_centers_r200.parquet').T_region_id:
 mm=m[m.T_region_id.eq(tid)];eids=bp.loc[bp.T_region_id.eq(tid),'Epi_region_id'];cids=cp.loc[cp.T_region_id.eq(tid),'cell_region_id'];pools={'supplied_bin8':em.loc[em.region_id.isin(eids),'bin_id'].unique(),'whole_cell':cm.loc[cm.region_id.isin(cids),'cell_id'].unique()}
 for source,xx,lookup in [('supplied_bin8',x,di),('whole_cell',cx,ci)]:
  v=mm[mm.source.eq(source)];src=pools[source];used=v.member_id.unique();assert set(used)<=set(src);total=np.asarray(xx[lookup.loc[src]].sum(0)).ravel();included=np.zeros(18085,np.int64);good=np.zeros(18085,np.int64)
  for uid,q in v.groupby('unit_id'):
   y=np.asarray(xx[lookup.loc[q.member_id]].sum(0)).ravel().astype(np.int64);included+=y
   if units.loc[uid,'paired_eligible']:good+=y
  excluded=np.asarray(xx[lookup.loc[list(set(src)-set(used))]].sum(0)).ravel();expected=np.asarray(xx[lookup.loc[used]].sum(0)).ravel();assert np.array_equal(included,expected)and np.array_equal(total,included+excluded)
  out.append(dict(T_region_id=tid,source=source,old_source_raw_UMI=int(total.sum()),geometry_retained_UMI=int(included.sum()),geometry_excluded_UMI=int(excluded.sum()),paired_score_retained_UMI=int(good.sum()),paired_score_excluded_UMI=int(included.sum()-good.sum()),per_gene_conserved=True,scope='per T neighborhood; overlapping neighborhoods are not additive independent samples'))
put(pd.DataFrame(out),'common_support_count_accounting.tsv');print('PER_GENE_ACCOUNTING_OK',len(out),flush=True)
