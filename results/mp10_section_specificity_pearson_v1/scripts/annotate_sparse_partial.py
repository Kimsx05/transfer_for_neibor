from aggregate import *
def main():
 a=pd.read_parquet(O/'tables/pearson_all.parquet')
 for c in ['partial_epi_depth_status','partial_front_status']:a[c]=a[c].astype(str)
 a.loc[a.n_partial.between(5,9),'partial_epi_depth_status']='sparse_n_lt10_descriptive_only'
 a.loc[a.n_front_subset.between(5,9),'partial_front_status']='sparse_n_lt10_front_adjustment_not_stable'
 a.to_parquet(O/'tables/pearson_all.parquet',index=False)
 for sec,z in a.groupby('section_id',observed=True):z.to_parquet(O/'tables'/f'{sec}_pearson_q.parquet',index=False)
 save(a[a.program.isin(CORE)],O/'tables/pearson_core.tsv.gz')
 print('Sparse partial-correlation labels updated; numerical results unchanged',flush=True)
if __name__=='__main__':main()
