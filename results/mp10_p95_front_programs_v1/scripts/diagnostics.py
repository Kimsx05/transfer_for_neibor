from analyze import *
summary=[];balance=[];seg=[];maprows=[];ringcov=[]
for sec in M.section_id:
 d=pd.read_parquet(O/'spots'/f'{sec}.parquet');r=M.set_index('section_id').loc[sec]
 dt=pd.read_csv(O/'scores'/f'{sec}_signature_detection.tsv.gz',sep='\t').merge(d[['barcode','query','tumor_epi','MP10_P95','front_band']],on='barcode',validate='m:1')
 for pop,mask in [('tissue',np.ones(len(dt),bool)),('Epi',dt['query']),('tumor_Epi',dt.tumor_epi),('P95',dt.MP10_P95)]:
  for p,g in dt[mask].groupby('program'):
   summary.append(dict(section_id=sec,population=pop,program=p,n=len(g),n_detected_median=g.n_detected_genes.median(),fraction_detecting_le2=(g.n_detected_genes<=2).mean(),fraction_detecting_zero=(g.n_detected_genes==0).mean(),top1_count_share_median=g.top1_count_share.median(),interpretation='Count contribution, not UCell rank contribution'))
 try:pa=pd.read_csv(O/'tables'/f'{sec}_matching.tsv.gz',sep='\t')
 except pd.errors.EmptyDataError:pa=pd.DataFrame()
 if 'matching' in pa:
  for k,g in pa.groupby('matching'):
   hi=d.set_index('barcode').loc[g.high_barcode];co=d.set_index('barcode').loc[g.control_barcode];t=d[d.tumor_epi]
   row=dict(section_id=sec,matching=k,n_pairs=len(g),median_match_distance=g.match_distance.median(),p90_match_distance=g.match_distance.quantile(.9),median_abs_front_difference=g.front_difference.abs().median())
   for col in ['mean__Epithelial','nCount_RNA']:
    sd=np.log1p(t[col]).std();row[col+'_SMD']=(np.log1p(hi[col]).mean()-np.log1p(co[col]).mean())/sd if sd>0 else np.nan
   balance.append(row)
 # Segment identity comes from the unchanged v3.4 table; no new boundary is generated.
 old=pd.read_parquet(F/'spots'/f'{sec}.parquet',columns=['barcode','reference_boundary_segment_id','reference_boundary_support','region_version','distance_unit']);a=d.merge(old,on='barcode',validate='1:1')
 for sid,g in a[a.tumor_epi&a.front_band.eq('front')].groupby('reference_boundary_segment_id',dropna=False):
  seg.append(dict(section_id=sec,reference_boundary_segment_id=sid,n_front=len(g),n_active=g.pEMT_high_front.sum(),active_fraction=g.pEMT_high_front.mean(),n_P95=g.MP10_P95.sum(),MP10_median=g[RAW].median(),pEMT_median=g[PE].median(),Flores_median=g.FLORES_CRISPR_PRO_INVASION_10.median()))
 c=pd.read_parquet(O/'centers'/f'{sec}.parquet')
 for rr in RINGS:
  n=c[rr+'_n_neighbors'];dist=c[rr+'_mean_distance_um'];ringcov.append(dict(section_id=sec,ring=rr,n_centers=len(c),n_covered=(n>0).sum(),median_neighbors=n[n>0].median(),min_neighbors=n.min(),median_mean_distance_um=dist.median(),q05_mean_distance_um=dist.quantile(.05),q95_mean_distance_um=dist.quantile(.95)))
 # Exact mapped-symbol availability for every member, retaining original symbol.
 genes=set(Path(r.genes_path).read_text().splitlines());det=pd.read_csv(O/'scores'/f'{sec}_gene_detection.tsv',sep='\t');detect=det.groupby('gene_symbol').n_detected_spots.max().to_dict()
 for row in G[G.gs_name.str.startswith(('PURAM_','FLORES_'))].itertuples():maprows.append(dict(section_id=sec,program=row.gs_name,original_symbol=row.original_symbol,mapped_symbol=row.gene_symbol,mapped_available=row.gene_symbol in genes,original_available=row.original_symbol in genes,n_detected_spots=detect.get(row.gene_symbol,0)))
for a,n in [(summary,'signature_detection_summary'),(balance,'matching_balance'),(seg,'active_front_segments'),(maprows,'literature_symbol_coverage'),(ringcov,'ring_coverage')]:save(a,O/'tables'/f'{n}.tsv')
# Concise module coverage, including explicit unavailability reasons.
a=pd.read_csv(O/'tables/effects.tsv.gz',sep='\t');c=a.groupby(['section_id','analysis','adjustment']).agg(n_programs=('program','nunique'),n_effects=('effect','count'),n_spatial_p=('p','count'),max_n=('n','max')).reset_index();c['status']=np.where(c.n_effects==0,'not_estimable',np.where(c.n_spatial_p==0,'effect_only_or_hierarchical_inference','estimated'));save(c,O/'tables/module_coverage.tsv')
print('diagnostics complete')
