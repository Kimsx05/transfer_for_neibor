from shared import *
four=read(T/'four_way.tsv')
decomp=[];prof=[]
for (s,lv),g in four.groupby(['section_id','level']):
 for metric in ['pearson','resid_pearson','spearman','resid_spearman']:
  w=g.pivot(index='cell_type',columns='combination',values=metric)
  a,b,h,m=[w[c] for c in ['A','B','H05_c2l','MP10_transfer']]
  w['change_function_at_c2l']=h-b;w['change_function_at_transfer']=a-m
  w['change_TME_at_MP10']=m-b;w['change_TME_at_H05']=a-h
  w['function_mean_abs_change']=(abs(h-b)+abs(a-m))/2
  w['TME_mean_abs_change']=(abs(m-b)+abs(a-h))/2
  w['function_signed_component']=((h-b)+(a-m))/2;w['TME_signed_component']=((m-b)+(a-h))/2
  w['A_minus_B']=a-b;w['same_sign_AB']=(a*b>0).where(a.notna()&b.notna())
  assert np.allclose(w.A_minus_B,w.function_signed_component+w.TME_signed_component,equal_nan=True)
  decomp.extend(w.reset_index().assign(section_id=s,level=lv,metric=metric).to_dict('records'))
  z=w.dropna(subset=['A','B','H05_c2l','MP10_transfer'])
  prof.append(dict(section_id=s,level=lv,metric=metric,n_types=len(z),A_B_type_profile_spearman=corr(rankdata(z.A),rankdata(z.B)),same_sign=int(z.same_sign_AB.sum()),median_function_swap=z.function_mean_abs_change.median(),median_TME_swap=z.TME_mean_abs_change.median(),mean_function_swap=z.function_mean_abs_change.mean(),mean_TME_swap=z.TME_mean_abs_change.mean()))
save(decomp,'tables/endpoint_swaps.tsv.gz');save(prof,'tables/endpoint_swap_summary.tsv')
f=read(T/'four_way.tsv');fs=read(T/'functional_TME.tsv.gz')
ranks=[]
for t,g in f[f.level=='subtype'].groupby('cell_type'):
 w=g.pivot(index='section_id',columns='combination',values='resid_pearson').dropna()
 ranks.append(dict(cell_type=t,n_sections=len(w),A_positive=(w.A>0).sum(),B_positive=(w.B>0).sum(),A_B_positive_same=((w.A>0)&(w.B>0)).sum(),A_B_negative_same=((w.A<0)&(w.B<0)).sum(),A_min=w.A.min(),A_max=w.A.max(),B_min=w.B.min(),B_max=w.B.max(),A_median=w.A.median(),B_median=w.B.median()))
save(ranks,'tables/direction_by_subtype.tsv')


ref=read(O/'reference_manifest.tsv.gz')
save(ref.groupby(['model','label']).agg(n_reference=('barcode','size'),n_specimens=('sample_uid','nunique'),n_datasets=('dataset','nunique')).reset_index(),'tables/reference_label_support.tsv')
