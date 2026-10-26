exec(open('mp10_subtype_screen_round1/scripts/10_candidate_evidence.py').read().split('axes=pd.read_csv')[0])
r=pd.read_csv(T/'candidate_axis_evidence.tsv',sep='\t');d=pd.read_csv(T/'candidate_axis_dataset_support.tsv',sep='\t');n=pd.read_csv(T/'CellChat_epithelial_comparisons.tsv.gz',sep='\t')
choices=[('Treg','NECTIN2_TIGIT','priority_exploratory'),('FOLR2 TAM','MIF_CD74_CD44','priority_exploratory'),('C1QC TAM','HBEGF_EGFR_ERBB2','secondary_exploratory_discordant_expression')]
out=[]
for ty,lr,tier in choices:
 a=r[r.candidate_subtype.eq(ty)&r.LR_ID.eq(lr)].copy();assert len(a)==1;a['priority_tier']=tier;out.append(a)
 print('\n',ty,lr,tier)
 print(d[d.candidate_subtype.eq(ty)&d.LR_ID.eq(lr)][['dataset','network_difference','low_n_evaluable','high_n_evaluable','low_n_supported_patient_ids','high_n_supported_patient_ids','expression_product_difference','leave_one_recorded_patient_out_same_expression_direction_fraction']].to_string(index=False))
 print(n[n.LR_ID.eq(lr)&n.source.eq(a.source.iloc[0])&n.target.eq(a.target.iloc[0])&n.scope.eq('ALL')][['comparison','strength0','strength1','difference','status0','status1']].to_string(index=False))
pd.concat(out).to_csv(T/'selected_axes.tsv',sep='\t',index=False,na_rep='NA')
