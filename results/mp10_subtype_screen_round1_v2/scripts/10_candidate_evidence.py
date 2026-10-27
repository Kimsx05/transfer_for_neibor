from common import *
d=read('B_Plasma_LR_contrasts.tsv.gz');e=read('candidate_sample_expression.tsv.gz');s=read('samples_v2.tsv');keys=['network_level','TME_type','direction','LR_id'];meta=['ligand','receptor','ligand_genes','receptor_genes','MP10_gene_overlap']
a=e[e.coverage_eligible & e.epi_delta.notna()].copy();a['expression_supported']=a.target_detection.ge(.1)&(a.epi_low_detection.ge(.1)|a.epi_high_detection.ge(.1));a['dataset']=a.dataset.astype(str)
pat=a.groupby(keys+['dataset','patient_uid_recorded'],dropna=False).agg(n_specimens=('sample_uid','size'),epi_delta=('epi_delta','mean'),aux_product_delta=('aux_product_delta','mean'),expression_supported=('expression_supported','max'),low_mean=('epi_low_mean','mean'),high_mean=('epi_high_mean','mean')).reset_index();save(pat,'candidate_patient_differences.tsv')
rows=[]
for kk,z in pat.groupby(keys+['dataset']):
 v=z.aux_product_delta.to_numpy();supported=z[z.expression_supported];w=supported.aux_product_delta.to_numpy();n=len(z);leave=[np.mean(np.delete(v,i)) for i in range(n)] if n>1 else []
 rows.append(dict(zip(keys+['dataset'],kk))|dict(n_patients=n,n_supported_patients=len(supported),n_positive=int((v>0).sum()),n_negative=int((v<0).sum()),n_zero=int((v==0).sum()),median_patient_delta=np.median(v),mean_patient_delta=np.mean(v),supported_median_delta=np.median(w) if len(w) else np.nan,leave_one_patient_min=min(leave) if leave else np.nan,leave_one_patient_max=max(leave) if leave else np.nan,leave_one_patient_direction_stable=bool(leave and (min(leave)>0 or max(leave)<0))))
q=pd.DataFrame(rows);net=d[d.comparison_type.eq('epi_state_within_specimen')&d.scope.ne('ALL')];q=q.merge(net[keys+['scope','raw_delta','status_low','status_high']],left_on=keys+['dataset'],right_on=keys+['scope'],how='left',validate='one_to_one');q['network_expression_concordant']=np.sign(q.raw_delta).eq(np.sign(q.supported_median_delta))&q.raw_delta.ne(0)&q.raw_delta.notna();q['coverage_concordant']=q.n_supported_patients.ge(2)&q.network_expression_concordant&(q.status_low.eq('detected')|q.status_high.eq('detected'));q['support_direction']=np.where(q.coverage_concordant,np.sign(q.raw_delta),0);save(q,'candidate_dataset_evidence.tsv')
rows=[]
for kk,z in q.groupby(keys):
 pos=z[z.support_direction.eq(1)];neg=z[z.support_direction.eq(-1)];direction=1 if len(pos)>len(neg) else -1 if len(neg)>len(pos) else 0;best=pos if direction==1 else neg if direction==-1 else z.iloc[:0]
 rows.append(dict(zip(keys,kk))|dict(n_datasets_evaluable=int(z.n_patients.ge(2).sum()),n_datasets_expression_supported=int(z.n_supported_patients.ge(2).sum()),n_concordant_positive_datasets=len(pos),n_concordant_negative_datasets=len(neg),n_datasets_same_direction=max(len(pos),len(neg)),consensus_direction=direction,n_supported_patients=int(z.n_supported_patients.sum()),n_patients_total=int(z.n_patients.sum()),n_direction_stable_datasets=int(best.leave_one_patient_direction_stable.sum()),support_datasets=';'.join(best.dataset),median_abs_paired_difference=float(np.nanmedian(abs(z.supported_median_delta))),evidence_tier='cannot_rank_no_dataset_with_two_patients' if not z.n_patients.ge(2).any() else 'multi_dataset_coverage_descriptive' if max(len(pos),len(neg))>=2 else 'single_dataset_or_discordant_descriptive'))
rank=pd.DataFrame(rows).merge(e[keys+meta].drop_duplicates(),on=keys);rank=rank.sort_values(['n_datasets_same_direction','n_direction_stable_datasets','n_supported_patients','median_abs_paired_difference'],ascending=False);save(rank,'B_Plasma_all_candidate_ranking.tsv')
selected=[]
# Prespecified coverage-based slots: major B/Plasma, each communication direction; then B subtypes where evaluable.
for lev,targets in [('major',['B','Plasma']),('fine',['GC B','Memory B','Naive B','IFN-stimulated B'])]:
 for target in targets:
  for direction in ['epi_to_TME','TME_to_epi']:
   z=rank[rank.network_level.eq(lev)&rank.TME_type.eq(target)&rank.direction.eq(direction)&rank.n_supported_patients.ge(2)&rank.n_datasets_evaluable.ge(1)]
   if len(z):selected.append(z.iloc[0])
sel=pd.DataFrame(selected).drop_duplicates(keys).reset_index(drop=True);sel['axis_id']=['AX%02d'%(i+1) for i in range(len(sel))];save(sel,'selected_review_axes.tsv')
# Primary and secondary are distinct contrasts. Compare direction descriptively, retain every dataset and NA.
context=d.merge(sel[keys+['axis_id']],on=keys);save(context,'selected_axes_all_network_context.tsv')
# State networks contain all common epi specimens; actual matched target subset is smaller. Quantify this explicitly.
con=read('network_cell_contributions.tsv.gz');con=con[con.network_id.str.contains('epi_state_within_specimen')];rows=[]
for ar in sel.itertuples():
 for ds in sorted(s.dataset.unique()):
  eligible=e[(e.network_level.eq(ar.network_level))&(e.TME_type.eq(ar.TME_type))&(e.direction.eq(ar.direction))&(e.LR_id.eq(ar.LR_id))&e.dataset.eq(ds)&e.coverage_eligible];ids=set(eligible.sample_uid);network=s[s.dataset.eq(ds)&s.epi_state_eligible];rows.append(dict(axis_id=ar.axis_id,network_level=ar.network_level,TME_type=ar.TME_type,dataset=ds,n_state_network_specimens=len(network),n_matched_target_specimens=len(ids),matched_specimen_fraction=len(ids)/len(network) if len(network) else np.nan,matched_patient_ids=eligible.patient_uid_recorded.nunique()))
save(rows,'candidate_matched_coverage.tsv');print(sel[['axis_id']+keys+['n_datasets_same_direction','n_supported_patients','support_datasets']].to_string(index=False))
