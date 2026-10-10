import bootstrap
from bootstrap import ROOT
import pandas as pd
import numpy as np
import json
def save(d,p):d.to_csv(p,sep='\t',index=False,na_rep='NA')
man=pd.read_csv(ROOT/'section_input_summary.tsv',sep='\t');dic=pd.read_csv(ROOT/'variable_dictionary.tsv',sep='\t');labels=dic.set_index('variable').source_label.to_dict()
params=json.loads((ROOT/'parameters.json').read_text());targetaudit=[]
for row in man.itertuples():
 folder=ROOT/'by_dataset'/row.dataset/row.section_id
 d=pd.read_csv(folder/'model_input.tsv',sep='\t');st=pd.read_csv(folder/'section_status.tsv',sep='\t')
 save(man[man.section_id==row.section_id],folder/'input_summary.tsv')
 bv=['MP10']+params['technical_variables']+params['tme_variables']
 b=d[d.analysis_eligible_legacy & d[bv].notna().all(axis=1)]
 records=[]
 for t in params['major_variables']+['MP10']:
  z=b if t=='MP10' else d[d[params['major_variables']].notna().all(axis=1)]
  var=z[t].var(ddof=1);reason=[]
  if len(z)<params['min_spots']:reason.append(f'n_complete_cases={len(z)} below {params["min_spots"]}')
  if not np.isfinite(var) or var==0:reason.append('target variance zero or undefined')
  records.append(dict(section_id=row.section_id,dataset=row.dataset,target=t,target_label=labels[t],n_tissue=len(d),n_legacy_eligible=int(d.analysis_eligible_legacy.sum()),n_complete=len(z),target_variance=var,n_target_zero=int(z[t].eq(0).sum()),n_unique_target=z[t].nunique(),status='SKIPPED_NOT_ESTIMABLE' if reason else 'PASS',reason='; '.join(reason) if reason else 'estimable target and sample size'))
 save(pd.DataFrame(records),folder/'target_variance_and_sample_audit.tsv');targetaudit.extend(records)
 f=folder/'mp10_baseline_vs_tme.tsv'
 if not f.exists() or f.stat().st_size<10:
  s=st[st.branch=='B'].iloc[0]
  save(pd.DataFrame([dict(model=tag,n_complete_case=len(b),n_predictions=0,target_variance=b.MP10.var(ddof=1),spatial_r2=np.nan,spatial_rmse=np.nan,spatial_mae=np.nan,spatial_q2_train_mean=np.nan,spatial_cv_complete=False,delta_r2_tme_minus_baseline=np.nan,delta_rmse_baseline_minus_tme=np.nan,importance_interpretable=False,status=s.status,reason=records[-1]['reason'] if s.status=='SKIPPED_NOT_ESTIMABLE' else s.reason) for tag in ['baseline','tme']]),f)
save(pd.DataFrame(targetaudit),ROOT/'target_variance_and_sample_audit.tsv')
files=['section_status.tsv','cellcell_importance.tsv','mp10_importance.tsv','mp10_baseline_vs_tme.tsv','internal_performance_unclipped.tsv','abundance_mp_technical_associations.tsv','mp10_spatial_fold_metrics.tsv','predictor_exclusion_audit.tsv','spatial_fold_audit.tsv']
all_tables={}
for name in files:
 parts=[]
 for row in man.itertuples():
  f=ROOT/'by_dataset'/row.dataset/row.section_id/name
  if not f.exists() or f.stat().st_size<3:continue
  d=pd.read_csv(f,sep='\t')
  for c in ['section_id','dataset','cohort','patient_uid','specimen']:d[c]=getattr(row,c)
  for c in ['Predictor','Target','predictor','target']:
   if c in d:d[c+'_label']=d[c].map(labels)
  parts.append(d)
 d=pd.concat(parts,ignore_index=True);save(d,ROOT/name);all_tables[name]=d
 for ds,g in d.groupby('dataset'):save(g,ROOT/'by_dataset'/ds/name)
for name,prefix in [('cellcell_importance.tsv','cellcell'),('mp10_importance.tsv','mp10')]:
 d=all_tables[name];agg=[]
 for (ds,t,pr),g in d.groupby(['dataset','Target','Predictor']):
  valid=g.standardized_importance.notna(); report=valid & g.importance_interpretable
  agg.append(dict(dataset=ds,Target=t,Predictor=pr,Target_label=labels.get(t,t),Predictor_label=labels.get(pr,pr),n_sections=len(g),n_estimable=int(valid.sum()),n_reportable=int(report.sum()),median_standardized_importance_all=g.loc[valid,'standardized_importance'].median(),median_standardized_importance_reportable=g.loc[report,'standardized_importance'].median(),median_raw_importance=g.raw_importance.median(),status='PASS' if report.any() else 'SKIPPED_NOT_ESTIMABLE'))
 a=pd.DataFrame(agg);save(a,ROOT/f'{prefix}_dataset_median_standardized_importance.tsv')
 for ds,g in a.groupby('dataset'):save(g,ROOT/'by_dataset'/ds/f'{prefix}_dataset_median_standardized_importance.tsv')
perf=all_tables['mp10_baseline_vs_tme.tsv'];agg=[]
for ds,g in perf[perf.model=='tme'].groupby('dataset'):
 b=perf[(perf.dataset==ds)&(perf.model=='baseline')]
 agg.append(dict(dataset=ds,n_sections=len(g),n_spatial_cv_complete=int(g.spatial_cv_complete.sum()),n_not_estimable=int((g.status=='SKIPPED_NOT_ESTIMABLE').sum()),n_positive_increment=int((g.delta_r2_tme_minus_baseline>0).sum()),n_reportable=int(g.importance_interpretable.sum()),median_baseline_spatial_r2=b.spatial_r2.median(),median_tme_spatial_r2=g.spatial_r2.median(),median_delta_r2=g.delta_r2_tme_minus_baseline.median(),median_delta_rmse=g.delta_rmse_baseline_minus_tme.median()))
save(pd.DataFrame(agg),ROOT/'dataset_performance_summary.tsv')
print(pd.DataFrame(agg).to_string(index=False))
