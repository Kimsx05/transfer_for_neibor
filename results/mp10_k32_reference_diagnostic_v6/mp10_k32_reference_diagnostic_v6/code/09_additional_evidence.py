from common import *
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score,average_precision_score
E=R/'06_evaluation';D=pd.read_csv(E/'per_spot_all.tsv.gz',sep='\t');rows=[]
def weights(d):
 w=1/d.groupby('sample')['sample'].transform('size').to_numpy(float);w/=d.groupby('source_group')['sample'].transform('nunique').to_numpy();return w/w.sum()
def weighted(d,col):
 good=d[col].notna();w=weights(d)[good];return np.average(d.loc[good,col],weights=w) if good.any() else np.nan
for keys in [['model','dataset'],['model','scenario','depth_condition'],['model','N_K32','depth_condition','N_epithelial'],['model','dataset','scenario','depth_condition']]:
 for vals,d in D.groupby(keys,dropna=False):
  base=dict(zip(keys,vals));base.update(stratification='|'.join(keys[1:]),n_spots=len(d),n_samples=d['sample'].nunique(),n_groups=d.source_group.nunique())
  for ep in ['epithelial','total_cells','epithelial_fraction_all','K32','other','K32_fraction_all','K32_fraction_epithelial']:
   if ep+'_mean' not in d or not d[ep+'_mean'].notna().any():continue
   good=d[ep+'_mean'].notna();w=weights(d)[good];w/=w.sum();err=d.loc[good,ep+'_error'];rows.append({**base,'endpoint':ep,'MAE':float(np.sum(w*abs(err))),'bias':float(np.sum(w*err)),'RMSE':float(np.sqrt(np.sum(w*err**2))),'coverage90':weighted(d,ep+'_coverage90'),'width90':weighted(d,ep+'_width90'),'n_unstable_denominator':int(d.epithelial_denominator_unstable.sum()),'weighting':'sample_equal_within_group_then_group_equal'})
wt(pd.DataFrame(rows),E/'metrics_group_equal_all_strata.tsv')
# Exact control bag identities and counts pairing independent of predictions.
a=pd.read_csv(R/'04_paired_depth_control/TARGET_UMI/manifest.tsv.gz',sep='\t');b=pd.read_csv(R/'04_paired_depth_control/FIXED_CAPTURE/manifest.tsv.gz',sep='\t');j=a.merge(b,on=['bag_id','depth_condition'],suffixes=('_A','_B'),validate='one_to_one');checks=[]
for col in ['member_hash','sample','source_group','dataset','pre_full_UMI','N_total','N_K32','N_other','N_epithelial']:
 checks.append(dict(check='paired_'+col,passed=bool((j[col+'_A']==j[col+'_B']).all())))
checks.append(dict(check='all20distinct_cells',passed=bool((a.N_total==20).all())));assert all(x['passed'] for x in checks);wt(pd.DataFrame(checks),E/'control_pairing_validation.tsv');wt(j[['bag_id','depth_condition','member_hash_A','sample_A','source_group_A','pre_full_UMI_A','post_full_UMI_A','post_full_UMI_B','thinning_probability_A','thinning_probability_B']],E/'control_actual_depth_pairs.tsv')
# Frozen threshold within each source and condition; historical threshold is never re-fit.
cls=[]
for keys in [['model','scenario','depth_condition'],['model','N_K32','depth_condition','N_epithelial'],['model','dataset']]:
 for vals,d in D[D.labelset=='S31'].groupby(keys,dropna=False):
  base=dict(zip(keys,vals));r={**base,'stratification':'|'.join(keys[1:]),'n_samples':d['sample'].nunique(),'n_groups':d.source_group.nunique(),'n_spots':len(d),'simulated_positive_fraction':np.average(d.N_K32>0,weights=weights(d))}
  for name,mask,value in [('FPR',d.N_K32==0,d.K32_mean>8.367393),('recall',d.N_K32>0,d.K32_mean>8.367393),('precision',d.K32_mean>8.367393,d.N_K32>0)]:
   sub=d.loc[mask,['sample','source_group']].copy();sub['value']=value[mask];r[name]=sub.groupby(['source_group','sample']).value.mean().groupby('source_group').mean().mean() if len(sub) else np.nan;r[name+'_eligible_samples']=sub['sample'].nunique();r[name+'_eligible_groups']=sub.source_group.nunique()
  cls.append(r)
wt(pd.DataFrame(cls),E/'detection_group_equal_strata.tsv')
# Existing fullDEV state intensity diagnostic only; no new bags or redefined K32 core.
s4=D[D.scenario.str.contains('S4',na=False) & D.model.isin(['HIST_MAIN_DEV','HIST_NO_P07T_DEV','HIST_MAIN_TEST','HIST_NO_P07T_TEST','S30_FULL_DEV'])].copy();pairs=[]
if len(s4):
 for tag,d in s4.groupby('model'):
  lo=d[d.state_intensity=='low'];hi=d[d.state_intensity=='high'];keys=['intensity_pair_id','depth_condition'];p=lo.merge(hi,on=keys,suffixes=('_low','_high'),validate='one_to_one')
  for col in ['epithelial_mean','K32_mean','total_cells_mean','K32_fraction_epithelial_mean','K32_mean_UCell','K32_mean_original_UMI','post_full_UMI']:
   if col+'_low' in p:p['delta_high_minus_low__'+col]=p[col+'_high']-p[col+'_low']
  p['model']=tag;pairs.append(p)
 if pairs:
  p=pd.concat(pairs,ignore_index=True);wt(p,E/'intensity_paired_existing_spots.tsv.gz');cols=[x for x in p if x.startswith('delta_')];q=p.groupby(['model','dataset_low','source_group_low','sample_low'])[cols].mean().reset_index();wt(q,E/'intensity_by_sample.tsv')
# Pure sources and hard negatives never removed due to poor results.
negative=D[(D.N_K32==0)&(D.labelset=='S31')];wt(negative.groupby(['model','dataset','source_group','sample','scenario','depth_condition']).agg(n=('K32_mean','size'),mean=('K32_mean','mean'),median=('K32_mean','median'),q95=('K32_mean',lambda x:x.quantile(.95)),maximum=('K32_mean','max'),ratio_mean=('K32_fraction_epithelial_mean','mean'),FPR=('K32_mean',lambda x:(x>8.367393).mean())).reset_index(),E/'negative_distributions_by_source.tsv');wt(D[(D.N_epithelial==20)&D.N_K32.isin([0,20])],E/'pure_K32_other_diagnostic.tsv.gz')
# Low abundance, pure controls and source-specific negatives are all retained; no spot-wise p-values.
js(E/'ADDITIONAL_EVIDENCE_DONE.json',{'status':'COMPLETED','pairing_checks':len(checks),'all_passed':True,'S4':'Existing fullDEV/TEST only; no new intensity bags; paired source effect sizes, no spot-wise p-values','confidence_intervals':'None; all sample/group effects supplied to avoid false precision from technical repeats'})
say('ADDITIONAL_EVIDENCE_DONE')
