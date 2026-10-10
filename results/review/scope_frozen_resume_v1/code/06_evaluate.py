from common import *
from scipy.stats import spearmanr,wilcoxon
import sys
parts=["DEV"] if "--DEV-only" in sys.argv else ["DEV","TEST"]
E=R/"06_evaluation"/("DEV_interim" if "--DEV-only" in sys.argv else "");E.mkdir(exist_ok=True)
cfg=json.loads((R/'CONFIG.yaml').read_text());frozen=json.loads((R/'06_evaluation/DEV_FROZEN.json').read_text());th=frozen['threshold_on_posterior_mean'];datasets=[]
for part in parts:
 truth=pd.read_csv(R/f'04_pseudospots/{part}/manifest.tsv.gz',sep='\t');truth['true_fraction_all']=truth.N_K32/truth.N_total;truth['true_fraction_epithelial']=truth.N_K32/truth.N_epithelial
 for branch in ['MAIN','NO_P07T']:
  O=R/f'05_deconvolution/{branch}_{part}';stat=json.loads((O/'MAPPING_DONE.json').read_text());d=truth.copy()
  for summ in ['mean','median','q05','q95']:
   p=pd.read_csv(O/f'abundance_{summ}.tsv.gz',sep='\t')
   if summ=='mean':p['estimated_total_cells']=p.drop(columns='pseudospot_id').sum(axis=1)
   p=p[['pseudospot_id','Epi_K32','Epi_other']+(['estimated_total_cells'] if summ=='mean' else [])].rename(columns={'Epi_K32':'K32_'+summ,'Epi_other':'other_'+summ});d=d.merge(p,on='pseudospot_id',validate='one_to_one')
  ratios=pd.read_csv(O/'posterior_ratios_and_epithelial.tsv.gz',sep='\t');d=d.merge(ratios,on='pseudospot_id',validate='one_to_one');d['reference']=branch;d['model_converged']=stat['converged'];d['detected']=d.K32_mean>th;d['true_positive']=d.N_K32>0;d['K32_error']=d.K32_mean-d.N_K32;d['other_error']=d.other_mean-d.N_other;d['epi_error']=d.total_epithelial_mean-d.N_epithelial;d['K32_abs_error']=abs(d.K32_error);datasets.append(d)
d=pd.concat(datasets,ignore_index=True);wt(d,E/'per_pseudospot_evaluation.tsv.gz')
depth=d[d.reference.eq('MAIN')].groupby(['split','depth_condition']).agg(n_spots=('pseudospot_id','size'),target_full_UMI=('target_full_UMI','first'),median_post_full_UMI=('post_full_UMI','median'),median_post_shared_UMI=('post_shared_UMI','median'),n_unattainable=('depth_target_unattainable','sum')).reset_index();wt(depth,E/'actual_pseudospot_depth_summary.tsv')
q=d[(d.reference=='MAIN')&(d.split=='DEV')];wt(q.groupby(['dataset','sample','scenario','depth_condition']).agg(n_spots=('pseudospot_id','size'),K32_MAE=('K32_error',lambda x:abs(x).mean()),K32_bias=('K32_error','mean'),epithelial_MAE=('epi_error',lambda x:abs(x).mean()),estimated_total_cells=('estimated_total_cells','mean')).reset_index(),E/'MAIN_DEV_source_diagnostic.tsv')
def wm(x,w):return float(np.sum(x*w)/w.sum())
def conditional_source_mean(z,mask,value):
 q=z.loc[mask].copy()
 if q.empty:return np.nan
 return float(q.assign(_metric_value=np.asarray(value)[np.asarray(mask)]).groupby('sample')._metric_value.mean().mean())
def metrics(z,equal_source=True):
 w=(1/z.groupby('sample')['sample'].transform('size')).to_numpy() if equal_source else np.ones(len(z));x=z.N_K32.to_numpy(float);y=z.K32_mean.to_numpy(float);err=y-x;xm=wm(x,w);ym=wm(y,w);den=wm((x-xm)**2,w);slope=wm((x-xm)*(y-ym),w)/den if den>0 else np.nan;pos=z.N_K32.to_numpy()>0;det=z.detected.to_numpy();neg=~pos
 o={'n_spots':len(z),'n_samples':z['sample'].nunique(),'n_source_groups':z.source_group.nunique(),'K32_MAE':wm(abs(err),w),'K32_MAE_positive':conditional_source_mean(z,pos,abs(err)),'K32_bias_positive':conditional_source_mean(z,pos,err),'K32_RMSE':np.sqrt(wm(err**2,w)),'K32_bias':wm(err,w),'K32_median_MAE':wm(abs(z.K32_median-x),w),'calibration_slope':slope,'calibration_intercept':ym-slope*xm,'K32_spearman_auxiliary':float(spearmanr(x,y).statistic) if len(np.unique(x))>1 and len(np.unique(y))>1 else np.nan,'Epi_other_MAE':wm(abs(z.other_error),w),'total_epithelial_MAE':wm(abs(z.epi_error),w),'total_cell_MAE':wm(abs(z.estimated_total_cells-z.N_total),w),'estimated_total_cell_mean':wm(z.estimated_total_cells,w),'simulated_positive_prevalence':wm(pos,w),'FPR':conditional_source_mean(z,neg,det),'recall':conditional_source_mean(z,pos,det),'precision':conditional_source_mean(z,det,pos),'ratio_all_MAE':wm(abs(z.K32_fraction_all_mean-z.true_fraction_all),w),'ratio_unstable_fraction':wm(z.epithelial_denominator_unstable,w),'posterior_90_interval_coverage_K32':wm((z.K32_q05<=x)&(z.K32_q95>=x),w),'posterior_90_interval_coverage_K32_positive':conditional_source_mean(z,pos,(z.K32_q05<=x)&(z.K32_q95>=x)),'mean_K32_90_interval_width':wm(z.K32_q95-z.K32_q05,w),'all_models_converged':bool(z.model_converged.all())}
 v=z.K32_fraction_epithelial_mean.notna().to_numpy();o['ratio_epithelial_MAE']=wm(abs(z.loc[v,'K32_fraction_epithelial_mean']-z.loc[v,'true_fraction_epithelial']),w[v]) if v.any() else np.nan;return o
specs={'overall':['reference','split'],'by_dataset':['reference','split','dataset'],'by_sample':['reference','split','dataset','sample','source_group'],'by_source_group':['reference','split','dataset','source_group'],'by_scenario_depth':['reference','split','scenario','depth_condition'],'by_dataset_scenario_depth':['reference','split','dataset','scenario','depth_condition'],'by_quantity_depth_epithelium':['reference','split','N_K32','depth_condition','N_epithelial'],'by_source_scenario_depth':['reference','split','dataset','sample','source_group','scenario','depth_condition']}
for name,keys in specs.items():
 out=[]
 for vals,z in d.groupby(keys,observed=True,dropna=False):out.append({**dict(zip(keys,vals)),**metrics(z)})
 wt(pd.DataFrame(out),E/f'metrics_{name}.tsv')
# Keep regular mixtures separate from pure-epithelial and intensity diagnostics.
regular=[]
for vals,z in d[d.scenario=='S1_regular'].groupby(['reference','split','N_K32','depth_condition'],observed=True):
 regular.append({**dict(zip(['reference','split','N_K32','depth_condition'],vals)),**metrics(z)})
wt(pd.DataFrame(regular),E/'regular_quantity_depth_sample_equal.tsv')
neg=d[d.N_K32==0];summ=[]
for keys,z in neg.groupby(['reference','split','dataset','sample','source_group','scenario','depth_condition']):
 summ.append({**dict(zip(['reference','split','dataset','sample','source_group','scenario','depth_condition'],keys)),'n_spots':len(z),'predicted_K32_mean':z.K32_mean.mean(),'predicted_K32_median':z.K32_mean.median(),'predicted_K32_q90':z.K32_mean.quantile(.9),'predicted_K32_q95':z.K32_mean.quantile(.95),'predicted_K32_q99':z.K32_mean.quantile(.99),'FPR':z.detected.mean()})
wt(pd.DataFrame(summ),E/'negative_abundance_distributions_by_source.tsv')
# Paired comparisons: same counts, truth, source, depth and spot IDs.
base=['pseudospot_id','split','dataset','sample','source_group','scenario','depth_condition','N_K32','N_epithelial'];cols=['K32_mean','K32_abs_error','detected','K32_fraction_epithelial_mean'];p=d[d.reference=='MAIN'][base+cols].merge(d[d.reference=='NO_P07T'][['pseudospot_id']+cols],on='pseudospot_id',suffixes=('_MAIN','_NO_P07T'),validate='one_to_one');p['delta_K32_mean']=p.K32_mean_NO_P07T-p.K32_mean_MAIN;p['delta_abs_error']=p.K32_abs_error_NO_P07T-p.K32_abs_error_MAIN;wt(p,E/'paired_reference_per_spot.tsv.gz');wt(p.groupby(['split','dataset','sample','source_group','scenario','depth_condition'],observed=True).agg(n_spots=('pseudospot_id','size'),delta_K32_mean=('delta_K32_mean','mean'),delta_MAE=('delta_abs_error','mean'),MAIN_detection=('detected_MAIN','mean'),NO_P07T_detection=('detected_NO_P07T','mean')).reset_index(),E/'paired_reference_by_source.tsv')
# State-intensity differences use paired backgrounds and report source-level effects, not pseudo-replicate p values.
i=d[d.scenario=='S4_intensity'];keys=['reference','split','sample','source_group','dataset','intensity_pair_id','depth_condition'];cols=['K32_mean','K32_mean_UCell','K32_mean_original_UMI','post_shared_UMI'];pairs=i[i.state_intensity=='low'][keys+cols].merge(i[i.state_intensity=='high'][keys+cols],on=keys,suffixes=('_low','_high'),validate='one_to_one');pairs['estimated_K32_high_minus_low']=pairs.K32_mean_high-pairs.K32_mean_low;wt(pairs,E/'intensity_paired_spots.tsv');isrc=pairs.groupby(['reference','split','dataset','sample','source_group','depth_condition'],observed=True).agg(n_pairs=('intensity_pair_id','size'),mean_estimated_high_minus_low=('estimated_K32_high_minus_low','mean'),mean_UCell_low=('K32_mean_UCell_low','mean'),mean_UCell_high=('K32_mean_UCell_high','mean'),mean_original_UMI_low=('K32_mean_original_UMI_low','mean'),mean_original_UMI_high=('K32_mean_original_UMI_high','mean')).reset_index();wt(isrc,E/'intensity_effect_by_source.tsv')
# No pseudo-replicate significance test; optional exploratory signed-rank at source-group level only.
tests=[]
for keys,z in isrc.groupby(['reference','split','depth_condition']):
 v=z.groupby('source_group').mean_estimated_high_minus_low.mean();pval=float(wilcoxon(v,method='exact').pvalue) if len(v)>=3 and (v!=0).all() and abs(v).nunique()==len(v) else np.nan;tests.append({**dict(zip(['reference','split','depth_condition'],keys)),'n_source_groups':len(v),'mean_effect':v.mean(),'p_two_sided_exploratory':pval,'interpretation':'few unverified source groups; not independent-patient confirmatory testing'})
wt(pd.DataFrame(tests),E/'intensity_source_group_exploratory.tsv');js(E/'EVALUATION_DONE.json',{'evaluation_scope':parts,'n_rows':len(d),'threshold':th,'TEST_models_converged':bool(d[d.split=='TEST'].model_converged.all()) if 'TEST' in parts else None,'no_truth_tuning':'threshold frozen before TEST mapping; analysis includes all generated scenarios'});say('EVALUATION_DONE',len(d))
