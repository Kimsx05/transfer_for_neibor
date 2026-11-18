from common import *
m=pd.read_parquet(F/'tables/expanded_mixtures.parquet');cfg=yaml.safe_load((V/'frozen_config.yaml').read_text());old=cfg['Treg_like_threshold_training_competitor_q95'];keys=['background','competitor','depth'];train=m[m.partition.eq('train')&m.Treg_RNA_fraction.eq(0)];rows=[];required=[]
for key,v in train.groupby(keys):
 n=int(v.T_eligible.sum());pre=v.T_eligible&v.multigene&v.core;scores=np.sort(v.loc[pre,'Treg_main'].to_numpy());k=int(np.floor(.05*n+1e-10));cut=0. if (scores>0).sum()<=k else float(scores[-k-1]);cut=max(0.,cut);required.append(cut)if n else None
 rows.append(dict(zip(keys,key),T_eligible_negatives=n,pre_gate_negatives=int(pre.sum()),allowed_FP=k,min_cutoff=cut if n else np.nan,empty_condition=n==0))
cut=max(required);rule=dict(cutoff=cut,comparison='score > cutoff',denominator='all T_eligible negatives within condition',selection='minimum nonnegative threshold satisfying every nonempty train condition full joint FPR<=.05',old_strict_threshold=old,previous_v2_candidate=float(yaml.safe_load((F/'evidence_decision.yaml').read_text())['candidate_threshold']),train_only=True,holdout_previously_seen=True,adopted_for_HD=False)
(R/'calibration_fix.yaml').write_text(yaml.safe_dump(rule,sort_keys=False));put(pd.DataFrame(rows),'calibration_training_constraints.tsv');m['fixed_joint']=m.T_eligible&m.multigene&m.core&m.Treg_main.gt(cut);perf=[]
for key,v in m.groupby(['partition']+keys+['Treg_T_RNA_fraction']):
 n=int(v.T_eligible.sum());truth=bool(v.Treg_RNA_fraction.iloc[0]>0);base=dict(zip(['partition']+keys+['Treg_T_RNA_fraction'],key));base.update(total_mixtures=len(v),T_eligible_denominator=n,Treg_total_RNA_fraction=v.Treg_RNA_fraction.iloc[0],metric='detection'if truth else'FPR',small_stratum=n<40)
 for rule_name,col in [('historical','joint_old'),('denominator_fixed','fixed_joint')]:
  k=int(v[col].sum());p=k/n if n else np.nan;z=1.96;den=1+z*z/n if n else np.nan;cen=(p+z*z/(2*n))/den if n else np.nan;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den if n else np.nan;perf.append(dict(**base,rule=rule_name,numerator=k,rate=p,low=cen-h,high=cen+h))
perf=pd.DataFrame(perf);assert (perf[perf.partition.eq('train')&perf.metric.eq('FPR')&perf.rule.eq('denominator_fixed')].rate.dropna()<=.05+1e-10).all();put(perf,'calibration_fixed_performance.tsv')
# Minimality: for each positive threshold the immediate previous attainable score must violate at least one constraint.
values=np.unique(np.r_[0.,train.loc[train.T_eligible&train.multigene&train.core,'Treg_main']]);prev=values[values<cut].max()if cut>0 else np.nan
viol=0
if np.isfinite(prev):
 for _,v in train.groupby(keys):
  n=v.T_eligible.sum()
  if n and(v.T_eligible&v.multigene&v.core&v.Treg_main.gt(prev)).sum()/n>.05+1e-10:viol+=1
assert cut==0 or viol>0
m[['profile_id','partition','background','competitor','depth','rep','Treg_RNA_fraction','Treg_T_RNA_fraction','T_eligible','multigene','core','Treg_main','joint_old','fixed_joint']].to_parquet(R/'tables/calibration_fixed_predictions.parquet',index=False)
pd.DataFrame([dict(status='recomputed_from_frozen_v2_cache_no_prior_fix_found',n=len(m),cutoff=cut,old_threshold=old,previous_v2_candidate=rule['previous_v2_candidate'],minimum_verified=True,previous_attainable=prev,violating_train_conditions_at_previous=viol,training_and_evaluation_denominator_same=True,no_empty_denominator_reported_zero=True,strict_labels_changed=False,FOXP3_definition_uses_cutoff=False,holdout_new_independent=False)]).to_csv(R/'calibration_fix_audit.tsv',sep='\t',index=False);print('FIXED_CUTOFF',cut,flush=True)
