from fu_common import *
F=yaml.safe_load((R/'frozen_config.yaml').read_text());d=pd.read_parquet(R/'tables/expanded_mixtures.parquet');group=['background','competitor','depth'];neg=d[d.partition.eq('train')&d.Treg_RNA_fraction.eq(0)];thresholds=[]
for key,v in neg.groupby(group):
 valid=v[v.T_eligible&v.multigene&v.core];thr=float(valid.Treg_main.quantile(.95,interpolation='higher'))if len(valid)else 0.;thresholds.append(dict(zip(group,key),eligible_negatives=len(valid),all_negatives=len(v),candidate_threshold=thr))
cut=max(q['candidate_threshold']for q in thresholds);put(pd.DataFrame(thresholds),'training_negative_thresholds.tsv')
# Freeze the single candidate threshold before reading heldout performance. No heldout retuning.
rule={'frozen_at':stamp(),'threshold':cut,'old_threshold':F['old_threshold'],'T_eligible_unchanged':True,'multigene_min':3,'core_min':1,'maxRank':1300,'signature_unchanged':True,'training_only':True,'selection':'maximum training negative conditional95th higher quantile across all background/competitor/depth strata','applicability_total_UMI':[500,5000]};rp=R/'candidate_evidence_rule.yaml'
if rp.exists():
 oldrule=yaml.safe_load(rp.read_text());rule['frozen_at']=oldrule['frozen_at'];assert rule==oldrule,'Candidate definition changed'
else:rp.write_text(yaml.safe_dump(rule,sort_keys=False))
d['high_candidate']=d.Treg_main>cut;d['joint_candidate']=d.T_eligible&d.multigene&d.core&d.high_candidate;d.to_parquet(R/'tables/mixture_rule_results.parquet',index=False)
rows=[];flows=[]
for key,v in d.groupby(['partition']+group+['Treg_T_RNA_fraction']):
 base=dict(zip(['partition']+group+['Treg_T_RNA_fraction'],key));base.update(n=len(v),Treg_total_RNA_fraction=v.Treg_RNA_fraction.iloc[0])
 for name in ['joint_old','joint_candidate','T_eligible','multigene','core','high_old','high_candidate','score_zero']:
  k=int(v[name].sum());lo,hi=wilson(k,len(v));rows.append(dict(**base,metric=name,count=k,rate=k/len(v),low=lo,high=hi,uncertainty='Wilson conditional simulations; shared finite reference cells, not patient CI'))
 vv=v.copy();checks=[('T_qualification',vv.T_eligible),('multigene',vv.multigene),('core',vv.core),('rank_window',~vv.score_zero),('old_high_score',vv.high_old)];remaining=np.ones(len(vv),bool)
 for name,mask in checks:
  fail=remaining&~mask.to_numpy();flows.append(dict(**base,stage=name,first_failed=int(fail.sum()),rule_failed_alone=int((~mask).sum()),passed_after_stage=int((remaining&mask.to_numpy()).sum())));remaining&=mask.to_numpy()
put(pd.DataFrame(rows),'mixture_per_condition_performance.tsv');put(pd.DataFrame(flows),'mixture_first_failure.tsv')
# Pooled heldout gain plus per-background FPR and positive-condition protection, frozen acceptance.
val=d[d.partition.ne('train')];vn=val[val.Treg_RNA_fraction.eq(0)];vp=val[val.Treg_RNA_fraction.gt(0)];nr=vn.groupby(group).joint_candidate.agg(['sum','count','mean']);nr['wilson_upper']=[wilson(k,n)[1]for k,n in zip(nr['sum'],nr['count'])];pr=vp.groupby(group+['Treg_T_RNA_fraction'])[['joint_old','joint_candidate']].mean();gain=float(vp.joint_candidate.mean()-vp.joint_old.mean());accepted=bool((nr['mean']<=.05).all()and(nr.wilson_upper<=.15).all()and gain>=.05 and((pr.joint_candidate-pr.joint_old)>=-.02).all());decision=dict(candidate_threshold=cut,adopt=accepted,heldout_old_positive_sensitivity=float(vp.joint_old.mean()),heldout_new_positive_sensitivity=float(vp.joint_candidate.mean()),gain=gain,max_heldout_negative_FPR=float(nr['mean'].max()),max_negative_Wilson_upper=float(nr.wilson_upper.max()),min_condition_sensitivity_change=float((pr.joint_candidate-pr.joint_old).min()),selection_never_used_HD_MP10=True,reason='all frozen acceptance criteria met'if accepted else 'one or more heldout FPR/uncertainty/gain/nonloss criteria failed; retain v1 rule',continuous_score_changed=False,analysis_population_changed=False)
(R/'evidence_decision.yaml').write_text(yaml.safe_dump(decision,sort_keys=False));put(nr.reset_index(),'heldout_negative_FPR.tsv');put(pr.reset_index(),'heldout_positive_detection.tsv')
# Original replay and actualHD full flow accounting.
flow=[];transition=[];coverage=[];examples=[]
for source,frame in [('original70_replay',pd.read_parquet(R/'tables/original70_replayed.parquet'))]:
 for key,v in frame.groupby(['competitor','depth','Treg_RNA_fraction']):
  remain=np.ones(len(v),bool)
  for name,mask in [('T_qualification',v.T_eligible),('multigene',v.multigene),('core',v.core),('rank_window',~v.score_zero),('old_high_score',v.high_old)]:
   flow.append(dict(section=source,competitor=key[0],depth=key[1],Treg_total_RNA=key[2],stage=name,n_total=len(v),first_failed=int((remain&~mask.to_numpy()).sum()),rule_failed_alone=int((~mask).sum()),passed_after_stage=int((remain&mask.to_numpy()).sum())));remain&=mask.to_numpy()
for s in SECTIONS:
 t=pd.read_parquet(V/s/'T32_regions.parquet');t['old_evidence']=t.evidence;t['calibration_in_range']=t.nCount.between(500,5000);t['candidate_joint']=t.measurement_eligible&(t.Treg_main_detected>=3)&(t.Treg_core_detected>=1)&(t.Treg_main>cut);t['followup_evidence']=t.old_evidence
 if accepted:
  mask=t.measurement_eligible&t.calibration_in_range;t.loc[mask,'followup_evidence']=np.where(t.loc[mask,'candidate_joint'],'Treg_like_calibrated','T_lineage_or_activated_calibrated')
 t['calibration_status']=np.where(~t.measurement_eligible,'measurement_insufficient',np.where(~t.calibration_in_range,'outside_reference_depth_range','applied_common_rule'if accepted else 'candidate_not_adopted'))
 t['label_change']=np.where(t.old_evidence.eq('Treg_like_multigene')==t.followup_evidence.str.startswith('Treg_like'),'unchanged',np.where(t.followup_evidence.str.startswith('Treg_like'),'increased','decreased'));t['formal_communication_ready']=False
 ranks=pd.read_parquet(V/s/'T32_measure.parquet',columns=['region_id']+['rank_'+g for g in SIG]).set_index('region_id');nwindow=(ranks.loc[t.region_id,['rank_'+g for g in SIG]].to_numpy()<1300).sum(1);zero=t.Treg_main.eq(0);det=t.Treg_main_detected>0;t['rank_window_detected_n']=nwindow;t['score_zero_reason']=np.where(~zero,'nonzero',np.where(det,'detected_but_outside_rank_window','no_signature_gene_detected'));assert (nwindow[zero]==0).all()
 cols=['region_id','section','x_um','y_um','component','nCount','nFeature','measurement_eligible','old_evidence','followup_evidence','calibration_status','label_change','candidate_joint','calibration_in_range','Treg_main','T_lineage','Tfh','Activated_T','Treg_main_detected','Treg_core_detected','Treg_main_UMI','rank_window_detected_n','score_zero_reason','formal_communication_ready'];t[cols].to_parquet(R/'tables'/f'{s}_evidence.parquet',index=False)
 remain=np.ones(len(t),bool)
 for name,mask in [('T_qualification',t.measurement_eligible),('multigene',t.Treg_main_detected>=3),('core',t.Treg_core_detected>=1),('rank_window',~zero),('old_high_score',t.Treg_main>F['old_threshold'])]:
  flow.append(dict(section=s,stage=name,n_total=len(t),first_failed=int((remain&~mask.to_numpy()).sum()),rule_failed_alone=int((~mask).sum()),passed_after_stage=int((remain&mask.to_numpy()).sum())));remain&=mask.to_numpy()
 transition.append(dict(section=s,all_T32=len(t),T_eligible=int(t.measurement_eligible.sum()),old_Treg_like=int(t.old_evidence.eq('Treg_like_multigene').sum()),new_Treg_like=int(t.followup_evidence.str.startswith('Treg_like').sum()),candidate_joint_all=int(t.candidate_joint.sum()),candidate_joint_in_calibration_range=int((t.candidate_joint&t.calibration_in_range).sum()),detected_but_zero_score=int((det&zero).sum()),eligible_detected_but_zero=int((det&zero&t.measurement_eligible).sum()),increased=int(t.label_change.eq('increased').sum()),decreased=int(t.label_change.eq('decreased').sum()),eligible_outside_range=int((t.measurement_eligible&~t.calibration_in_range).sum()),adopted=accepted))
 for q in t[det&zero&t.measurement_eligible].head(10).itertuples():
  row=dict(section=s,region_id=q.region_id,nFeature=q.nFeature,nCount=q.nCount,UCell=q.Treg_main)
  for g in SIG:row[g+'_raw']=getattr(q,'raw_'+g);row[g+'_rank']=float(ranks.loc[q.region_id,'rank_'+g])
  examples.append(row)
 n=pd.read_parquet(V/s/'T32_neighborhoods.parquet');n=n[n.source_rule.eq('exclude_shared')&n.association_eligible].merge(t[cols],left_on='T_region_id',right_on='region_id',suffixes=('','_followup'),validate='many_to_one')
 for rad,v in n.groupby('radius_um'):coverage.append(dict(section=s,radius_um=rad,continuous_analysis_centers=len(v),old_like_centers=int(v.old_evidence.eq('Treg_like_multigene').sum()),new_like_centers=int(v.followup_evidence.str.startswith('Treg_like').sum()),candidate_like_centers=int(v.candidate_joint.sum()),score_and_continuous_population_unchanged=True))
put(pd.DataFrame(flow),'screening_first_failure.tsv');put(pd.DataFrame(transition),'five_section_evidence_changes.tsv');put(pd.DataFrame(coverage),'evidence_neighborhood_coverage.tsv');put(pd.DataFrame(examples),'detected_zero_UCell_examples.tsv');print('EVIDENCE_DECISION',decision,flush=True)
