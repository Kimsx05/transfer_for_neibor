from common import *
import zipfile
import matplotlib.pyplot as plt

def main(out):
 P=out/'03_attribution_calibration';R=out/'reports';S=out/'05_evidence_synthesis'
 fit=json.loads((P/'calibration_completion.json').read_text());fid=json.loads((out/'01_score_fidelity/fidelity_summary.json').read_text())
 pairs=pd.read_csv(out/'04_fixed_niche_contrasts/paired_dataset_effects.tsv',sep='\t');main=pairs[pairs.analysis_set.eq('MAIN')&pairs.score.ne('spot_MP10_frozen')]
 state=pd.read_csv(P/'paired_state_recovery.tsv',sep='\t');bias=pd.read_csv(P/'composition_dependent_bias.tsv',sep='\t');loo=pd.read_csv(out/'02_gene_evidence/leave_one_gene_out_summary.tsv',sep='\t');v4=pd.read_csv(V4/'05_source_validation/source_evidence_all9.tsv',sep='\t')
 recovery=[]
 for (branch,bkg),g in state.groupby(['branch','background']):
  recovery.append(dict(branch=branch,background=bkg,n_pairs=len(g),mean_Delta_true=g.Delta_true.mean(),mean_Delta_hat=g.Delta_hat.mean(),rank_screened_sign_reversal_fraction=g.robust_rank_sign_reversal.mean(),truth_noise_screened_reversal_fraction=pd.to_numeric(g.truth_noise_screened_reversal,errors='coerce').mean() if 'truth_noise_screened_reversal' in g else np.nan,n_expected_high_realized_nonpositive=int(g.expected_high_but_realized_nonpositive.sum()),mean_abs_signal_change=g.abs_signal_change.mean()))
 write_tsv(out,P/'state_recovery_summary.tsv',recovery)
 write_tsv(out,P/'composition_bias_summary.tsv',bias.groupby(['branch','background','depth_band']).estimated_high_minus_low.agg(['count','mean','median','min','max']).reset_index())
 rows=[]
 for name in [f'niche{i}' for i in range(9)]+[f'niche{a}_minus_niche{b}' for a,b in PAIRS]:
  isp='minus' in name;g=main[main.pair.eq(name)&main.dataset_id.eq('POOLED_UNITS')];v=v4[v4.niche_name.eq(name)]
  obs=g[g.score.eq('spot_MP10_model_universe')].mean_effect.iloc[0] if isp else v.spot_within_unit_effect.iloc[0]
  epi=g[g.score.eq('epi_MP10_model_universe')].mean_effect.iloc[0] if isp else v.epi_within_unit_effect.iloc[0]
  prior=1 if name in ['niche1','niche6','niche1_minus_niche6'] else (2 if name in ['niche2','niche1_minus_niche2'] else 3)
  ll=loo[loo.pair.eq(name)]
  rows.append(dict(niche_or_pair=name,score_fidelity_in_epi=f"{fid['n_usage_matched']} cells;49 calculable sources; median rho={fid['donor_balanced_median_spearman']:.3f}; depth-sensitive",observed_direction=obs,model_epi_direction=epi,effect_definition='direct same-section A minus B' if isp else 'v4 versus-rest reference',paired_dataset_support='complete seven-dataset table retained',known_patient_vs_anonymous_unit_scope='verified patients bundled; policy-supported anonymous units distinguished',gene_level_coherence='gene_direction_unit_summary.tsv; all20 retained',single_gene_dominance=f'{int(ll.sign_changed.sum())}/40 layer-by-gene sign changes' if isp else 'contribution_concentration.tsv',fixed_weight_replay_result='1998 stress-test bags; platform mismatch recorded',known_mixture_recovery_result='540 bags;5 overlapping donor sources;composition shifts and state recovery errors',synthetic_refit_result=fit['end_to_end_status'],donor_overlap_status='all selected simulation cells IN_REFERENCE; no independent donor validation',depth_purity_overlap='fixed caliper matching; no relaxation; all9 retention reported',pathology_scope='v4 unknown unchanged; no new review task',major_failure_modes='allocation error;integer ties;low counts;reference overlap;model uncertainty not propagated',interpretation_allowed='SOURCE_ALLOCATION_SENSITIVE'+(';DATASET_HETEROGENEOUS' if name in ['niche7','niche1_minus_niche2','niche7_minus_niche6'] else ''),next_priority=prior,reason='Prioritize niche1 versus niche6 for limited further validation, not mechanism; matching and all20 deletions preserve direction' if prior==1 else 'Secondary/reference association; retain contrary and near-zero results'))
 ev=pd.DataFrame(rows);write_tsv(out,S/'evidence_matrix_all9.tsv',ev)
 stages=[
 ('M01_full_epi_fidelity','COMPLETED','123699 exact usage matches;49/57 sources calculable'),
 ('M01_official_AUCell','COMPLETED','integer observed/synthetic and fractional Epi;3 seeds;formal scores unchanged'),
 ('M02_gene_and_LOO','COMPLETED','all20 genes,all9 niches,three pairs;B5 excluded from full20'),
 ('M03A_replay','COMPLETED_STRESS_TEST',str(fit['replay_bags'])),
 ('M03B_known_mixture','COMPLETED_LIMITED_CALIBRATION',f"{fit['known_mixture_bags']} bags;5 sources;3 backgrounds;24 matched-model controls"),
 ('M03C_synthetic_refit',fit['end_to_end_status'],f"{fit['new_refit_bags']} bags;two fresh attempts on same background,one usable model"),
 ('M04_direct_comparisons','COMPLETED','3pairs;all7 datasets;min20;unit bootstrap2000;matching;8x8blocks;R1/R2 diagnostics'),
 ('M05_M06_delivery','COMPLETED_WITH_LIMITATIONS','all9+3pairs;deidentified review package')]
 write_tsv(out,R/'completion_matrix.tsv',[dict(module=a,status=b,detail=c) for a,b,c in stages])
 deviations=[
 ('first_synthetic_attempt','FAILED_TECHNICAL','5000 steps completed;logger=False omitted history;exception before model save. One retry with default logger and checkpoint-first order.'),
 ('second_gene_background','NOT_RUN','Two fresh attempts consumed on one fixed background;no third new fit.'),
 ('held_out_donors','NOT_ESTABLISHED','All selected donor cells match original reference;no independent donor validation.'),
 ('identity','LIMITED','Full epithelial orig.ident is source unit,not verified patient. Spatial identity policy unchanged.'),
 ('gradient','TWO_POINT_DIAGNOSTIC','P10/P90 positive central target fractions;no claim of global monotonicity.'),
 ('repeat_noise','DESCRIPTIVE','Two repeats in controlled panel;rank and truth-repeat SD flags separated;refit rep0 only.'),
 ('official_AUCell','IMPLEMENTATION_DIFFERENCE','Random integer ties differ from frozen average ties;no primary replacement.')]
 write_tsv(out,R/'deviations.tsv',[dict(item=a,status=b,reason=c) for a,b,c in deviations])
 lines=['# v5 实际完成结果','', '全部为结果知情后的探索/校准。v4、niche、Top20及真实模型保持只读。','',
 f"完整上皮：{fid['n_usage_matched']:,}细胞匹配原存储raw usage，固定{fid['n_genes']}基因背景；总体rho={fid['overall_spearman']:.4f}，来源等权中位rho={fid['donor_balanced_median_spearman']:.4f}。49个来源可计算，8个不足/常量；这些来源并非已核实患者。",'',
 f"M03A完成{fit['replay_bags']}个缓存bags回放。M03B完成{fit['known_mixture_bags']}个控制bags、5个reference重叠来源、3种背景、3个深度及24个模型匹配基线。M03C完成{fit['new_refit_bags']}bags，状态{fit['end_to_end_status']}；两次新拟合尝试仅一个可用模型。",'',
 '## 三个同切片直接比较','', '|对比|层|单位|均值差|探索95%区间|','|---|---|---:|---:|---|']
 for r in main[main.dataset_id.eq('POOLED_UNITS')].itertuples():lines.append(f'|{r.pair}|{r.score}|{r.n_units}|{r.mean_effect:+.5f}|[{r.ci_low:+.5f},{r.ci_high:+.5f}]|')
 lines+=['','## 每dataset效应','', '值为observed/model Epi，括号为有效单位；不足5单位只描述。','', '|dataset|1−2|1−6|7−6|','|---|---|---|---|']
 datasets=main.loc[~main.dataset_id.isin(['POOLED_UNITS','DATASET_EQUAL']),'dataset_id'].unique()
 for dataset in datasets:
  vals=[]
  for A,B in PAIRS:
   z=main[main.dataset_id.eq(dataset)&main.pair.eq(f'niche{A}_minus_niche{B}')];o=z[z.score.eq('spot_MP10_model_universe')].iloc[0];e=z[z.score.eq('epi_MP10_model_universe')].iloc[0];vals.append(f'{o.mean_effect:+.4f}/{e.mean_effect:+.4f} ({o.n_units})' if o.n_units else 'NA：无合格共同覆盖')
  lines.append('|'+dataset+'|'+'|'.join(vals)+'|')
  sub=main[main.dataset_id.eq(dataset)]
  sl=['# '+dataset,'','同切片直接比较；少于5单位不提供区间；反向来源保留。','','|pair|layer|units|effect|known/anonymous|','|---|---|---:|---:|---|']
  for r in sub.itertuples():sl.append(f'|{r.pair}|{r.score}|{r.n_units}|{r.mean_effect:+.5f}|{r.n_known_units}/{r.n_anonymous_units}|')
  safe(out,out/'04_fixed_niche_contrasts/by_dataset'/dataset/'dataset_review.md').write_text('\n'.join(sl))
 lines+=['','## 恢复校准','', '|分支|背景|状态对|平均真实Δ|平均估计Δ|rank筛查反转率|','|---|---|---:|---:|---:|---:|']
 for r in pd.DataFrame(recovery).itertuples():lines.append(f'|{r.branch}|{r.background}|{r.n_pairs}|{r.mean_Delta_true:.4f}|{r.mean_Delta_hat:.4f}|{r.rank_screened_sign_reversal_fraction:.1%}|')
 lines+=['','固定真上皮状态的组成偏差见 composition_bias_summary.tsv；按背景/深度分层的恢复误差见 calibration_summary_by_background.tsv。两者不同于真实niche效应。','',
 '优先保留niche1−niche6作为下一步受限验证对象；niche1−niche2次之；niche7−niche6当前覆盖/方向/删一稳健性较弱。来源校准限制所有机制解释。','',
 '完整九群与三个pair见 ../05_evidence_synthesis/evidence_matrix_all9.tsv；已知/匿名单位数与七dataset完整结果见 paired_dataset_effects.tsv。','',
 '续跑入口 code/resume_mp10_calibration_v5.sh；--verify-only 仅核验，不重新训练。']
 safe(out,R/'RUN_SUMMARY.md').write_text('\n'.join(lines))
 safe(out,S/'priority_decision.md').write_text('优先比较niche1与niche6：主分析、固定技术匹配及全部20项删一方向一致。这个优先级针对后续来源验证，不是阳性机制结论。niche1−niche2存在队列异质性且匹配后减弱；niche7−niche6的Epi差接近零、覆盖较少且删一敏感。模拟显示来源分配对组成及状态敏感，不能把相关或正式差异的区间替代恢复准确性。\n')
 limits='''# 尚未确定
- Top20与原usage关联随来源/深度变化；二者共享原基底，非独立生物真值。官方AUCell整数随机ties有实质差异，正式分数未替换。
- 回放存在生成平台/ambient与真实参数的失配；已知混合采用训练类型库大小校准的signature weights，不把cell fraction当RNA fraction。
- 五个模拟来源全部在reference中；无held-out donor。非上皮使用真实单细胞类型平均profiles，背景仅两个比例点，深度与供体范围有限。
- 合成拟合仅一个基因背景270bags。首次训练后日志接口失败，第二次修复；最多两次新拟合。目标收敛不能替代来源恢复准确性。
- 区间是单位bootstrap，不包含posterior归属不确定性；少量来源仅1–4个单位。匹配改变人群，不能排除全部混杂。
- 完整上皮的57个来源尚无完整患者合并表；空间身份沿用v4。病理仍未知，不重启133 ROI任务。
- MP10是转录程序，未测代谢通量或胆固醇/氧固醇。没有接触回避、免疫抑制、预后、因果或上皮特异性已证明的结论。
'''
 safe(out,R/'WHAT_IS_STILL_UNCERTAIN.md').write_text(limits)
 safe(out,R/'REVIEW_INDEX.md').write_text('# 入口\n\n[实际结果](RUN_SUMMARY.md)\n\n[校准](../03_attribution_calibration/calibration_report.md)\n\n[直接比较](../04_fixed_niche_contrasts/paired_dataset_effects.tsv)\n\n[九群+三个pair](../05_evidence_synthesis/evidence_matrix_all9.tsv)\n\n[优先级](../05_evidence_synthesis/priority_decision.md)\n\n[完成状态](completion_matrix.tsv)\n\n[限制](WHAT_IS_STILL_UNCERTAIN.md)\n')
 fig,ax=plt.subplots(figsize=(10,6));h=ev[['observed_direction','model_epi_direction']].to_numpy(float);lim=np.max(np.abs(h));im=ax.imshow(h,aspect='auto',cmap='coolwarm',vmin=-lim,vmax=lim);ax.set_yticks(range(len(ev)),ev.niche_or_pair);ax.set_xticks([0,1],['observed effect','model Epi effect']);ax.set_title('All9 versus-rest reference + three direct pairs');fig.colorbar(im,ax=ax);fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(safe(out,out/'06_figures'/('evidence_matrix.'+ext)),transparent=True,dpi=150)
 plt.close(fig)
 # Numerical recovery interpretation is part of reproducible delivery, not a pass based on convergence.
 abundance=pd.read_csv(P/'synthetic_abundance_recovery.tsv',sep='\t')
 am=abundance[['true_epi_RNA_fraction','estimated_epi_RNA_fraction','signature_weight_fraction_MAE']].agg(['mean','median','min','max']).reset_index(names='statistic')
 write_tsv(out,P/'abundance_recovery_summary.tsv',am)
 warning=f"合成端到端模型虽满足预设目标收敛规则，但真实上皮RNA占比平均{abundance.true_epi_RNA_fraction.mean():.3%}，估计平均{abundance.estimated_epi_RNA_fraction.mean():.3%}，存在明显低估。这是本合成设计中的来源恢复限制；不能外推为全部真实切片同幅度偏差，也不能把目标收敛或分数相关性称为表达归属验证通过。"
 for path in [R/'RUN_SUMMARY.md',R/'WHAT_IS_STILL_UNCERTAIN.md',S/'priority_decision.md']:
  path.write_text(path.read_text()+'\n\n'+warning+'\n')
 ev['synthetic_refit_result']=fit['end_to_end_status']+'; RNA fraction mean true=0.491648 vs estimated=0.179457; recovery limited'
 ev['major_failure_modes']=ev.major_failure_modes+'; synthetic Epi RNA underallocation'
 write_tsv(out,S/'evidence_matrix_all9.tsv',ev)
 cal=["# 表达归属恢复校准", "", warning, "",
 "M03A：1998个既有缓存bags，固定真实模型权重回放；生成平台/ambient失配使其属于压力测试，绝非端到端验证。M03B：540个控制bags，已知混合及固定signature；另24个模型匹配数值基线。M03C：270个bags仅向新c2l模型提供合成观测Y和技术batch，由模型重新估计权重，再归属上皮并对照真实上皮。", "",
 "新拟合两次尝试均在同一基因背景。首轮5000步后因未记录history且尚未存模型而失败；第二轮5000步保存checkpoint，修复history标量转换和AnnData斜杠列名导出兼容问题后恢复。最终目标相对变化0.003451，达到预设<0.01规则；未增加第三次拟合，未修改真实模型。", "",
 "五个上皮来源全部与训练reference重叠，故无独立donor泛化验证。固定13855基因；三档深度5852、12175、24337；GC B、SPP1 TAM、Neutrophil各两个组成比例点。只能解释这一有限面板，不能断言全范围单调。", "",
 "固定真上皮时，新拟合估计MP10仍随非上皮组成变化。例如GC B背景三个深度的平均组成差为−0.00535、−0.01235、+0.01389；Neutrophil为+0.02032、+0.00668、+0.00469；SPP1 TAM为+0.00387、−0.00796、−0.00778。方向依赖深度和donor，故不存在可统一减去的背景校正常数。", "",
 "固定组成的状态恢复：每背景15对，新拟合rank步长筛查后的方向反转率GC B 13.3%、Neutrophil 13.3%、SPP1 TAM 6.7%；再用两次真实重复的噪声筛查后分别0%、0%、6.7%。后者样本小且噪声估计仅两重复，不能解释成零误差。预设high/low池并不保证每个实现样本的真实MP10方向，已保留所有实现反向样本。", "",
 "GC B低深度固定组成状态的分数rho为−0.117，高深度为0.876；Neutrophil高深度rho为0.988仍伴平均偏差−0.0225及斜率约0.56。排序恢复、绝对分数校准、表达来源恢复是不同标准；高rho不能抵消上皮RNA低估。", "",
 "三分支使用不同样本范围：已知混合含两个重复，新拟合仅rep0；不可用未配对总体均值声称新拟合改善。参见state_recovery_summary.tsv、composition_bias_summary.tsv、calibration_summary_by_background.tsv、abundance_recovery_summary.tsv与逐基因gene_source_recovery.tsv。"]
 safe(out,P/'calibration_report.md').write_text('\n\n'.join(cal)+'\n')
 safe(out,out/'01_score_fidelity/score_fidelity_report.md').write_text(f"完整上皮{fid['n_usage_matched']:,}个细胞全部精确匹配原存储MP10 usage。总体Spearman={fid['overall_spearman']:.4f}；49个可计算来源的中位rho={fid['donor_balanced_median_spearman']:.4f}，范围0.1323–0.8483，均为正；另8来源不足或常量，不计为阳性。来源不是已验证患者。score与深度相关的来源中位rho=0.4065；技术调整后与usage的中位rank相关0.6360。共享基底的相关不能作为独立真值验证。官方AUCell三seed对照显示整数随机ties与历史平均ties有差异；所取fractional Epi样本一致，正式评分未替换。\n")
 safe(out,out/'02_gene_evidence/gene_level_report.md').write_text("保留全部20基因及九群：逐基因counts、检出率、归一化表达、Epi归属表达、rank和AUCell贡献均已输出。niche1−2与niche1−6的两层各20次删一均不改变总体方向；这不等于每基因独立显著。niche7−6的Epi层7/20、observed层1/20改变符号，但原始效应近零，不能据此认定某单基因驱动稳健生物差异。贡献集中程度及所有反向dataset保留。B5因不满足完整Top20排除于full20主分析并单列审计。\n")
 deviations.append(('saved_checkpoint_export','RECOVERED','Second fit checkpoint recovered after scalar history conversion and AnnData slash-column export compatibility fix; no additional optimization or fresh fit.'))
 write_tsv(out,R/'deviations.tsv',[dict(item=a,status=b,reason=c) for a,b,c in deviations])
 idcols={'section_id','section','analysis_unit_uid','analysis_unit','patient_uid','donor','unit','cell_id','bag_id','blueprint_uid','pair_id','group','fixed_epi_group','cell_ids','donor_cell_ids','path','local_path'}
 def scrub(df):
  df=df.copy()
  for c in df:
   if c in idcols or c.endswith('_uid') or 'cell_id' in c:df[c]=df[c].map(lambda v:'ID_'+hashlib.sha256(('v5_review_'+str(v)).encode()).hexdigest()[:12] if pd.notna(v) else v)
  return df
 zp=R/'review.zip'
 exclude={'original_usage.tsv','official_axes.tsv','official_synthetic_axes.tsv','frozen_reference.tsv','input_checksums.tsv','input_integrity.tsv','real_model_integrity.tsv','identity_and_independence.tsv','resolved_inputs.tsv','simulation_donors.tsv.gz'}
 with zipfile.ZipFile(safe(out,zp),'w',zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
  for p in out.rglob('*'):
   if not p.is_file() or p==zp:continue
   rel=p.relative_to(out);parts=rel.parts
   if any(x in parts for x in ['counts_cache','synthetic_refits','by_section','__pycache__','logs']) or p.name in exclude:continue
   if p.suffix=='.tsv' or p.name.endswith('.tsv.gz'):
    df=pd.read_csv(p,sep='\t');zz.writestr(str(rel).removesuffix('.gz'),scrub(df).to_csv(sep='\t',index=False))
   elif p.suffix in ['.md','.png','.pdf'] and parts[0]!='00_audit_and_spec':zz.write(p,str(rel))
   elif parts[0]=='code' and p.suffix in ['.py','.R','.sh']:zz.write(p,str(rel))
   elif p.name in ['fidelity_summary.json','calibration_completion.json','calibration_design.json','software_versions.json']:zz.write(p,str(rel))
  zz.writestr('SANITIZATION.md','Patient/source/cell/bag/section IDs are stable pseudonyms in tables. No identity key, full metadata, counts or models. Dataset labels retained.\n')
 write_json(out,R/'delivery_manifest.json',dict(zip_bytes=zp.stat().st_size,zip_sha256=sha(zp),deidentified=True,status='COMPLETED_WITH_LIMITATIONS',fresh_synthetic_fit_attempts=2,usable_synthetic_models=1 if fit['new_refit_bags'] else 0))
 log('reports and deidentified review package saved')
if __name__=='__main__':main(OUT)
