"""Finalize evidence only after all default mappings have actually exported posterior."""
from common import *
import subprocess,shutil

cfg=config(); names=['DEV']+cfg['sections']; states={}
for name in names:
 p=outmap(name)/'MAPPING_DONE.json'
 assert p.exists(),f'{name}: final posterior still pending'
 states[name]=json.loads(p.read_text())
 assert states[name]['actual_n_draws']==1000
 assert states[name]['stop_reason'] in ['NUMERICAL_STABLE','EPOCH_CAP','WALL_TIME']
assert (R/'00_manifest/DEFAULT_MAPPING_FINISHED.json').exists()
assert not json.loads((R/'00_manifest/DEFAULT_MAPPING_FINISHED.json').read_text())['errors']
launcher=json.loads((R/'00_manifest/RUN_ID.json').read_text())['launcher_pid']
ps=subprocess.run(['ps','-p',str(launcher),'-o','stat='],capture_output=True,text=True)
assert not ps.stdout.strip() or ps.stdout.strip().startswith('Z'),'Wait for scheduler final report refresh to exit before packaging'
active=R/'00_manifest/active_pids.tsv'
if active.exists() and not (R/'00_manifest/last_active_pids_at_scheduler_completion.tsv').exists():shutil.copy2(active,R/'00_manifest/last_active_pids_at_scheduler_completion.tsv')
wt(pd.DataFrame(columns=['name','pid','mode','start_epoch']),active)
val=json.loads((R/'00_manifest/VALIDATION_STATUS.json').read_text());assert val['n_failed']==0

metrics=pd.read_csv(R/'05_evaluation/DEV_recovery_metrics.tsv',sep='\t')
def metric(scope,target):
 return metrics[(metrics.scope==scope)&(metrics.level=='ALL')&(metrics.target==target)].iloc[0]
k=metric('A_C_unique','K32_count');e=metric('A_C_unique','epithelial_count');r=metric('A_C_unique','K32_fraction_epi');pos=metric('A_C_K32_positive','K32_count');inter=metric('A_C_ratio_interior','K32_fraction_epi')
assoc=pd.read_csv(R/'05_evaluation/spatial_associations.tsv',sep='\t');rows=[];paragraphs=[]
for sec in cfg['sections']:
 o=outmap(sec);st=states[sec];a=pd.read_csv(o/'abundance_mean.tsv.gz',sep='\t',index_col=0);d=pd.read_csv(o/'posterior_derived.tsv.gz',sep='\t');cv=pd.read_csv(o/'abundance_convergence.tsv',sep='\t');last=cv[cv.epoch==cv.epoch.max()];missing=int(d.I_state_mean.isna().sum());n=len(d)
 qs={c:float(a[c].median()) for c in [K,L,'vCAF']};rho={x.predictor:x.spearman for x in assoc[assoc.section==sec].itertuples()}
 fail=last.loc[~last.passed,'label'].tolist();rows.append({'section':sec,'epochs':st['epochs'],'numerical_status':st['numerical_status'],'actual_n_draws':st['actual_n_draws'],'n_spots':n,'n_ratio_NA':missing,'ratio_NA_fraction':missing/n,'median_K32_abundance':qs[K],'median_lower80_abundance':qs[L],'median_vCAF_abundance':qs['vCAF'],'median_I_state':float(d.I_state_mean.median()),'loss_relative_change':st['convergence']['loss_relative_change'],'failed_abundance_labels':';'.join(fail)})
 detail='；'.join(f"{x.label}: median={x.median_relative_change:.3%}, q95={x.q95_relative_change:.3%}" for x in last.itertuples())
 paragraphs.append(f"{sec}：保留{n:,}个spots，实际{st['epochs']:,} epochs，{st['actual_n_draws']:,}个丰度后验draws；{st['numerical_status']}，停止原因{st['stop_reason']}。末次丰度比较为{int(last.previous_epoch.iloc[0]):,}→{int(last.epoch.iloc[0]):,} epochs，{detail}。I_state因低分母不可稳定计算{missing}/{n}（{missing/n:.1%}），这些spots未从数量图删除。K32/lower80丰度中位数{qs[K]:.3f}/{qs[L]:.3f}。I_state与混合spot MP10、nCount、总上皮的Spearman分别为{rho['MP10_Top20_UCell']:.3f}、{rho['nCount_full']:.3f}、{rho['Epi_total_mean']:.3f}。这些描述不等于真值验证；未稳定时连图形模式本身也只是初版。")
 st['biological_evidence']='不可评价'
 st['biological_evidence_detail']='真实切片无状态/细胞数量真值；仅探索性空间表达一致性。'+('未满足预设丰度稳定性，不能据此作可靠状态定位或细胞定量。' if st['numerical_status']!='STABLE' else '数值稳定不代表生物学定量通过。')
 for fn in ['RESULT_STATUS.json','MAPPING_DONE.json']:js(o/fn,st)
 states[sec]=st
for name,st in states.items():
 p=outmap(name)/'progress.json';pg=json.loads(p.read_text()) if p.exists() else {};pg['last_training_stage']=pg.get('stage');pg.update({'stage':'FINAL_POSTERIOR_EXPORTED','actual_epochs':st['epochs'],'numerical_status':st['numerical_status'],'active_pid':None,'finalized_epoch_time':time.time()});js(p,pg)
wt(pd.DataFrame(rows),R/'05_evaluation/spatial_final_summary.tsv')

optional=R/'05_evaluation/PAIRED_BASELINE_STATUS.json'
if optional.exists():
 b=pd.read_csv(R/'05_evaluation/paired_baseline_summary.tsv',sep='\t');bb=b[(b.scope=='A_C_unique')&(b.level=='ALL')&(b.target=='K32')].iloc[0]
 baseline=f"已实际追加旧TRAIN-only S31对完全相同436个DEV bag的配对映射。A/C K32 MAE新{bb.new_MAE:.4f}、旧{bb.old_MAE:.4f}，配对差（新−旧）{bb.paired_MAE_delta_new_minus_old:.4f}。不同背景定义的表达重估与尺度变化均包含在这个整体比较内，不能单独归因为某一个基因或阈值。详情见paired_baseline_summary.tsv及DEV_OLD数值状态。"
else:
 assert not all(s['numerical_status']=='STABLE' for s in states.values()),'All defaults stable: reconsider authorized optional paired baseline before finalization'
 baseline='未追加旧S31同输入配对基线：附件要求默认三个mapping均已稳定后才启动可选基线，本轮该前提未满足。不是因DEV生物学效果不佳而跳过。旧v6数值仅为历史背景，不计算“标签修改带来的MAE改善”。第三张切片也未追加。'
js(R/'00_manifest/OPTIONAL_TASKS_STATUS.json',{'old_DEV_baseline':'COMPLETED' if optional.exists() else 'NOT_RUN_DEFAULT_STABILITY_PREREQUISITE_UNMET','third_section':'NOT_RUN_DEFAULT_DELIVERY_PRIORITY','reason':baseline})

text=f"""# 中文审阅结论

建议：先解决参考分配与数量校准问题，暂不扩大到全量切片，也不把当前高参考丰度解释为可靠的K32细胞数或全部MP10-high上皮细胞数。新参考存在多基因表达对比和DEV剂量方向响应，因此对“可供探索的K32-like参考状态”仅部分支持；数量恢复不支持可靠定量。真实空间结果的数值状态逐片列出，与生物学证据分开。

A. 范围、阈值和支持

TRAIN全部57,008个冻结上皮的float64 linear q80为0.22016683621566624。完整K32=3,282：3,176高于阈值（96.7703%），106不高于阈值仍保留；K32只覆盖全部24,242个高分上皮的13.1012%。这两个覆盖关系都不是伪spot检出性能。
TRAIN/DEV/TEST的K32人数为2,395/634/253，lower80为45,544/23,943/29,864，withheld高分非K32为9,069/4,385/7,612。全量lower80=99,351、withheld=21,066。旧范围外18,182另列outside_frozen_discovery_scope，既不标低状态也不进入本轮阴性真值。新NB使用126,791个TRAIN细胞、31类、11,230个同序G0基因；原Cycling来源合法冻结成员保留。25/8/8个保守来源组未独立核验为患者；DEV实有15个样本，不为匹配预计14而删样本。

B. Reference差异

NB在250 epochs达到预设loss稳定条件；K32/lower80 log1p Pearson=0.8545，旧K32/Epi_other=0.8837。相关性变化只作描述。NB signature log2FC（pseudocount=0.1）：ELOVL5=4.920、INSIG1=4.715、LDLR=3.882，另有ACSS2、HMGCR、FASN等同向。TRAIN两组各≥30细胞的6个样本、5个来源组、3个dataset中，ELOVL5/INSIG1方向6/6一致，LDLR为5/6。其余支持不足来源列NA；不是逐细胞显著性检验。
使用UCell定义lower80使程序差异部分成为构建预期。真正来源留出的证据来自下述DEV混合恢复，不能用TRAIN效应或同程序空间共定位代替。展示基因ELOVL5/INSIG1/LDLR在新空间预测前冻结。K32与旧K32 signature相关0.9941，表达总幅度比1.1111；lower80与旧other相关0.9981、幅度比0.9311，vCAF幅度比1.0299。背景变化未自动解决绝对数量尺度问题。

C. DEV两端回收、低状态背景与竞争分配

436个唯一同样本20细胞bag，8个来源组、4个dataset，full-RNA目标7005，2个技术重复；所有bag细胞包括withheld均贡献counts。DEV在9,000 epochs稳定，导出1,000 draws。以下主指标先组内样本等权、再来源组等权，不把bag当新增患者。
A/C共230个唯一bag：K32 MAE={k.MAE:.4f}、RMSE={k.RMSE:.4f}、bias={k.bias:.4f}，校准斜率={k.calibration_slope:.4f}、截距={k.calibration_intercept:.4f}。总体小MAE受大量低数量/零值影响：真实K32=12仅6个bag、3个样本/来源组，平均预测5.0957，MAE6.9043；不能包装成高数量准确回收。真实K32=0预测均值0.8582；低丰度也与零背景重叠，本轮未调二元阈值。
总上皮MAE={e.MAE:.4f}、bias={e.bias:.4f}、斜率={e.calibration_slope:.4f}。上皮内比例在200个可评价bag上的MAE={r.MAE:.4f}、bias={r.bias:.4f}，不支持精确比例解释。90%后验区间整体覆盖偏低；连续正值后验有零/一边界限制，但排除边界后问题仍在：阳性K32覆盖{pos.coverage90:.1%}（{int(pos.n_evaluable)}个bag），内部比例覆盖{inter.coverage90:.1%}（{int(inter.n_evaluable)}个bag）。这些区间是模型后验不确定度，不能当作已校准的频率学置信区间。
C面板纯非上皮bag：K32伪丰度均值0.1726、来源加权q95=0.5412。20个纯lower80上皮bag：K32伪丰度均值1.8115、q95=5.6740，同时总上皮只预测13.0711；vCAF真值0却预测0.7090。错误并非仅来自混入TME；全29类分配简表和逐来源极端值均保留，未删竞争类型。q95是模拟丰度分布的分位数，不是阳性概率或新检测阈值。
最强低背景误分配实际出现在InhouseData|T4380845和HRA003620|HRA003620_BC-30：真实K32均为0、纯lower80上皮20个，技术重复中K32最高预测分别5.6740和5.4155。A/C K32 MAE按dataset为CNP0000460 1.5994（1组）、HRA003620 0.8474（4组）、InhouseData 1.1677（2组）、PRJNA662018 0.7403（1组）；组数少且可构建条件不同，不作dataset优劣显著性排名。

D. 高分非K32泛化与目标解释

同具体背景可配对的A/B端点（高剂量减零剂量）：上皮4时K32-reference丰度增量A=1.4455、B=0.6060，I_state增量A=0.4324、B=0.1225，7个样本/6组/14配对重复；上皮12时丰度增量A=3.2796、B=1.4220，I_state增量A=0.4509、B=0.1380，3个样本/3组/6配对重复。B存在响应但弱于A；在exact K32目标下B真值始终为0，在广义状态目标下则是有限泛化证据，两种解释同时保留。
A/B高池的原始UCell分布、RNA贡献并未完全匹配，所以较强A响应不能单独证明严格K32特异性。逐来源剂量曲线及对真实高分比例/原始UCell均值的响应均提供。模型更接近带来源和RNA量依赖的K32-like参考分配，不是全部MP10活跃细胞的通用计数器；不能临时换成宽泛状态目标来掩盖exact K32数量误差。
{baseline}

E. 两张真实切片的实际状态与图形解释

两张均沿用已接受历史prior N_cells_per_location=20、detection_alpha=20；20不是实际核计数。保留原1017/2399个QC spots和完整31类。只有可靠物理坐标，没有可验证嵌入H&E，故用坐标图。原完整矩阵Top20均20/20可测，G0为19/20，HMGCS1不在G0；B17 KRT18不可测而非零表达。
{chr(10).join(paragraphs)}
按图中原始坐标方向观察：B1的K32绝对丰度在左侧/左上部较强，而I_state在右下部更高，说明“绝对分配多”与“上皮内部比例高”并非同一空间模式；该片未稳定，不能固定为可靠生物学区域。B17的高K32分配主要见于上部左侧连接带及部分下部边缘，与较高总上皮区域有重叠；右侧若干低上皮区域的比例为灰色NA。B17的K32与lower80分布大体共享上皮区域，而I_state对比明显较弱。原始EPCAM/KRT8与ELOVL5图可作表达参照，INSIG1/LDLR更为稀疏；无组织学依据给这些区域标肿瘤、正常、niche或某种细胞。上述图形描述与数值摘要共同提供，非独立空间验证。
两片共同丰度色标采用pooled q99，比例固定0–1，裁切数和低分母缺失见display_limits.tsv；完整值未截断。I_state是逐draw比例再汇总，不能用两个q05相除。空间强度不能直接贴成一种细胞，也不按scRNA q80给spot判阳性。混合spot MP10与模型I_state的相关/反相关均不能单独证明上皮内在状态；共享程序基因共定位仅表达一致性。空间观察必须结合DEV定量失败和当前数值稳定性解读。

F. 完成与未完成分别记录

成员、counts、G0、来源隔离、唯一bag、posterior比例与数值校验已执行，当前{val['n_checks']}项检查无失败。唯一新NB与DEV已稳定；两张真实切片都完成实际拟合、检查点、1,000 draws、重建QC和空间图，但“已实际训练/导出”不等于“数值稳定/生物学通过”。具体stop_reason、epochs、稳定变化量、模型路径和命令见MODEL_INDEX.tsv及各RESULT_STATUS.json。
如达到30,000 epochs仍未稳定，状态为PROVISIONAL_NOT_CONVERGED，不因loss平坦改写为稳定；这属于既定epoch上限，不是8小时墙钟耗尽。当前block命令遵守原30,000上限，不会暗中继续训练。需要另行决定后续策略时可从服务器保存模型/优化器恢复，不应重新初始化冒称续训。技术无效与有限但未稳定分开，本轮没有用伪造图补技术错误。

G. 下一步建议与边界

先检查数量尺度和参考竞争分配，利用本轮同来源bag、全类型误差和后验校准证据定位问题；在解决前不将这些图作可靠K32细胞密度或全MP10-high丰度，不扩大63张。可保留当前K32-like状态图作带清楚限制的探索材料。若真实片达到epoch上限仍未稳定，应先解决其丰度漂移，再讨论扩展。
K32来自既有全量发现集，DEV/TEST已有历史查看，本轮是冻结候选的探索性迁移而非全新独立发现或临床验证。伪spot不能完整模拟scRNA到真实Visium的平台差异；保守来源组不是已核验患者。P07T占TRAIN K32的1461/2395=61.0%，本轮未新增NO_P07T模型，不能排除该来源依赖。范围、数量和来源限制不因技术训练完成而消失。真实两片探索已按本轮新授权执行；没有自动扩展G1、多阈值、sample-balance、CNV、预后、CellChat或niche。
"""
(R/'07_review/INTERPRETATION.md').write_text(text)
stable=all(s['numerical_status']=='STABLE' for s in states.values())
js(R/'STATUS.json',{'execution_status':'COMPLETED' if stable else 'FINISHED_WITH_PROVISIONAL_SPATIAL_OUTPUTS','numerical_status':'ALL_STABLE' if stable else 'NB_AND_DEV_STABLE_SPATIAL_MODEL_SPECIFIC','biological_evidence':'部分支持','target_specific_evidence':{'exact_K32_count':'不支持可靠定量','K32_like_directional_state':'部分支持','all_MP10high_count':'不支持等价解释','real_spatial_validation':'不可评价'},'models':states,'active':{},'validation':val,'training_cap':30000,'training_deadline':cfg['training_deadline'],'delivery_deadline':cfg['delivery_deadline'],'updated':time.time(),'elapsed_hours':(time.time()-cfg['t0_epoch'])/3600,'report':str(R/'REPORT.html'),'review_zip':str(R/'07_review/review.zip'),'optional_tasks_note':baseline})
method=R/'METHODS.md';body=method.read_text().split('\n## 最终实际运行量\n')[0];body+='\n## 最终实际运行量\n\nNB 250 epochs；'+ '；'.join(f"{n} {s['epochs']} epochs、{s['actual_n_draws']} draws、{s['numerical_status']}（{s['stop_reason']}）" for n,s in states.items())+'。\n\n'+baseline+'\n';method.write_text(body)
say('EVIDENCE FINALIZED',stable)
