from common import *
import subprocess,sys,html,shutil,zipfile
mode=sys.argv[1] if len(sys.argv)>1 else 'incremental';cfg=config();q=json.loads((R/'00_manifest/Q80.json').read_text());models=['GSE171351_B1','DEV','GSE238145_B17']+(['DEV_OLD'] if (outmap('DEV_OLD')/'RESULT_STATUS.json').exists() else []);index=[]
for name in models:
 o=outmap(name);st=json.loads((o/'RESULT_STATUS.json').read_text()) if (o/'RESULT_STATUS.json').exists() else {};conv=json.loads((o/'convergence.json').read_text()) if (o/'convergence.json').exists() else {};index.append({'model':name,'path':str(o),'execution_status':st.get('execution_status','RUNNING' if (o/'model/model.pt').exists() else 'PENDING'),'numerical_status':st.get('numerical_status','NOT_YET_EVALUABLE'),'biological_evidence':st.get('biological_evidence','EVALUATION_PENDING'),'epochs':st.get('epochs',conv.get('epochs',0)),'actual_n_draws':st.get('actual_n_draws',0),'n_spots':st.get('n_spots',None),'resume_command':f'bash "{R}/code/python.sh" "{R}/code/04_mapping.py" {name} block'})
nb=json.loads((R/'01_reference/REFERENCE_DONE.json').read_text()) if (R/'01_reference/REFERENCE_DONE.json').exists() else {};index.insert(0,{'model':'NB_LOWER80','path':str(R/'01_reference'),'execution_status':'COMPLETED' if nb else 'RUNNING','numerical_status':nb.get('status','TRAINING'),'biological_evidence':nb.get('biological_evidence','DESCRIPTIVE_ONLY'),'epochs':nb.get('epochs',json.loads((R/'01_reference/progress.json').read_text()).get('epochs') if (R/'01_reference/progress.json').exists() else 0),'actual_n_draws':1000 if nb else 0,'n_spots':None,'resume_command':f'bash "{R}/code/python.sh" "{R}/code/01_reference.py"'});idx=pd.DataFrame(index);wt(idx,R/'MODEL_INDEX.tsv')
# Regenerate evaluation/figures only after a new result status or final request.
finger=key([str(x) for x in index]);old=(R/'07_review/.render_fingerprint').read_text() if (R/'07_review/.render_fingerprint').exists() else ''
if finger!=old or mode=='final':
 for script in ['06_evaluate.py','11_technical_summary.py','07_figures.py']:
  with open(R/'logs/report_render.log','a') as log:subprocess.run(['bash',str(R/'code/python.sh'),str(R/'code'/script)],stdout=log,stderr=subprocess.STDOUT,check=True)
 (R/'07_review/.render_fingerprint').write_text(finger)
def table(path,mask=None):
 p=R/path
 if not p.exists():return '<p>尚未生成，随模型进度更新。</p>'
 try:d=pd.read_csv(p,sep='\t')
 except pd.errors.EmptyDataError:return '<p>尚无可评价结果。</p>'
 d=mask(d) if mask else d
 return d.to_html(index=False,escape=True,float_format=lambda x:f'{x:.4g}',na_rep='NA',border=0)
def image(path,title):return '<figure><img src="'+path+'" loading="lazy"><figcaption>'+html.escape(title)+'</figcaption></figure>' if (R/path).exists() else '<p>'+html.escape(title)+'：等待模型摘要。</p>'
sections=[]
sections.append('<h2>A. 标签与来源范围</h2><p>TRAIN 全部 57,008 个冻结上皮的 q80 = <b>'+format(q['threshold'],'.17g')+'</b>（float64，linear / R type7）。K32 定义为冻结 FULL_Z_r1.00 的 C18，完整 3,282 个保留，其中高于阈值 3,176 个、低于或等于阈值 106 个；K32 内高分占比 '+f"{q['K32_high_fraction']:.2%}"+'，对全部高分上皮的覆盖率 '+f"{q['K32_coverage_all_high']:.2%}"+'。后者不是伪 spot 检出率。</p><p>新 NB：126,791 个 TRAIN 细胞，31 类，固定 G0 11,230 基因。非 K32 高分细胞从 reference 排除，保留用于 DEV 泛化。18,182 个旧范围外上皮仍标 outside_frozen_discovery_scope，不标作低状态或阴性。原 Cycling 来源合法成员保留。沿用 25/8/8 个保守来源组，不是已核验患者独立验证。</p>'+table('00_manifest/support_by_split_label_v7.tsv',lambda d:d[d.label_v7.isin([K,L,H])]))
sections.append('<h2>B. Reference 表达特征</h2><p>所有 31 类均用本轮 TRAIN 重新拟合，不拼接历史全数据 signature。展示基因 ELOVL5、INSIG1、LDLR 已在空间预测前按 TRAIN 支持与效应固定。按 UCell 选择 lower80 后程序基因差异部分属于设计预期，不独立证明生物学可区分性。计数均值 log2FC 使用 pseudocount=0.1，无逐细胞 DEG 显著性检验。</p>'+image('06_figures/reference_MP10_gene_expression.png','新参考中的冻结 MP10 基因表达')+image('06_figures/reference_K32_lower80_scatter.png','K32 与 lower80 的全部 G0 基因散点')+table('01_reference/old_new_signature_amplitudes.tsv',lambda d:d[d.new_label.isin([K,L,'vCAF'])]))
sections.append('<h2>C. DEV 已知数量恢复与低背景</h2><p>436 个唯一 bag，15 个实际 sample、8 个来源组；实际来源人数按清单核验，未为了附件预计的14而删来源。每 bag 同一样本20个不同细胞；A/B零剂量与C重叠仅映射一次。full RNA 锚点7005，整数二项稀释。所有20个细胞（含 withheld-high）贡献 counts；31列真值和可小于20。A/C用于 exact K32 回收；N_epi=0的比例是NA。</p>'+table('05_evaluation/DEV_recovery_metrics.tsv',lambda d:d[(d.scope=='A_C_unique')&(d.level=='ALL')])+table('05_evaluation/DEV_recovery_metrics.tsv',lambda d:d[(d.level=='ALL')&(((d.scope=='A_C_K32_positive')&(d.target=='K32_count'))|((d.scope=='A_C_ratio_interior')&(d.target=='K32_fraction_epi')))])+'<p>区间覆盖另列阳性K32与内部比例：连续正丰度后验对精确0、比例对精确0/1的覆盖存在边界限制；不以q05&gt;0定义检出。</p>'+image('06_figures/DEV/DEV/ALL/quantity_calibration.png','A/C：数量与比例恢复')+image('06_figures/DEV/DEV/ALL/C_low_background.png','低状态背景和纯非上皮的 K32-anchor 分配')+table('05_evaluation/DEV_low_background_distribution.tsv',lambda d:d[d.level=='n_epi_true'])+table('05_evaluation/DEV_29type_allocation.tsv',lambda d:d[d.label=='vCAF']))
sections.append('<h2>D. 高分非 K32 泛化与来源依赖</h2><p>B 的 exact K32 真值始终为0；对高参考的响应在广义 MP10 状态目标下另作泛化描述。二者并列，不因指标好坏切换目标。按上皮量4/12、同一具体背景比较剂量响应，完整逐来源结果见05_evaluation。主汇总先组内样本等权再来源组等权；技术重复不是新增患者。不使用二元检测阈值，不沿用8.367393或q05&gt;0。</p>'+image('06_figures/DEV/DEV/ALL/A_B_dose_response.png','A/B 相同背景中的剂量响应')+image('06_figures/DEV/DEV/ALL/A_B_true_state_response.png','对真实高分比例与原始 UCell 均值的响应')+table('05_evaluation/DEV_matched_endpoint_groupweighted.tsv')+table('05_evaluation/DEV_recovery_metrics.tsv',lambda d:d[(d.scope=='B')&(d.level=='ALL')&(d.target.isin(['K32_count','epithelial_count']))]))
sections.append('<h2>E. 两张真实切片探索</h2><p>切片由既定 QC/规模/assay 选择。两张均沿用已接受历史 prior：N_cells_per_location=20、detection_alpha=20；这是模型先验，不是真实核计数。所有原 QC spots 和完整31类保留，没有按预测上皮过滤。无可靠嵌入图像，采用原坐标图。MP10在完整实测基因宇宙中评分，B1精确复用原始Top20分数，B17按同一既定参数补算。</p>'+table('00_manifest/SELECTED_SECTIONS.tsv',lambda d:d.drop(columns=['counts_path']))+table('MODEL_INDEX.tsv',lambda d:d[d.model.isin(cfg['sections'])].drop(columns=['resume_command','path'])))
for sec in cfg['sections']:
 pref='06_figures/by_dataset/'+sec.split('_')[0]+'/'+sec+'/';sections.append('<h3>'+sec+'</h3>'+image(pref+'mapping_overview.png','共同 abundance 色标；I_state固定0–1；未稳定结果显著标记')+image(pref+'raw_expression_overview.png','完整原始表达、MP10连续分数和预先固定展示基因')+image(pref+'state_expression_associations.png','I_state与表达/UMI/上皮数量的描述性关系'))
sections.append(table('05_evaluation/spatial_associations.tsv')+'<p>MP10是混合 spot 转录组分数，不能自动解释为上皮内在代谢活性；与相同程序基因共定位仅是表达一致性，不是统计独立验证。q05不是阳性概率。低分母比例缺失标灰，数量估计仍保留。</p>')
sections.append('<h2>F. 运行状态、路径与续跑</h2>'+idx.to_html(index=False,escape=True,na_rep='NA',border=0)+'<p>loss检查：相邻500 epochs均值变化&lt;0.005；丰度检查：至少相隔2000实际epochs的256抽样均值，K32/lower80/总上皮相对变化中位数&lt;0.05且q95&lt;0.20；至少5000 epochs。时间限制与数值稳定分开。参考冻结后所有映射共享相同signature hash。实际运行配置、PID、种子与哈希见00_manifest。</p>')
sections.append(table('05_evaluation/model_technical_summary.tsv')+'<p>RNA重建质量与真实细胞数量恢复是不同评价对象；良好的表达重建不能抵消DEV数量偏差。达到epoch上限的模型保留检查点，当前block命令不会突破冻结上限。</p>')
for name in models:
 sections.append(image('06_figures/numerical/'+name+'_abundance_stability.png',name+'：预先固定的丰度稳定标准；平坦loss不等于丰度稳定'))
if (R/'05_evaluation/PAIRED_BASELINE_STATUS.json').exists():
 sections.append('<h3>同一DEV counts的旧S31配对基线</h3>'+table('05_evaluation/paired_baseline_summary.tsv',lambda d:(d[(d.scope=='A_C_unique')&(d.level=='ALL')]))+'<p>两次映射共享同一436个bag、counts、spot顺序与随机种子；这才是本轮可直接配对的新旧比较。不是相减历史v6指标。</p>')
interpret=(R/'07_review/INTERPRETATION.md').read_text() if (R/'07_review/INTERPRETATION.md').exists() else '模型仍在运行，生物学判断待完整数值输出后填写；原始空间表达图已经可查看。'
sections.append('<h2>G. 综合结论与下一步范围</h2><pre>'+html.escape(interpret)+'</pre><p>本轮K32 TRAIN中P07T贡献1461/2395（61.0%），没有新增去P07T参考，因此不能单独排除P07T依赖。本轮为既有冻结候选的探索性参考迁移，DEV/TEST此前已被查看；不是从头独立发现或临床验证。伪spot不能完整模拟scRNA与Visium平台差异。K32丰度不等于所有MP10-high细胞丰度。本轮明确授权两张真实空间探索，无全量63张扩展、G1、多阈值、CNV、CellChat或niche。</p>')
hero='<h2>核心结果总览</h2>'+idx[['model','execution_status','numerical_status','biological_evidence','epochs','actual_n_draws']].to_html(index=False,na_rep='NA',border=0)
if (R/'07_review/INTERPRETATION.md').exists():
 hero='<p style="background:#f5f7fa;padding:16px"><b>审阅建议：先解决参考分配和数量校准问题。</b>仅部分支持K32-like状态的方向响应；不支持把当前丰度作为可靠K32细胞计数或全部MP10-high细胞计数。真实切片的未稳定模型仅作探索初版，具体证据与限制见A–G。</p>'+hero
for sec in cfg['sections']:
 pref='06_figures/by_dataset/'+sec.split('_')[0]+'/'+sec+'/'
 hero+=image(pref+('mapping_overview.png' if (R/(pref+'mapping_overview.png')).exists() else 'raw_expression_overview.png'),sec+'：实际可用空间图；模型状态以上表为准')
hero+=table('05_evaluation/DEV_recovery_metrics.tsv',lambda d:d[(d.scope=='A_C_unique')&(d.level=='ALL')&(d.target.isin(['K32_count','epithelial_count','K32_fraction_epi']))])
body='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>K32 lower80 v7 审阅报告</title><style>body{font:16px/1.65 system-ui;margin:30px auto;max-width:1250px;padding:0 24px;color:#24313c}h1,h2{color:#173B75}table{border-collapse:collapse;font-size:12px;display:block;overflow-x:auto}td,th{padding:5px 8px;border-bottom:1px solid #ddd;text-align:left}img{max-width:100%;height:auto}figure{margin:20px 0}figcaption{font-size:13px;color:#555}pre{white-space:pre-wrap;font:inherit;background:#f5f7fa;padding:16px}</style><h1>K32 lower80 v7：参考、留出混合与真实空间探索</h1><p>更新时间：'+datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')+'；'+('最终审阅版' if mode=='final' else '增量运行版')+'。目标交付北京时间08:16:54；训练软截止07:16:54。</p>'+hero+''.join(sections)+'</html>'; (R/'REPORT.html').write_text(body)
(R/'README.md').write_text('# K32 lower80 spatial v7\n\n打开 REPORT.html 查看中文总览。轻量上传目录：07_review；最终包：07_review/review.zip。服务器完整模型路径及续跑命令见 MODEL_INDEX.tsv。\n\n本轮已授权：新31类NB参考、精简DEV、两张真实切片；旧v5/v6只读。本轮训练完成与数值稳定、生物学证据分开记录。\n')
with (R/'README.md').open('a') as f:
 f.write('\n## 实际模型状态\n\n'+ '\n'.join(f"- {x.model}: {x.epochs} epochs；{x.numerical_status}；{x.actual_n_draws} draws；生物学证据：{x.biological_evidence}。" for x in idx.itertuples())+'\n\n## 阅读入口\n\n- [中文总览](REPORT.html)\n- [实际方法](METHODS.md)\n- [模型路径与命令](MODEL_INDEX.tsv)\n- 核心已知真值结果：05_evaluation/DEV_recovery_metrics.tsv；逐来源、低背景及完整29类分配表在同目录。\n- 空间PDF/300dpi PNG、英文图注和plot_data：06_figures/by_dataset。\n- 实际耗时、稳定性与生物学结论分别见STATUS.json及中文报告；不以技术完成代表生物学通过。\n')
if mode=='final':
 review=R/'07_review';selected=[]
 roots=['00_manifest','01_reference','02_pseudospots','03_dev_mapping','04_spatial_mapping','05_evaluation','06_figures','code','logs']
 for root in roots:
  for p in (R/root).rglob('*'):
   if not p.is_file() or p.suffix in ['.h5ad','.npy','.npz','.pt','.ckpt'] or 'model' in p.parts or '__pycache__' in p.parts or p.name.startswith('checkpoint_mean_') or p.stat().st_size>20_000_000:continue
   selected.append(p)
 selected += [R/n for n in ['REPORT.html','README.md','METHODS.md','CONFIG.yaml','STATUS.json','MODEL_INDEX.tsv'] if (R/n).exists()]
 for p in selected:
  dest=review/p.relative_to(R);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
 manifest=pd.DataFrame([{'path':str(p.relative_to(review)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in review.rglob('*') if p.is_file() and p.name not in ['review.zip','PACKAGE_MANIFEST.tsv'] and not p.name.startswith('.')]);wt(manifest,review/'PACKAGE_MANIFEST.tsv')
 with zipfile.ZipFile(review/'review.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in review.rglob('*'):
   if p.is_file() and p.name!='review.zip' and not p.name.startswith('.'):z.write(p,p.relative_to(review))
 say('REVIEW ZIP', (review/'review.zip').stat().st_size)
say('REPORT UPDATED',mode)
