from base_data import *
def table(df):
 df=df.copy().fillna('—').astype(str);return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(row)+' |' for row in df.itertuples(index=False,name=None)])
def read(k):return pd.read_csv(O/'tables'/k,sep='\t')
report=['# HD 与 Xenium 注释及平台内共同 niche','**暂定计算结果，Codex 初审，病理/人工确认待返回。** 本轮止于身份及表达型共同 niche。','## 执行范围','Cervilla 和 D1_Johnson 完成原 2 µm 网格重建。UC_SCC、MPBC1、MPBC2 缺完整配准的 2 µm 输入。五张 GSE326226 Xenium 均沿用原分割。两个平台独立建模，各片输入及未纳入对象见 [输入审计](tables/input_manifest.tsv)。','## HD 重建、对照及限制']
rc=[];effect=[]
for s in HD:
 q=read(s+'_reconstruction.tsv');data=q.set_index('source');cons=json.loads((paths(s)/'conservation.json').read_text());img=pd.read_parquet(paths(s)/'secondary_image_support.parquet');n=data.loc['GEX_secondary'];rc.append((s,int(data.loc['HE_primary','objects']),int(n.objects),int(n.umi),int(n.median_UMI),int(n.median_genes),int(img.image_suspect.sum()),cons['unassigned_umi']))
 z=pd.read_parquet(paths(s)/'analysis_index.parquet');unknown=z.old_major.isin(['Unknown','Mixed']);fit=z.other_major.ne('no_fit');agree=z.ct_major_first.eq(z.other_major);resc=unknown&z.major_identity.eq(z.other_major)&z.major_identity.ne(z.ct_major_first)&z.review_level.isin(['fine_supported','major_only']);resolved=z[resc].major_identity.value_counts();effect.append((s,int(fit.sum()),round(agree[fit].mean(),3),int(unknown.sum()),int(resc.sum()),', '.join(f'{k} {v}' for k,v in resolved.head(3).items()),int((unknown&z.review_level.isin(['mixed_or_suspect','major_unresolved'])).sum())))
report.extend([table(pd.DataFrame(rc,columns=['切片','H&E 有成员对象','新增 GEX','新增 UMI','新增中位 UMI','新增中位基因','新增形态可疑','未分配 UMI'])),
'所有历史保留 primary 均按真实对象 ID 对齐并逐元素核对，见 [相同 primary 的 counts 和坐标](tables/primary_exact_comparison.tsv)。十个预定位置的 H&E-only 与旧 4 µm 扩张使用同一映射；GEX 只在 Cervilla 密集上皮核心窗口补入 2 个对象。全片补漏有 RNA 支持，且 primary 分配保持原样。六个按深度选取的 GEX 原图中，多个信号团块位于浅染区或组织/空腔交界，独立核归属存疑。对象面积来自观测格点。成分估计、标签和 niche 均不修复分割。没有逐核人工真值，无法估计漏检/错误合并真实发生率。[同位 ROI 对照](tables/hd_ROI_comparison.tsv) · [观测格点碎裂审计](tables/hd_geometry_QC.tsv)',
'## HD 参考路线与 TACCO',table(pd.DataFrame(effect,columns=['切片','TACCO 有拟合','两路线 dominant 一致比例','旧 Mixed/Unknown','其中采纳与分类器不同的 TACCO 大类','主要去向','仍混合或未定'])),
'TACCO 权重和第二成分保存在每片 tacco_weights.parquet。上表的采纳数需要独立原始大类 marker 支持，范围只限历史有标签对象。无 old_id 或旧对象未保留时另记 no_prior_label。充分的成分支持允许从不确定恢复大类；跨谱系冲突与无法确立物理边界的对象依然单列。各谱系/来源详情见 [TACCO 路线比较](tables/hd_TACCO_comparison.tsv)。',
'## Xenium 新旧身份与面板限制'])
z=read('xenium_major_transition.tsv');u=z[z.old_major.eq('Unknown')];uu=u.groupby(['new_major','new_level']).n.sum().reset_index().sort_values('n',ascending=False)
report.extend([f'旧 major 为 Unknown 共 {u.n.sum():,} 个；按共同 scVI 中的谱系内簇、CellTypist、原 Seurat 候选和原始 marker 支持更新。保留 major-only，防止稀疏细胞被迫成为细亚群。',table(uu.rename(columns={'new_major':'去向','new_level':'层级','n':'数量'})),
'[留出参考评估](tables/reference_evaluation.tsv) · [混淆矩阵](tables/reference_confusion.tsv) · [原高支持候选与新判读冲突](tables/xenium_old_confident_conflicts.tsv) · [簇级证据](tables/xenium_cluster_review.tsv)。原旧门控失败原因、数值差距和未测 anchor 见每片 old_gate_audit.parquet。无可用细类 anchor 的 SPP1 TAM、Venous EC、SMC 保留 panel-limited。参考表达可分性不代表组织内个体真值。','旧高支持大类也经过重新审阅。以旧支持 ≥0.8 且领先幅度 ≥0.25 为界，有 2,301 个对象的新大类或层级与旧判定不符；其中 455 个旧 Epithelial 转入 Mixed，324 个旧 Myeloid 转为 Unknown。它们是分歧复核候选。','## 对象覆盖'])
a=read('analysis_objects.tsv');report.append(table(a[['sample','original_units','niche_qualified','fine_supported','major_only','mixed','unknown']]))
report.append('## 两个平台的共同 niche')
for pl in ['hd','xenium']:
 tab=read(pl+'_niche_recurrence.tsv');ct=read(pl+'_niche_composition.tsv');names=[];st=read(pl+'_niche_stability.tsv');selected=len(tab.niche.unique())
 for k,sub in tab.groupby('niche'):
  f=ct[ct.niche_id.eq(k)].groupby('major_identity').n.sum().sort_values(ascending=False);total=f.sum()
  names.append((k,', '.join(f'{v} {n/total:.0%}' for v,n in f.head(3).items()),', '.join(sub.loc[sub.robust_presence,'sample']) or '—',', '.join(sub.loc[~sub.robust_presence,'sample']) or '—'))
 graph=read(pl+'_graph.tsv').query('hop==3');q=read(pl+'_niche_confounding.tsv').set_index('field')
 report.extend([f'### {pl}',f'同一平台模型选 K={selected}，种子间中位 ARI={st.loc[st.K==selected,"median_seed_ARI"].iloc[0]:.3f}。三跳邻域的逐片中位物理距离为 {graph.median_um.min():.1f}–{graph.median_um.max():.1f} µm。样本与 niche 的均衡抽样 NMI={q.at["sample","normalized_mutual_information_balanced_sample"]:.3f}，主要身份 NMI={q.at["major_identity","normalized_mutual_information_balanced_sample"]:.3f}。样本来源、未定身份或深度主导的区只作有限解释。',table(pd.DataFrame(names,columns=['共同编号','前三大类（全量）','同片 ≥50 且 ≥1%','检出但未达覆盖门槛'])),f'[逐片组成](tables/{pl}_niche_composition.tsv) · [技术指标与占比](tables/{pl}_niche_support.tsv) · [稳定性](tables/{pl}_niche_stability.tsv) · [表达嵌入与技术协变量](tables/{pl}_embedding_QC.tsv)'])
report.extend(['## 优先复核与下轮读取','1. 复核两片 GEX-only 的六张原图，重点 D1_Johnson 空腔/组织界面，确认有无同位置独立核。\n2. 查看 A2R2 冲突 ROI 的原多边形、A1R4 未定区及 A1R6 的低覆盖组织条带。T/NK 与 CAF/内皮的若干判定依赖面板中少数 marker。\n3. 优先比较高 Unknown 占比区及单片充足支持区，附带查看深度/面积。跨片检出不等同为跨患者验证。','[读入路径和字段](README.md) · [逐片分析对象](tables/analysis_objects.tsv) · [证据图](figures/) · [实际看图记录](shared/visual_review.md) · [方法和作者版本](METHODS.md) · [资源记录](logs/resources.tsv) · [交付检查](tables/delivery_checks.tsv)。当前挂载无 Git 工作树或经核实的本轮同步入口，交付于本地。'])
report.extend(['## 版本、恢复与交付状态',
'bin2cell 0.3.5、StarDist 0.9.2、TACCO 0.5.0、CellTypist 1.7.1、scvi-tools 1.5.1、CellCharter 0.3.9。旧 H&E 的复用、GEX-only 新重建、新谱系注释及共同 niche 的边界见 [方法](METHODS.md)。HD 最终 scVI 以每片 35,000 个对象均衡训练后全片编码；首次全量训练退出 143。TensorFlow 私有线程池和 TACCO 只读数组错误均已修复并重新运行。',
'全部计算使用 CPU，重计算使用同一组四个逻辑 CPU。共享资源组启用后记录到 19 tasks、7.50 GB 峰值，未发生 OOM。初段仅有逐阶段峰值，不能给出同口径全程总峰值。[模块状态](tables/module_status.tsv) · [逐阶段资源](tables/stage_resources.tsv) · [恢复日志说明](shared/runtime_notes.md)。',
'[全图索引](tables/figures.tsv) · [最少优先复核清单](tables/priority_review.tsv) · [HD 谱系内描述性 marker 排名](tables/hd_lineage_marker_ranks.tsv) · [Xenium 谱系内描述性 marker 排名](tables/xenium_lineage_marker_ranks.tsv)。'])
(O/'REPORT.md').write_text('\n\n'.join(report)+'\n')

