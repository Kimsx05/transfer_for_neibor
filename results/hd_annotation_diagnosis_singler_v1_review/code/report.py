from pathlib import Path
import pandas as pd,json,datetime,zipfile
root=Path(__file__).resolve().parents[1]
def md(d):
 d=d.copy().fillna('NA');return '| '+' | '.join(map(str,d.columns))+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(str(v).replace('|',' / ') for v in row)+' |' for row in d.itertuples(index=False,name=None))+'\n'
summary=pd.read_csv(root/'summary.tsv',sep='\t');lines=['# Unknown 来源诊断与 SingleR 对照','',f'更新：{datetime.datetime.now().isoformat(timespec="seconds")} (服务器本地时间)。两片全部原 qc_keep_analysis 对象已完成；本轮候选不替换原正式注释。','',f'结果根目录：`{root}`。主输入位于 `../hd_major_annotation_v1/<section>/seurat.rds`；预处理保持 `../segmentation_assignment_v2/`。','',md(summary[['section','n','UMI','original_unknown','unknown_UMI','unknown_BLCA_nonpruned','unknown_HPCA_nonpruned','unknown_both_agree','unknown_both_disagree','unknown_neither_nonpruned']]),'一致数仅作描述，不设“两方法一致”的新接纳门槛。未被 pruning 的候选也不等于高置信度或已救回身份。','', '## 原规则追踪与主要结论','',
'原 Seurat 使用冻结 BLCA 44,896 个参考细胞，训练前已归并为九大类，直接 TransferData 大类预测；并非细亚型先预测再归并。原脚本没有调用 GetTransferPredictions，因此不存在实际 score.filter，也不是被其默认 0.75 截断。原始预测均有类别，Unknown 产生于随后 marker/一致性规则。',
'主接纳要求逐对象：(1) marker 强度最大谱系与参考原预测一致；(2) 预测谱系至少检出两枚 panel marker；(3) 不触发跨谱系竞争。竞争定义为至少两组各检出≥2枚且强度≥最强组的70%，仅 Fibroblast+Mural 共同收缩表达时豁免。Plasma 另要求 MZB1/JCHAIN/TNFRSF17/DERL3 至少一个 core marker，否则其 panel 证据置零。该竞争指标是诊断规则，不是已验证的混合细胞或双细胞真值。',
'参考支持带取 marker 一致且 nFeature≥30 对象的类内分数/分差第10百分位，稀少类用合并分布；主要决定 supported 与 low_confidence。主接纳对象即使低于带仍有暂定标签。另有“一枚 marker + 强参考支持 + 簇参考一致率≥0.8”的低置信度补充，以及至少两枚对象 marker 和两枚 cluster-enriched marker 的低置信度参考不匹配修订。原未通过 marker 一致性和数量规则的对象，不能简单归因于参考置信阈值。',
'原 marker 计算使用完整分析集合的 RNA 原始 counts 检出与 RNA LogNormalize data 强度；不是 sketch、HVG 或 scale.data。两片全部标签、置信状态与冲突标志已根据原代码和缓存精确重现。原 annotation_reason 为历史记录；diag_first_rejection_reconstructed 按书写的合取顺序“谱系一致→数量→冲突”重建，仅对最终 Unknown 赋值。矢量代码没有历史短路执行日志，故不将此字段冒充原始运行日志。原 reason 的赋值覆盖顺序与首次失败不同。',
'共有基因总体匹配正常，但BLCA共有panel缺KRT18及Mast的TPSAB1/TPSB2；HPCA还缺KRT8、TRBC2及部分免疫球蛋白基因。覆盖不足是具体局限，不能称所有marker完整。原RNA的完整panel统计用于独立复核这些候选。',
'原报告区分了支持性 marker 不足和竞争，但两者都会产生 Unknown。关闭单规则的结果只评估敏感性，保持原支持带不变，不自动改正式标签。重叠规则不能直接相加。','', '## 两路 SingleR 设置与适用范围','',
'固定 seed=123，SingleR 2.12.0；4 个计算线程，BLAS 1，逐片处理，query 每块最多5,000对象。每个切片×参考的 classifySingleR(prune=FALSE) 结果按 cell_id 合并后仅统一调用一次默认 pruneScores(nmads=3,min.diff.med=-Inf,min.diff.next=0)。验证时由保存的完整 scores、labels、delta.next 再现全片 pruning。scores 是 fine-tuning 前相关性分数，不是概率；其最大分列可能不同于最终 labels。保留官方 delta.next 与相对中位数 delta，未以普通 top1−top2 代替。',
'BLCA：完全复用原参考缓存的44,896细胞及标签；9类均覆盖，14,325个可靠共有 Ensembl 基因。trainSingleR(de.method="wilcox",de.n=50,aggr.ref=TRUE)，查询仍逐对象，非簇级。其与原 Seurat 的一致性是同参考下的方法稳健性证据，不是独立准确率验证。原参考来源与细标签到大类映射见 reference/blca_label_mapping.tsv；原 reference.rds 保留在上一轮。',
'HPCA：celldex::HumanPrimaryCellAtlasData()，隔离安装的 celldex 1.12.0，ExperimentHub HPCA 1.0.0 (EH3492/3508 logcounts)，复用官方下载缓存。网络中断通过断点续传恢复，无替换或合成参考。作者指定10个 label.main 子集共448样本，训练前映射到6大类；19,363参考基因中15,916与HD symbol共有。trainSingleR(de.method="classic",aggr.ref=FALSE)。不覆盖 Mast、Plasma、Mural；这些 BLCA 候选的 HPCA 分歧须单列，不能解读为类型不存在或成功救回。作者代码对部分T/B使用细标签，本轮按大类；本轮是重建细胞的RNA对照，并非8um方格+SCT的严格复现。',
'查询来自完整18,085基因 RNA LogNormalize data；没有额外QC、HVG限制或重新整合。BLCA按Ensembl ID匹配。HPCA临时视图将TBCE、HSPA14、TMSB15B各自的重复symbol行先对原始counts求和 (总计合并3行)，再以1e4总量LogNormalize；不改写原对象。模型内部为分类选取的marker基因与“只用HVG输入”不同。各片 gene_marker_coverage.tsv 记录完整共有基因和模型panel覆盖。模型基因表达全零/恒定者保留原ID，标为不可评估；不会强行赋值或悄悄掉行。','']
for _,r in summary.iterrows():
 s=r.section;o=root/s;d=pd.read_csv(o/'diagnosis.tsv.gz',sep='\t',low_memory=False);reasons=pd.read_csv(o/'original_unknown_reasons.tsv',sep='\t');qc=pd.read_csv(o/'QC_distributions.tsv',sep='\t');abl=pd.read_csv(o/'single_rule_ablation.tsv',sep='\t');first=pd.read_csv(o/'first_rejection_reconstructed.tsv',sep='\t');cand=pd.read_csv(o/'unknown_candidate_summary.tsv',sep='\t');roi=pd.read_csv(o/'roi_comparison.tsv',sep='\t');cl=pd.read_csv(o/'cluster_diagnosis.tsv',sep='\t');u=d.final_major.eq('Unknown')
 lines += [f'## {s}','',f'固定分母 {r.n:,} 对象 / {r.UMI:,} 原始UMI。原Unknown {r.original_unknown:,} ({r.original_unknown/r.n:.2%})，承载 {r.unknown_UMI:,} UMI ({r.unknown_UMI/r.UMI:.2%} 全片输入)；原已注释低置信度 {(d.annotation_status=="low_confidence").sum():,}。','', '原历史原因（互斥最终原因）：',md(reasons.round(4)),'重建首次拒绝（不能与上表相加）：',md(first),'单规则关闭的敏感性（只诊断）：',md(abl[['removed_rule','original_unknown_get_label','original_unknown_UMI_get_label','status_changed']]),f'全panel强度为0的原Unknown有 {int((u & d.diag_all_marker_strength_zero).sum()):,}，这类的argmax并列落在首个谱系；其重建的marker agreement失败并不代表存在相反谱系证据。原Unknown中任一组有≥2枚marker者 {int((u & d.diag_any_group_two_markers).sum()):,}，仍应区分缺证据与竞争。','信息量中位数：',md(qc[['status','metric','50%','25%','75%']]),'Unknown最多的五个簇（按对象数；原簇不改动）：',md(cl.sort_values('unknown',ascending=False).head(5)[['cluster_projected','n','unknown','unknown_fraction','median_nFeature','major_raw','main_original_rejection']].round(3))]
 for route in ['blca','hpca']:
  m=json.loads((o/f'{route}_prediction_manifest.json').read_text()); lines += [f'{route.upper()}：可评估 {m["evaluated"]:,}；不可评估 {m["not_assessable"]:,}；全片统一 pruning 排除 {m["pruned"]:,}；初始分数最高类别与最终 fine-tuning 类别不同 {m["score_top_not_final"]:,}。',md(cand[cand.route==route].drop(columns='route').round(4))]
 lines += [f'原Unknown中两路一致 {r.unknown_both_agree:,}，两路都有候选但不同 {r.unknown_both_disagree:,}；两路均无未pruned候选 {r.unknown_neither_nonpruned:,}。BLCA候选属于HPCA未覆盖三类的有 {r.unknown_BLCA_not_HPCA_covered:,}，参见 coverage_stratified_comparison.tsv。后者与“双方都有候选但不同”可重叠，不能相加。',f'原Unknown的未pruned BLCA候选承载 {r.unknown_blca_candidate_UMI:,} UMI；HPCA候选承载 {r.unknown_hpca_candidate_UMI:,} UMI。所有 method_composition/label_crosstabs 均同时报告对象数、UMI与固定全片分母比例。','', 'ROI同坐标复核统计：',md(roi),'关键图：',f'![原标签与BLCA交叉]({s}/figures/cross_Original_final_BLCA_pruned.png)',f'![Unknown候选marker]({s}/figures/markers_unknown_blca.png)',f'![BLCA全片候选]({s}/figures/spatial_candidate_blca.png)',f'![间质ROI三路对照]({s}/figures/ROI_stroma.png)','']
review=root/'spatial_review.md'
if review.exists():lines += [review.read_text(),'']
lines += ['## 解读与下一轮建议','',
'两片 Unknown 是信息不足、参考与marker不一致及竞争规则共同作用，不能单独归因于默认 Seurat score.filter。严格逐对象marker要求确实压低标签覆盖，但关闭冲突规则的边际影响最大；这不证明被拒绝的竞争对象是真混合，也不证明全部应恢复身份。低检测量造成“未检出”与足够检测量下其他谱系成组表达，应继续分开标记。',
'默认SingleR pruning在这批低深度、分布宽的数据上没有提供有效的绝对置信门槛（实际排除数见各片manifest）。阈值为每类delta中位数减3MAD，可为负；这是默认统计机制的实际结果，不为得到目标Unknown比例而改阈值。大量未pruned候选只应作为下一轮修订证据，不能宣称Unknown已被解决或准确率提升。',
'建议下一轮先评估成组marker均值/检测率与原表达簇支持的候选，逐对象保留信息量、替代谱系证据、参考支持与候选来源。Mural维持独立；HPCA不裁决Mural/Mast/Plasma。B/Plasma及T/NK等稀少候选尤其检查成组表达和低信息比例。不要新增“必须两法一致且每个对象检出全部marker”的门槛；不要用整簇或空间邻居多数票填满。',
'只有在表达量充足、不同谱系成组表达持续竞争，且ROI上位于8um归属交界时，才值得下一轮局部核对原bin→cell映射；本轮证据不足以把所有Unknown归咎于核扩张。当前不需要重跑全片分割、归属或正式下注释。','', '## 文件与读取','',
'每片：`diagnosis.tsv.gz` 为唯一逐对象压缩诊断表，以 section+cell_id 为键；原字段保留，新诊断使用 diag_、singler_blca_、singler_hpca_ 前缀。`singler_*_scores.rds` 为轻量分类分数，行名即cell_id。`validation.json` 记录按ID附加字段、原counts/QC/簇/UMAP一致及全片pruning复现。参考细胞数、映射和基因覆盖在 reference/ 及各片 gene_marker_coverage.tsv。',
'汇总包括重叠规则、首次拒绝、单规则关闭、nFeature分层、QC分布、cluster诊断、15组两两标签交叉（含对象与UMI）、候选marker均值/阳性率、覆盖分层、ROI统计。candidate_marker_statistics.tsv 使用全分析集合原RNA data与counts>0，按原Unknown/已知分组；不是插补表达，也未重跑全量差异分析。',
'绘图复用 jins-routine-plotting；白底Arial约12pt，PNG300dpi；有图例的PDF另含legend_only、panel_only。统计图尽量矢量；H&E和高密度点层为栅格。UMAP完全沿用旧坐标，各片独立；Cervilla浅蓝遮罩表示无supplied bins，不能解释为Unknown/零表达组织。ROI青色线是原核边界，灰线是8um归属边界，圆点填色才代表标签，不能将青色核边界误认Mast。ROI三面板统一坐标，Original/SingleR BLCA/SingleR HPCA仅标明身份，无总结性标题。',
'```r\nroot <- "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis"\nsource(file.path(root,"hd_major_annotation_v1/code/load_result.R"))\nx <- load_major(file.path(root,"hd_major_annotation_v1/Cervilla"))\nsource(file.path(root,"hd_annotation_diagnosis_singler_v1/code/load_result.R"))\nx <- attach_diagnosis(x, file.path(root,"hd_annotation_diagnosis_singler_v1/Cervilla"))\n# x$final_major 保持原值；候选在 x$singler_blca_pruned_labels 等字段。\n```',
'`bash rerun.sh` 复用冻结模型、完成的预测和图表，更新报告与review.zip；若输入路径/大小/mtime或核心脚本/配置改变，缓存检查会停止，需有意重建相应步骤，不能静默混用旧统计。只保存小型预测、诊断与参考训练缓存，没有新的完整Seurat/h5ad/全基因矩阵或图像副本。','',
'`review.zip` 含报告、代码、汇总、压缩逐对象诊断、分类分数与核心图 (组织/UMAP的PNG及统计图PDF；其他PDF仍在服务器)；不含原始/标准化矩阵、主Seurat、参考/模型缓存。目录外原输入保持不变。','', '## 官方/作者接口依据','',
'- [作者 SingleR 脚本](https://github.com/scervilla/SpatialBenchmarking/blob/main/analysis/VisumHD/singleR.R)；本地只读副本 reference/author_singleR.R。',
'- [SingleR 单细胞参考与聚合](https://bioconductor.org/books/release/SingleRBook/sc-mode.html)。',
'- [SingleR 分数与 pruning](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html)。',
'- [GetTransferPredictions 接口](https://satijalab.org/seurat/reference/gettransferpredictions)；本轮历史审计确认原流程未使用该函数。']
(root/'REPORT.md').write_text('\n\n'.join(lines))
# One lightweight review archive; references and models explicitly excluded.
with zipfile.ZipFile(root/'review.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as zz:
 for p in root.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(root);parts=rel.parts
  if parts[0] in ['reference','Rlib','logs'] or p.name=='review.zip':continue
  if p.suffix=='.rds' and not p.name.endswith('_scores.rds'):continue
  if 'figures' in parts and (p.name.endswith('_legend_only.pdf') or p.name.endswith('_panel_only.pdf')):continue
  if 'figures' in parts and p.suffix=='.pdf' and p.name.startswith(('ROI_','spatial_','umap_')):continue
  zz.write(p,rel)
print('REPORT_AND_REVIEW_UPDATED')
