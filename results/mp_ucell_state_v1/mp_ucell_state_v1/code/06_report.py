from pathlib import Path
import csv, html, json, collections, datetime
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1')
def read(p):
 with (R/p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(rows,cols=None):
 if not rows:return '<p>无可评价记录。</p>'
 cols=cols or list(rows[0]); esc=lambda x:html.escape(str(x))
 return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(c)+'</th>' for c in cols)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+esc(fmt(r.get(c,'')))+'</td>' for c in cols)+'</tr>' for r in rows)+'</tbody></table></div>'
def fmt(v):
 try:
  x=float(v)
  return str(int(x)) if x.is_integer() else f'{x:.4f}'
 except:return v
parts=[]
def section(title,text):parts.append('<h2>'+html.escape(title)+'</h2>'+text)
def p(s):return '<p>'+s+'</p>'
def image(path,caption):return '<figure><img loading="lazy" src="'+html.escape(path)+'"><figcaption>'+html.escape(caption)+'</figcaption></figure>'
def link(path,title=None):return '<a href="'+html.escape(path)+'">'+html.escape(title or path)+'</a>'
figure_audit=read('00_input_audit/figure_file_validation.tsv')
assert figure_audit and all(x['pass_check']=='True' for x in figure_audit)
assert all(x['pass']=='TRUE' for x in read('00_input_audit/final_validation.tsv'))
assert all(x['pass']=='TRUE' for x in read('00_input_audit/PDF_parser_validation.tsv'))
assert (R/'logs/PALETTE_FINALIZATION_COMPLETE').exists()
assert not (R/'00_input_audit/missing_required_figures.txt').read_text().strip()
run=read('06_comparisons/run_summary.tsv');counts=read('00_input_audit/dataset_source_counts.tsv');cands=read('07_mp10_review/MP10_candidate_state_table.tsv');primary=sorted([r for r in cands if r['instance']=='FULL_Z_res0.4'],key=lambda r:float(r['raw_mean']),reverse=True)
comp=read('06_comparisons/ARI_NMI.tsv');rho=read('07_mp10_review/UCell_original_usage_spearman.tsv');qc=read('07_mp10_review/MP10_QC_spearman.tsv');rawqc=read('02_ucell_scores/raw_score_QC.tsv');vq=read('02_ucell_scores/vector_QC.tsv');validation=read('00_input_audit/final_validation.tsv')
section('技术执行状态',p('固定方案已执行：一次FULL UCell评分、8个直接14维建图/UMAP实例、10套Leiden标签。原RNA/NMF/注释/预后/空间输入只读。技术完成与生物学状态是否独立分离分开判断。')+p('共输出'+str(len(figure_audit))+'套PDF及300dpi透明PNG；文件配对、PNG透明通道/分辨率和全部必需dataset视图通过验证。')+table(run)+p('验证：'+str(sum(r['pass']=='TRUE' for r in validation))+'/'+str(len(validation))+'通过。详见'+link('00_input_audit/final_validation.tsv')+'。'))
section('1. FULL全集与样本/患者单位',p('FULL为原Seurat对象与冻结01_input_state.qs$cells完全匹配的123,699个细胞，57个生物样本、5个dataset。沿用原上皮谱系与Cycling输入，不依据恶性注释过滤；包含22,436个原Type=Normal细胞。冻结cell hash为36f7a452e661e558。原局部NMF仅保留50个样本，但本轮保留全部57个，其中7个少于20细胞。')+p('平衡单位为既有specimen审计中的orig.ident，使用与冻结analysis_sample逐细胞吻合的dataset|orig.ident全局key。每seed保留34,094细胞，每sample最多1000，未删小sample。BAL比例仅对应抽样集合，不代表FULL丰度。')+p('患者层面不可评价：现有映射中37个样本的患者ID为名称推导（derived），20个为临时（provisional），无独立确认的患者身份。保留原记录及多样本结构供审计，不把样本数当患者数，不声明患者等权或患者重复验证。')+table(counts)+p(link('00_input_audit/sample_patient_mapping.tsv')+'；'+link('00_input_audit/recorded_unverified_patient_structure.tsv')))
section('2. 14个冻结signature与评分质量',p('唯一来源为output_full_cc_seed5_v2/07_EpiMP_gene_basis.qs$sample_balanced_mean（3000×14；EpiMP01–14），与同目录CSV数值一致。按loading降序、基因名C/radix顺序处理并列；700个MP–gene条目与已有冻结top50逐项一致。每MP均为50个唯一基因、50/50精确存在于RNA；不做ID改名或补位，不使用NMF逐基因权重评分。MP间重叠保留。')+p('完整RNA统一counts背景为90,011个feature，稀疏输入；未缩为top50并集、HVG或3000个NMF基因。14个signature一次提交UCell，maxRank=1500，ties.method=average，missing_genes=impute，500-cell分块、2个worker，无平滑。M_Z使用FULL列均值和样本SD显式标准化，未截断。')+table(rawqc)+table(vq)+p('重复向量未被删除：若相同14维向量被图算法分到不同群，其界限不能解释为该14维分数支持的离散生物状态。')+p('各dataset nFeature<1500比例约1.9%–26.6%，低深度可能影响检出与活动分数；未因图形改变maxRank。详见'+link('02_ucell_scores/depth_by_dataset_sample.tsv')+'。'))
section('3. FULL_Z/res0.4群与MP10完整轮廓',p('全部群以C编号中性命名，不预设14群。以下按MP10原始UCell群均值排序；above_FULL_mean仅为描述性标记，不是验证亚群存在的阈值。MP10_rank_among_z比较同群14个MP的FULL全局z均值。')+table(primary,['cluster','n_cells','fraction_of_this_branch','raw_mean','raw_median','global_z_mean','MP10_rank_among_z','top3_z_MPs','n_samples','n_datasets','max_sample','max_sample_fraction'])+image('03_full_z/figures/FULL_Z__ALL__res0.4__cluster.png','FULL_Z主分支cluster UMAP')+image('03_full_z/figures/FULL_Z__ALL__res0.4__MP10.png','FULL_Z MP10原始UCell；精确0为黑色')+image('03_full_z/figures/FULL_Z__ALL__res0.4__cluster14MP_GLOBALz_heatmap.png','全部群的14MP全局z均值；z=0表示FULL均值，不表示无表达')+p(link('03_full_z/cluster_MP_profiles_res0.4.tsv','每群完整14MP raw均值/中位数/分位数与全局z均值')))
section('4. 来源重复性与混杂',p('样本和dataset分布、最大样本贡献以及每个样本内群vs其余细胞的MP10差异均已输出。少量细胞的“出现”不等于样本或患者重复验证。以下技术相关为Spearman，不推断因果。患者层面因映射未确认而不可评价。')+table([r for r in qc if r['unit']=='FULL'])+table(read('07_mp10_review/MP10_by_dataset_source.tsv'))+p(link('03_full_z/MP10_sample_within_contrasts_res0.4.tsv')+'；'+link('03_full_z/cluster_by_dataset_res0.4.tsv')+'；'+link('03_full_z/cluster_by_sample_res0.4.tsv')))
section('5. 尺度、抽样与resolution敏感性',p('ARI/NMI均仅按真实共同cell ID计算；NMI采用两分区熵算术平均归一化。不同seed补充比较仅使用交集，n_common明确报告。FULL_Z 0.2/0.6复用同一SNN及UMAP。群编号不跨分支等同。最佳Jaccard在共同细胞内计算，拆分与合并须同时看归一化重叠和14MP组合。未以任何ARI阈值判定亚群存在，也未挑选最佳seed。')+table(comp)+p(link('06_comparisons/best_cluster_Jaccard.tsv')+'；'+link('06_comparisons/matched_branch_14MP_profiles.tsv')))
section('6. 原NMF usage衔接',p('仅使用08_cell_EpiMP_scores.qs$usage_raw，原为14×cells，经维名核对转置后按cell ID匹配。常数或有效N不足时相关记NA。总体14MP及各sample相关均已保存。该一致性是新旧分析衔接，不是外部验证。')+table([r for r in rho if r['unit']=='FULL'])+image('03_full_z/figures/FULL_Z__original_MP10_usage_raw.png','原MP10 usage_raw叠加到本轮FULL_Z embedding'))
interpret=R/'07_mp10_review/interpretation_final.html'
section('7–8. 状态判断与下一轮marker优先顺序',interpret.read_text() if interpret.exists() else p('结果解释尚未最终审阅；请勿视为已完成的生物学结论。'))
section('方法、执行边界与未完成项',p('框架借鉴用户指定的Hu等2026描述“MP top50→UCell→基于分数的降维和无监督分群”；本项目预设作者未公开的参数，不宣称完全复现作者代码。直接14维欧氏距离：Seurat FindNeighbors(k=20,annoy,n.trees=50,l2.norm=FALSE,prune.SNN=1/15)；Leiden algorithm=4，leidenbase/RBConfigurationVertexPartition，seed=42，主resolution=0.4；uwot UMAP(n.neighbors=30,min.dist=0.3,2D,euclidean,seed=42)，单线程，pca=NULL。')+p('未运行单细胞CNV、bulk预后、CellChat、空间映射、正式marker命名或新亚群预后。探索性全转录组marker非本轮必需，留至人工讨论候选群后。患者统计/患者UMAP因现有患者映射未验证而未生成。未定位到bulk weighted_top50原表，已记录搜索范围；不影响canonical basis/top50一致性审计。')+table(read('00_input_audit/software_versions.tsv'))+p('参考：'+link('https://www.bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html','UCell参数')+'；'+link('https://satijalab.org/seurat/reference/findneighbors','FindNeighbors')+'；'+link('https://satijalab.org/seurat/reference/findclusters','FindClusters')+'；'+link('https://satijalab.org/seurat/reference/runumap','RunUMAP')))
keypaths=['README.txt','config.R','00_input_audit/input_manifest.tsv','00_input_audit/cell_manifest.tsv.gz','01_signatures/top50_signatures.tsv','02_ucell_scores/M_RAW.qs','02_ucell_scores/M_Z.qs','02_ucell_scores/global_scaling.tsv','05_balanced/seed42_cell_manifest.tsv.gz','05_balanced/seed43_cell_manifest.tsv.gz','05_balanced/seed44_cell_manifest.tsv.gz','06_comparisons/labels_all10.qs','06_comparisons/ARI_NMI.tsv','07_mp10_review/MP10_candidate_state_table.tsv','08_figures_by_dataset','code','logs']
section('关键输出绝对路径','<ul>'+''.join('<li>'+link(str(R/x))+'</li>' for x in keypaths)+'</ul>')
body='<!DOCTYPE html><html lang="zh-CN"><meta charset="utf-8"><title>14个冻结MP的UCell状态分群</title><style>body{background:#fff;font-family:system-ui,sans-serif;max-width:1300px;margin:40px auto;padding:0 24px;line-height:1.65;color:#152536}h1,h2{color:#183f61}h2{border-top:1px solid #ddd;padding-top:20px}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ddd;text-align:left}th{background:#edf4f8;position:sticky;top:0}.scroll{overflow:auto;max-height:760px}img{max-width:100%;height:auto}figure{max-width:1100px;margin:20px 0}figcaption{color:#5b6770}a{color:#006da8}code{overflow-wrap:anywhere}</style><h1>14个冻结MP的UCell状态分群</h1><p>运行目录：'+html.escape(str(R))+'</p>'+''.join(parts)+'</html>'
(R/'REPORT.html').write_text(body)
print(R/'REPORT.html')
