from pathlib import Path
import csv,gzip,html,collections,statistics,json,datetime
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_followup_v2')
def rd(p):
 p=R/p
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def wr(p,rows):
 if not rows:return
 with (R/p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def table(rows,cols=None):
 if not rows:return '<p>无可评估结果。</p>'
 cols=cols or list(rows[0]);return '<table><thead><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in cols)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(c,'')))+'</td>' for c in cols)+'</tr>' for r in rows)+'</tbody></table>'
def link(p,label=None):return '<a href="'+p+'">'+html.escape(label or p)+'</a>'
def num(v,d=3):
 try:return f'{float(v):.{d}g}'
 except:return str(v)
models=rd('02_pseudobulk_markers/model_summary.tsv');elig=rd('00_audit/comparison_eligibility.tsv');usage=rd('01_usage_bridge/MP10_sample_equal_summary.tsv');comp=rd('03_marker_review/source_composition.tsv');markers=rd('03_marker_review/display_markers.tsv');qc=rd('03_marker_review/QC_paired_summary.tsv');short=rd('03_marker_review/all_shortlists.tsv')
# Full 57 sample eligibility rows per comparison, including zero-group library sizes.
libraries=[]
for c in [r['comparison'] for r in models]:
 pb=rd('02_pseudobulk_markers/'+c+'_pseudobulk_metadata.tsv'); ix={(r['sample'],r['group']):r for r in pb}
 for r in elig:
  if r['comparison']==c:
   row=dict(r);row['library_left']=ix.get((r['sample'],'Left'),{}).get('library_size','0');row['library_right']=ix.get((r['sample'],'Right'),{}).get('library_size','0');libraries.append(row)
wr('00_audit/comparison_eligibility_with_libraries.tsv',libraries)
# Every counterexample among qualified pairs, with expected direction explicit.
pairs=rd('03_marker_review/display_marker_all_sample_directions.tsv');failed=[]
for m in markers:
 state=m['state'];allowed={'C8':['r04_C8_vs_rest','r04_C8_vs_C6','r04_C8_vs_C3'],'C12':['r06_C12_vs_rest','r06_C12_vs_C16'],'C16':['r06_C16_vs_rest','r06_C12_vs_C16']}[state]
 for p in pairs:
  if p['gene_id']!=m['gene_id'] or p['comparison'] not in allowed:continue
  sign=-1 if state=='C16' and p['comparison']=='r06_C12_vs_C16' else 1
  if float(p['delta_expression'])*sign<=0:
   failed.append(dict(state=state,gene_id=m['gene_id'],comparison=p['comparison'],sample=p['sample'],dataset=p['dataset'],source=p['source'],n_left=p['n_left'],n_right=p['n_right'],expected='Left lower' if sign<0 else 'Left higher',delta_expression=p['delta_expression']))
wr('03_marker_review/display_marker_counterexamples.tsv',failed)
methods='''# 方法与解释范围

## 冻结输入
所有原NMF、UCell、聚类和降维均只读。FULL严格沿用冻结01_input_state.qs$cells和旧cell_metadata；123699细胞、57个specimen、5个dataset，保留22436个Normal。按cell_id连接标签/分数/坐标。C8来自r0.4，C12/C16为r0.6完整群，两套resolution分别处理。sample=dataset|orig.ident。37个sample的患者映射为derived、20个为provisional；41个recorded patient key不是41位已核实独立患者，其中16个key各对应两个sample。

## 原usage衔接
只使用usage_raw，未变成行和比例。每群14MP的均值、中位数、Q25/Q75、零比例同时保存usage及UCell。六项比较均正向=左侧高；各组至少30细胞的sample可配对，使用组内全部细胞，不重采样。完整57样本资格及文库大小另存。pooled细胞加权结果与sample等权结果分开。样本均值双侧配对Wilcoxon，exact=FALSE，至少3对；全零/不足样本返回NA及原因。BH在每个comparison×scale内14MP单独校正；dataset/source分层表同样在各自层内14MP校正，不声称整个多表筛选的家族错误受控。冻结MP10-high复用10_high_usage_global_cutoffs.csv的原全局90% quantile阈值和usage>=阈值定义；此轮没有重新选择阈值。

## RNA与pseudobulk
原Seurat RNA单一counts layer，90011 features×123699 cells，388983646非零项，全部非负整数。未把完整表达矩阵转dense。已有RNA data在固定25个细胞的所有features上与log1p(counts/library×10000)一致，最大差8.88e-16，因此复用普通LogNormalize data做展示。先用稀疏指示矩阵聚合sample×原cluster的count sum、LogNormalize sum、detected cell sum，再按固定比较合成sample×group。未复制完整Seurat；原对象仅一次加载。

edgeR 4.8.2/limma 3.66.0，先保存并检查本地官方help。filterByExpr(y,group=group)默认min.count=10,min.total.count=15,large.n=10,min.prop=0.7；过滤按两大group，不用含sample阻断的高杠杆design放宽检测要求。过滤后重新计算文库大小、TMM归一化；estimateDisp(robust=TRUE)、glmQLFit(robust=TRUE,legacy=FALSE)、glmQLFTest。设计~sample+group，Right reference，检验groupLeft系数；满秩、文库非零、残余自由度>0均核对。dataset/source在每个sample内唯一，已被sample固定效应吸收，不再加共线主效应。各比较在全部通过过滤的基因内BH，完整结果/过滤记录保留。没有基于细胞P值的DE。即使全部模型运行成功，也不意味着统计独立性、混杂或分群后选择偏倚已经消除。

## marker审阅
完整基因ID不变；gene_symbol仅对已知SYMBOL或唯一ENSEMBL映射填写，歧义保留NA，不合并行。表达与检测率先按sample×group求均值，再对该比较合格样本等权。检测率即RNA counts>0的细胞比例，不等同蛋白阳性。短名单门槛为FDR<0.05、|log2FC|≥0.5、目标平均检测率≥0.10、样本内表达方向支持≥0.70；排序依次为dataset均值差方向支持、sample方向支持、朝目标方向的检测率差、|log2FC|、gene_id。每comparison至多20，C12 vs C16每方向至多10，其余只取左侧；这是一套透明的实用筛选规则，不是生物学验证阈值或综合验证分。

展示marker每state至多3个，在固定短名单与冻结基因诊断结果中人工复核局部对照及功能多样性，不更改DE或群。C8的ELOVL5来自冻结诊断/符合条件但未进入该比较前20的候选，其加入为展示脂质支持，而非改变短名单排序。输出其余原群中最接近/更高表达群，不新增全群两两DE。同分辨率所有非空sample×cluster（含小群）用于描述性特异性汇总，n_samples/n_cells均保留；因此此汇总不能替代严格配对局部对照。MP07-only/MP10-only/shared=38/38/12；逐基因表达与检测率，不重算任何UCell/删共享后分数。功能仅用org.Hs.eg.db 3.22.0及GO.db本地注释的有限关键词记录（含GO ID/term/evidence），没有通路富集网格。

## 图形与复现
冻结FULL_Z UMAP及原Seurat umap_allgene按cell_id取子集，不重拟合。点序seed42；精确零黑色并先画，正值用plasma；NA灰色。颜色范围跨dataset固定并单独记录；usage连续图为明确标注的log1p显示，统计仍是raw。RNA DotPlot为sample等权普通log-normalized RNA，不是usage；大小为平均检测率。QC counts/features展示log1p(mean)，percent.mt保留百分数；无缺失填0。PDF和300dpi透明PNG，每图保留实际数据和说明。GSE326225保留所有细胞/图及零资格记录，但六项比较无合格样本，不伪造配对或DotPlot。最多2线程，原Seurat读取/模型顺序执行。没有重新运行NMF/UCell/UMAP/聚类、空间/CNV/预后或其他任务。

## 参考
- [edgeR官方入口](https://bioconductor.org/packages/edgeR/)；执行API以00_audit中edgeR4.8.2官方help为准。
- [OSCA配对多样本比较](https://bioconductor.org/books/3.17/OSCA.multisample/multi-sample-comparisons.html)
- [生物学重复与单细胞DE](https://www.nature.com/articles/s41467-021-25960-2)
- [EGR1原始功能研究](https://pubmed.ncbi.nlm.nih.gov/8336701/)
'''
(R/'METHODS.md').write_text(methods)
decisions='''# 本轮决定与解释边界

优先保留完整C8作为下一轮独立marker评估的主候选范围，C12与C16同时保留为需要区分的表达状态。当前更支持一个有内部表达差异、边界可能连续的MP10相关程序集合；不能据此宣布两个独立且已验证的亚群。没有以UMAP岛形态作否定或通过依据，也没有根据结果选择stable core或新seed。

- C8：原MP10 usage在20/28样本更高；EGR1、GADD45B、ELOVL5优先展示。EGR1局部C6/C3有支持，但属于诱导应答信号；GADD45B在C13更高；ELOVL5为MP07/MP10共享输入基因，不能当独立验证。
- C12：HIPK3、GSE1用于观察与C16的区分，ELOVL5用于脂质支持。前两者在10/10配对样本高于C16，但其他群可更高；ELOVL5本身不能区分C12/C16。不将C12自动升级为独立脂质亚群。
- C16：NRARP、IDI1、LDLR保留展示。CAPG/RAB25/TMEM256虽区分C12，却不显著高于rest且其他群更高，不优先作为特异marker。C16近七成来自HRA两个sample，且配对深度更高；暂不优先作为单独的新亚群进入机制/预后研究。

所有输入MP基因参与了建图；非输入HIPK3/GSE1/NRARP仍来自同批发现RNA。患者映射未核实，16个recorded key各有两个sample，模型尚未充分处理潜在同患者相关性。配对样本少和来源覆盖不齐是证据限制，不能用细胞总量替代重复数。

技术问题已修正：QC汇总整数/浮点类型；data.table配置序列化的内部属性比较；短名单检测率差需按目标方向取符号；绘图色阶需覆盖各dataset均值。拼图容器的透明背景也经PNG像素检查后修正。上述修复没有修改既有输入/标签/基因集或重跑评分/图聚类。原NMF与原预后效应的定义保留，不能转写为C8/C12/C16的效应。
'''
(R/'DECISIONS.md').write_text(decisions)
summary=[]
for m in models:
 u=next(r for r in usage if r['comparison']==m['comparison'] and r['scale']=='usage_raw' and r['stratum']=='ALL')
 summary.append({'比较':m['comparison'],'配对sample':m['n_pairs'],'dataset':m['n_datasets'],'检验基因':m['n_tested'],'usage正/负':u['n_positive']+'/'+u['n_negative'],'usage差中位数[IQR]':f"{num(u['delta_median'])} [{num(u['delta_q25'])}, {num(u['delta_q75'])}]",'探索性usage FDR':num(u['FDR'])})
wr('01_usage_bridge/review_summary.tsv',summary)
bpair=[{'sample':r['sample'],'来源':r['source'],'C12':r['n_left'],'C16':r['n_right']} for r in elig if r['comparison']=='r06_C12_vs_C16' and r['eligible']=='TRUE']
composition=[{'state':x['target'],'dataset':x['level'],'细胞数':x['target_n'],'占该群%':num(float(x['fraction_of_target'])*100),'占dataset完整池%':num(float(x['fraction_of_source_pool'])*100)} for x in comp if x['field']=='dataset']
marktab=[{'state':x['state'],'gene':x['gene_id'],'角色/限制':x['role'],'更高其他群':x['highest_other_cluster'],'其他群≥目标数':x['n_equal_or_higher'],'参与14MP':x['in_any_14MP_top50']} for x in markers]
qtab=[{'比较':x['comparison'],'QC':x['metric'],'配对数':x['n_pairs'],'左减右中位差':num(x['median_delta']),'正/负':x['positive']+'/'+x['negative']} for x in qc if x['comparison']=='r06_C12_vs_C16']
figs=['ALL__FULL_Z__frozen_state_overlays','ALL__original_RNA__frozen_state_overlays','ALL__MP10_usage_paired_deltas','ALL__paired_QC','ALL__shortlist_DotPlot','ALL__MP07_MP10_gene_diagnostics']
body=f'''<h1>MP10候选状态后续审阅</h1><p>完成日期：2026-10-08。分析目录：<code>{R}</code></p>
<div class="key"><p><b>主要判断：</b>C8具有跨样本的原MP10 usage及部分非共享脂质/甾醇基因支持；C12/C16存在可重复表达差异，但深度和来源限制明显。优先保留完整C8作下一轮独立marker评估，C12/C16作为并行比较状态。目前不足以确认两个独立生物学亚群，更适合描述内部异质、可能连续的MP10相关程序。</p><p><b>技术状态：</b>六项固定配对pseudobulk全部完成。既有NMF、14MP、分数、标签和两套UMAP均未重算或修改。输出51套PDF/透明300dpi PNG，全部附绘图数据。证据不足不是程序失败。</p></div>
<h2>1. 输入与细胞身份</h2><p>FULL=123,699细胞、57个biological specimen、5个dataset，保留22,436个Normal。C8=8,603，C12=5,596，C16=2,779；交集分别5,444与2,674。C12/C16使用完整原群；与C8的重叠不是独立复现。两套resolution分别统计，没有三组互斥检验。</p><p>原RNA counts为90,011×123,699稀疏整数矩阵；普通RNA LogNormalize data验证通过。41个recorded patient key均未独立核实（37个sample derived，20个provisional）；16个key各有2个sample。只能称跨样本，不能称独立患者验证。</p><p>{link('00_audit/identity_checks.tsv')} · {link('00_audit/input_manifest.tsv')} · {link('00_audit/sample_mapping.tsv')} · {link('00_audit/counts_check.tsv')}</p>
<h2>2. C8原usage是否也高？</h2><p>C8 vs rest在28个合格sample中20正、8负；差中位数2.876，IQR −0.201至11.518，sample等权平均usage为12.546 vs 5.542。四个dataset的正向数分别CNP 6/10、HRA 6/8、Inhouse 3/4、PRJNA 5/6。来源分层Normal 5/9、MIBC 9/12、NMIBC 6/7。UCell为28/28正，因此原usage衔接支持是多数样本而非普遍一致。</p><p>C12 vs C16有10对，其中原usage和UCell均7对为C16更高；原usage差中位数−4.078，探索性FDR0.177。C16 vs rest有10/12正但FDR0.152；不要把方向支持与统计不确定性混为一谈。完整14MP概览可判断共同程序，未把其他usage归一化成MP10比例。</p>{table(summary)}<p>{link('01_usage_bridge/MP10_sample_equal_summary.tsv')} · {link('01_usage_bridge/pooled_cell_weighted_14MP.tsv')} · {link('01_usage_bridge/frozen_MP10_high_coverage.tsv')} · {link('00_audit/comparison_eligibility_with_libraries.tsv')}</p>
<h2>3. C12/C16表达差异及真正提供证据的来源</h2><p>C12 vs C16的10对来自CNP0000460 5、HRA003620 4、PRJNA662018 1；Normal3、MIBC5、NMIBC2。Inhouse和GSE没有合格C12/C16配对，不能从其少量细胞声称重复。该比较测试7,298基因，2,966在该表BH FDR&lt;0.05（1,427左高、1,539右高），均是分群后发现阶段筛选，不是确认性验证。</p><p>C12高于C16的HIPK3、GSE1样本方向均10/10，log2FC约0.79、0.98；C16高的CAPG、RAB25、TMEM256亦10/10同向，但它们相对全体rest并不显著更高，其他群也高。C8/C12的vs-rest短名单富含EGR1、ATF3、FOS/JUN等信号。EGR1具有即时诱导应答性质，因而“应答较强”是合理的候选解释，尚不能据此区分体内状态与处理过程影响。<a href="https://pubmed.ncbi.nlm.nih.gov/8336701/">EGR1原始研究</a></p>{table(bpair)}<p>{link('02_pseudobulk_markers/model_summary.tsv')} · {link('03_marker_review/display_marker_evidence.tsv')}</p>
<h2>4. 来源与QC不能忽略</h2><p>C16有2,190/2,779（78.8%）来自HRA003620；NC-24与NC-27合计1,931（69.5%）。C12最大单一样本P01T占22.9%，C8最大占14.9%。这些组成使用完整FULL细胞分母，没有把抽样比例当丰度。取材Normal/NMIBC/MIBC不是逐细胞恶性身份，不推断进展方向。</p><p>全部10个C12/C16配对样本中C16的平均nCount和nFeature均高；C12−C16的样本差中位数分别−5,983.9、−861.9。percent.mt有8/10为C12更高，中位差+1.07个百分点。TMM及配对模型并不自动消除这些与状态耦合的检测率/质量差异，故C16高检出基因不能直接视为特异marker。本轮按授权未回归QC、删除样本/细胞或改群。</p>{table(qtab)}{table(composition)}<p>{link('03_marker_review/source_composition.tsv')} · {link('03_marker_review/QC_sample_paired_differences.tsv')} · {link('00_audit/paired_source_coverage.tsv')}</p>
<h2>5. 候选marker与反例</h2><p>固定规则的每项短名单最多20；展示每state最多3。C8建议EGR1/GADD45B/ELOVL5并看：EGR1相对C6、C3的log2FC约0.57、1.70，方向13/15、21/21；ELOVL5为1.38、1.17，方向15/15、19/21。GADD45B在r0.4 C13更高，不能当C8特异标记。ELOVL5是MP07/MP10共享输入基因，只作脂质支持。</p><p>C12用HIPK3/GSE1辅助区分C16，但它们相对rest效应较小且其他群可更高；ELOVL5不区分C12/C16。C16的NRARP/IDI1/LDLR是支持性候选：NRARP与C12近似，IDI1在C18更高，LDLR在C14更高；均不足以单基因定义C16。CAPG/RAB25/TMEM256不优先作特异标记。具体反例包括：C8相对C6的EGR1在CNP_P04T、PRJNA_P1T1并不更高；C8相对C3的ELOVL5在HRA_NC-04、HRA_NC-27并不更高。完整样本key与反例另表列出，没有挑选只支持的样本。</p>{table(marktab)}<p>{link('03_marker_review/display_markers.tsv')} · {link('03_marker_review/display_marker_counterexamples.tsv')} · {link('03_marker_review/all_shortlists.tsv')} · {link('03_marker_review/shortlist_alternative_cluster_specificity.tsv')}</p>
<h2>6. MP10非共享基因的功能支持</h2><p>MP07独有38、MP10独有38、共享12，与冻结列表一致。限定于本地GO带脂质/甾醇关键词的MP10独有基因，C8的IDI1和INSIG1均27/28样本高、LDLR26/28、HMGCR24/28、FDFT1及MVD22/28；因此信号不只来自共享基因。但DBI仅8/28、DHCR24仅7/28更高，不能称整条通路整齐激活。C16相对rest有较多同向支持，C12/C16之间并非所有脂质基因均能区分。</p><p>本轮称“MP10相关脂质/甾醇代谢候选状态”。RNA表达不足以证明代谢通量、蛋白阳性或脂肪酸代谢已验证亚群。GO关键词审查是有限注释，不是富集检验；原MP10/MP07分数及top50均未改动。</p><p>{link('03_marker_review/MP07_MP10_gene_partition.tsv')} · {link('03_marker_review/MP07_MP10_per_gene_summary.tsv')} · {link('03_marker_review/MP07_MP10_lipid_sterol_annotation.tsv')}</p>
<h2>7. 下一轮优先级与未能回答的部分</h2><p>优先完整C8的多基因独立验证，同时保留C12/C16对照，尤其检查响应基因、脂质支持基因与深度/取材条件的关系。C16虽有较强MP10支持，来源集中和QC耦合使其暂不适合作为唯一主推的新亚群；C12也不宜直接替代C8。当前更支持有内部差异的MP10程序，连续性与稳定离散状态仍不能区分。此次没有运行空间、CNV、预后、机制或蛋白验证，不将原NMF预后效应转写给新群。</p><p>全部技术模块已执行；受限的是独立患者重复、GSE及部分来源的合格配对、QC混杂与外部生物学验证。原usage、输入marker、非输入marker均来自同批RNA；C8/C12/C16重叠和不同resolution不是外部复现。</p><p>{link('candidate_decision.tsv')} · {link('METHODS.md')} · {link('DECISIONS.md')}</p>
<h2>8. 推荐先看的图</h2>'''
for f in figs:body+=f'<figure><a href="04_figures_by_dataset/{f}.pdf"><img src="04_figures_by_dataset/{f}.png" loading="lazy"></a><figcaption>{f}；点击打开PDF。配套plotdata与说明见图索引。</figcaption></figure>'
body+='<p>'+link('04_figures_by_dataset/figure_index.tsv','全部图、dataset视图及数据索引')+'</p><h2>9. 绝对路径与复现</h2><ul>'
for p in ['README.md','METHODS.md','DECISIONS.md','candidate_decision.tsv','00_audit/input_manifest.tsv','00_audit/comparison_eligibility_with_libraries.tsv','02_pseudobulk_markers/model_summary.tsv','03_marker_review/all_shortlists.tsv','04_figures_by_dataset/figure_index.tsv','code/run_all.sh','logs/sessionInfo_edgeR.txt']:
 body+='<li>'+link(p,str(R/p))+'</li>'
body+='</ul><p>方法参考：<a href="https://bioconductor.org/packages/edgeR/">edgeR</a>；<a href="https://bioconductor.org/books/3.17/OSCA.multisample/multi-sample-comparisons.html">OSCA</a>；<a href="https://www.nature.com/articles/s41467-021-25960-2">生物学重复与DE</a>。实际API/版本保存在本地审计目录。</p>'
(R/'REPORT.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MP10 follow-up v2</title><style>body{max-width:1200px;margin:32px auto;padding:0 24px;font:16px/1.7 system-ui,sans-serif;color:#202330}h1,h2{line-height:1.35}h2{margin-top:2em}table{border-collapse:collapse;width:100%;font-size:12px;display:block;overflow:auto}td,th{border:1px solid #ddd;padding:7px;text-align:left}th{background:#f4f5f9}img{max-width:100%;height:auto;background:white}code{overflow-wrap:anywhere}.key{padding:16px 22px;background:#f3f0f9;border-left:4px solid #a4489e}a{color:#39479d}figcaption{font-size:13px;overflow-wrap:anywhere}li{overflow-wrap:anywhere}</style><body>'+body+'</body></html>')
readme=f'''# MP10状态后续分析 v2

已完成六项固定比较、原usage衔接、sample配对edgeR QL、marker/来源/QC复核及冻结坐标出图。

工作目录：`{R}`

首先打开 [REPORT.html](REPORT.html)，再看 [candidate_decision.tsv](candidate_decision.tsv) 和 [METHODS.md](METHODS.md)。

|比较|合格配对sample|dataset|检验基因|
|---|---:|---:|---:|
'''
for m in models:readme+=f"|{m['comparison']}|{m['n_pairs']}|{m['n_datasets']}|{m['n_tested']}|\n"
readme+='''
C8原MP10 usage在20/28合格样本更高；C12/C16有表达差异，但C16在全部10对中nCount/nFeature更高，且69.5%来自HRA两个sample。优先保留完整C8作独立marker评估，C12/C16作并行状态对照，不能宣布已验证的新亚群。

展示候选：C8 EGR1/GADD45B/ELOVL5；C12 HIPK3/GSE1/ELOVL5；C16 NRARP/IDI1/LDLR。各自局限和样本反例见03_marker_review。MP10非共享IDI1/INSIG1/LDLR等支持脂质/甾醇相关解释，但不是整条通路一致激活。

原NMF/Seurat/UCell/labels/UMAP全部只读。FULL含Normal，未删细胞。患者映射未核实，不能称独立患者验证。GSE326225保留完整视图但六项比较无合格配对；无配对图的模块不伪造结果。未运行空间、CNV、预后或额外聚类。

输出目录：00_audit输入/资格，01_usage_bridge原usage，02_pseudobulk_markers完整DE/过滤/pseudobulk缓存，03_marker_review短名单/功能/来源/QC，04_figures_by_dataset透明PDF+300dpi PNG及实际绘图数据，code与logs。

复现入口：`bash code/run_all.sh`。代码固定绝对输入路径，配置一致的聚合/DE检查点可复用，原UCell和UMAP始终读取。审阅包为`mp10_state_followup_v2_review.zip`，排除原Seurat、完整表达/usage矩阵及大型pseudobulk缓存，包含完整DE压缩表、报告、代码、关键图及图索引。
'''
(R/'README.md').write_text(readme)
print('Reports and counterexample/library audit tables written')
