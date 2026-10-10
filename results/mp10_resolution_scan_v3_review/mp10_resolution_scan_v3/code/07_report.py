from pathlib import Path
import csv,gzip,html,collections,statistics,datetime,json
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3')
def rd(p):
 p=R/p
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def wt(p,rows):
 if not rows:return
 with (R/p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def tab(rows,cols=None):
 if not rows:return '<p>无可评估结果；见状态表。</p>'
 cols=cols or list(rows[0]);return '<table><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in cols)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(r.get(c,'')))+'</td>' for c in cols)+'</tr>' for r in rows)+'</table>'
def link(p,t=None):return f'<a href="{p}">{html.escape(t or p)}</a>'
def pct(x):
 try:return f'{100*float(x):.2f}%'
 except:return x
def nfmt(x):
 try:return f'{float(x):.3g}'
 except:return x
can=rd('08_review/FINAL_CANDIDATES.tsv');rec=[x for x in can if x['recommended']=='TRUE'];models=rd('04_pseudobulk_markers/model_summary.tsv');met=rd('02_coverage/all_clusters_two_high.tsv');sz=rd('01_partitions/all_cluster_sizes.tsv');eq=rd('02_coverage/sample_equal_summary.tsv');improve=rd('02_coverage/source_improvement_summary.tsv');fg=rd('08_review/function_family_summary.tsv');fgall=rd('08_review/function_top_all_pathways.tsv');markers=rd('08_review/display_marker_evidence.tsv');directions=rd('08_review/display_marker_sample_directions.tsv.gz');counterexamples=[x for x in directions if float(x['delta'])<=0];wt('08_review/display_marker_counterexamples.tsv',counterexamples);stability=rd('06_stability/candidate_Jaccard.tsv');baseline=rd('02_coverage/old_baselines.tsv');qct=rd('08_review/paired_QC_summary.tsv');source=rd('02_coverage/candidate_source_coverage.tsv');elig=rd('04_pseudobulk_markers/all_comparison_eligibility.tsv')
# 24 partition overview with both champions and their four P/R values.
over=[]
for p in dict.fromkeys(s['partition'] for s in sz):
 z=[s for s in sz if s['partition']==p];cs=[c for c in can if c['partition']==p];row={'partition':p,'n_clusters':len(z),'min_cluster':min(int(x['n_cells']) for x in z),'max_cluster':max(int(x['n_cells']) for x in z)}
 for h in ['NMF','UCell']:
  c=next(x for x in cs if h in x['won_for'].split(';'));row[h+'_winner']=c['result_id'];row[h+'_candidate_id']=c['candidate_id']
  for k in ['precision_NMF','recall_NMF','precision_UCell','recall_UCell']:row[h+'_'+k]=c[k]
 over.append(row)
wt('08_review/24_partition_overview.tsv',over)
# Source-qualified plot status for every candidate/contrast/dataset, including omitted panels.
scopes=['ALL']+sorted({x['dataset'] for x in elig});avail=[]
for m in models:
 for ds in scopes:
  e=[x for x in elig if x['comparison_id']==m['comparison_id'] and x['eligible']=='TRUE' and (ds=='ALL' or x['dataset']==ds)];c=next(x for x in can if x['candidate_id']==m['candidate_id']);avail.append(dict(candidate_id=m['candidate_id'],comparison_id=m['comparison_id'],scope=ds,n_pairs=len(e),model_status=m['status'],display_markers=c['display_markers'],plot_status='available' if e and c['display_markers'] else 'no_eligible_source_pairs' if not e else 'no_rest_and_local_supported_display_gene'))
wt('08_review/marker_plot_availability.tsv',avail)
summary=[]
for b in [x for x in baseline if x['baseline']=='C8']:
 pass
bn=next(x for x in baseline if x['baseline']=='C8' and x['high']=='NMF');bu=next(x for x in baseline if x['baseline']=='C8' and x['high']=='UCell')
summary.append({'候选':'旧 FULL_Z r0.4 C8','细胞':bn['n_cluster'],'NMF P/R/F1':f"{pct(bn['precision'])} / {pct(bn['recall'])} / {nfmt(bn['F1'])}",'UCell P/R/F1':f"{pct(bu['precision'])} / {pct(bu['recall'])} / {nfmt(bu['F1'])}",'选择':'固定基线','rest/local pairs':'v2 28 / 按对照','展示marker':'见v2'})
for c in rec:summary.append({'候选':c['candidate_id']+' '+c['result_id'],'细胞':c['n_cells'],'NMF P/R/F1':f"{pct(c['precision_NMF'])} / {pct(c['recall_NMF'])} / {nfmt(c['F1_NMF'])}",'UCell P/R/F1':f"{pct(c['precision_UCell'])} / {pct(c['recall_UCell'])} / {nfmt(c['F1_UCell'])}",'选择':c['won_for']+'；全局F1冠军','rest/local pairs':c['rest_pairs']+' / '+c['local_pairs'],'展示marker':c['display_markers'] or '无符合展示规则的基因'})
wt('08_review/core_candidate_comparison.tsv',summary)
counts=collections.Counter(x['status'] for x in models);unique_models=len({x['model_key'] for x in models});notest=[{k:x[k] for k in ['comparison_id','partition','cluster','control_cluster','n_pairs','status','reason']} for x in models if x['status']!='COMPLETE']
# No fixed biological pass threshold; summarize 50/60/70 references descriptively only.
refs=[]
for h in ['NMF','UCell']:
 for t in [.5,.6,.7]:refs.append(dict(high=h,reference=t,n_single_clusters_both_ge=sum(float(x['precision'])>=t and float(x['recall'])>=t for x in met if x['high']==h)))
wt('08_review/PR_reference_counts.tsv',refs)
methods='''# v3方法与冻结规则

## 输入与范围
原NMF output_full_cc_seed5_v2的123699个上皮分析细胞，包含Normal及原已纳入的Cycling，不是TME atlas。沿用57个dataset|orig.ident biological samples、5 datasets、22436 Normal。原Seurat、NMF、MP top50、usage、v1/v2标签/降维只读；没有重跑NMF、UCell、PCA、Harmony、UMAP或BAL。患者映射仍derived/provisional，不声称独立患者重复。精确文件/版本、大小/mtime、cell ID、软件版本、图/分数校验见00_audit。

## 分区与high
M_RAW与M_Z按cell ID精确核验，14列固定。Z与FULL全局mean/sample-SD显式列运算一致。两套原SNN复用，原矩阵SHA256核对；k20、annoy/50trees、Euclidean、l2.normFALSE、prune1/15。Seurat algorithm4、leidenbase、RBConfigurationVertexPartition、modularity.fxn1、n.iter10、group.singletonsTRUE、seed42；固定12个resolution×2表示，每次前设seed。原RAW r.4和Z r.2/.4/.6都按成员集合复现，允许群号置换。只有resolution改变。

在新扫描前冻结high逐细胞逻辑及精确二进制阈值：NMF使用原usage_raw>=保存的16.820936842046201，共12370；UCell使用FULL原MP10 quantile(type7,p.90)=0.25970159376059676，>=并列纳入，共12372，交集7002。这是v3新定义的UCell-high，非v1/v2旧规则。旧C8 TP2014/3937核验通过。q80/q95只对预先确定的全局F1冠军描述，不进入主排行榜。

## 覆盖、候选与对照
每个非空单群计TP/FP/FN/TN，P=TP/n_cluster，R=TP/n_high，F1=2TP/(n_cluster+n_high)，enrichment=P/(n_high/n_FULL)。每分区分别选两种high F1冠军，相同者合并；排序F1、min(P,R)、P降序，最后canonical数字群号升序。全部群和两种high的分区内/全局P-R Pareto前沿保留。全局推荐按同样指标，最终representation/resolution确定性排序，在DE之前冻结；不按marker替换候选，不合并群冲榜。

局部竞争群从对应原SNN计算sum(SNN[C,K])/(n_C*n_K)最大者，排除自身、canonical群号破同分。所有候选均使用同一图规则，没有逐候选切换质心。另保存NMF/UCell最大群外TP承载群，并按冠军所属high（双冠军优先NMF）最多加一个不同于最近邻的carrier对照。对照资格不足不更换。

sample/dataset/source覆盖使用完整来源池，无30细胞过滤。high>0且候选=0的R为0；high=0的R为NA；候选=0的P为NA。pooled和sample等权分开，NA/零召回分母全保存。来源占候选比例与候选占来源比例分别输出。

## RNA marker与资源
原RNA90011×123699稀疏counts重新核验非负整数；已有data与普通LogNormalize(count/library×10000)核对，展示与计数模型分开。原Seurat只加载一次，按FULL sample总量和每partition所需sample×cluster顺序稀疏聚合；rest由同sample总量减candidate，核验counts非负及守恒。没有dense全基因×cell矩阵。线程/worker≤2；重型图、聚合、模型顺序执行，BLAS1/fgsea1；估计RNA峰值约30–40GiB、输出缓存约数GiB，实际资源见00_audit及报告。

模型与v2一致：每侧≥30细胞的完整sample对，至少3对，design~sample+group，Right reference、groupLeft正向；满秩/剩余自由度/非零库核验。filterByExpr(group=group,min.count10,min.total.count15,large.n10,min.prop.7)，TMM，estimateDisp robustTRUE，glmQLFit robustTRUE legacyFALSE，glmQLFTest。每比较所有通过过滤基因单独BH。NOT_TESTABLE保留原因、资格和描述；不是无差异。完整基因/设计/效应/检测率/方向/输入成员关系均保存。只按两侧完整成员hash及完全相同模型配置复用；同细胞集意味着同sample资格。Jaccard高不复用DE。

短名单：FDR<.05、log2FC≥.5、目标样本等权检测率≥.1、表达正方向≥.7；按dataset方向、sample方向、检测率差、logFC、ID排序，至多20。展示最多3，优先有已知gene_symbol、同时通过rest与固定最近竞争群门槛的非输入基因，兼顾已冻结功能类别和不同基因名前缀；不足不补。局部NOT_TESTABLE时若展示rest基因必须标为仅rest支持。额外carrier的完整效应保留，不能因此改冠军。基因前缀只是展示去冗余辅助，非正式基因家族注释。最终关键候选与旧C8还在真正共同合格sample上汇报描述性表达方向。

## 统一功能分析
复用冻结目录11A_msigdb_gene_sets_cache.qs，msigdbr26.1.0缓存签名7d3852357e591330；提取9377个人类GO:BP/Reactome集合。缓存未独立记录数据库发行版本，故按精确文件SHA256锁定，不宣称额外版本。排名sign(log2FC)*sqrt(F)，保存每个受检gene_id；映射到已知SYMBOL/唯一ENSEMBL符号，非有限/未映射记录排除，重复符号统一取logCPM最高、gene_id破同分，不选最大效应。fgseaMultilevel sampleSize101、nPermSimple1000、minSize10、maxSize500、eps1e-4、nproc1、seed42、scoreType std。P值低于1e-4不继续追求极小值，NA保留；每比较所有两库通路一起BH，非全候选搜索统一FDR。全表NES/FDR/leading edge及输入基因贡献保存。未做ORA。

展示固定区分fatty_acid、cholesterol_sterol、other_lipid、response_stress、proliferation、epithelial六类通路名关键词；正/负NES分别展示、每类最多3，类别可重叠，不能相加计独立通路数。完整9377集合的实际受检富集表一直保留。所有marker、富集和非输入基因都是同批发现RNA证据，不是外部验证；MP10参与建群/high选择带来定义性描述及选择偏倚。

## 图与参考
原C8重复模型测试基因集与v2完全相同（15902基因），log2FC最大差1.42e-9、FDR最大差9.06e-9，表达与检测率一致；微差来自行排序下的浮点计算。

同一分支共享原UMAP连续底图，RNA用已冻结umap_allgene；dataset取全局坐标子集，不重拟合。全部5dataset保留。非负连续plasma，精确0黑、NA灰，FULL范围固定；usage仅显示log1p，统计始终raw。RNA DotPlot是sample等权mean LogNormalize/detection，非usage；无网格。每图透明PDF和300dpi PNG，容器背景也核验；图数据/来源索引保存。成员流向的纵坐标为群大小堆积布局，不是表达或UMAP距离。

- [fgsea官方教程](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [edgeR官方入口](https://bioconductor.org/packages/edgeR/)，执行API以本地保存help为准。
- [Seurat FindClusters](https://satijalab.org/seurat/reference/findclusters)
'''
(R/'METHODS.md').write_text(methods)
body=f'<h1>MP10 FULL分辨率扫描 v3</h1><p>完成日期：{datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).date()}；目录：<code>{R}</code></p><div class="key"><p><b>主要结果：</b>NMF-high最优单群为FULL_RAW r0.40 C4（K07），UCell-high最优为FULL_Z r0.80 C13（K31）。RAW与Z的优劣依赖high定义。前者P/R约28.0%/28.5%；后者UCell P/R约62.5%/31.6%。没有扫描单群同时跨过50% P和50% R参考线，但该线不是生物学通过门槛。</p><p>有可解释的候选并不等于覆盖了完整MP10程序。来源偏倚、局部marker资格/功效以及成员变化仍限制单一亚群解释；不因UMAP岛形状下结论。</p></div>'
body+='<h2>1. 核心候选和旧C8</h2>'+tab(summary)+f'<p>24套分区全部完成，38个冻结候选、108个固定对照（{unique_models}个精确模型键）。模型状态：'+html.escape(str(dict(counts)))+'。'+link('08_review/24_partition_overview.tsv')+' · '+link('08_review/FINAL_CANDIDATES.tsv')+'</p>'
body+='<h2>2. 相对C8：改善与代价</h2><p>K07的NMF P、R均高于旧C8；但UCell P下降约11.75个百分点，R仅提高约2.77个百分点，UCell F1下降。K31的NMF P、R亦提高；UCell P提高约16.76个百分点，但R下降约0.24个百分点（3937降至3907个UCell-high细胞）。因此不能把F1提高描述成两项无代价改善，也不能把两名候选并集混入单群排行榜。</p><p>全群及未获冠军的Pareto群均保留：'+link('02_coverage/all_clusters_two_high.tsv')+' · '+link('02_coverage/PR_Pareto_fronts.tsv')+' · '+link('02_coverage/old_baselines.tsv')+'</p>'
body+='<h2>3. 来源覆盖与零召回</h2>'
for c in rec:
 k=c['candidate_id'];es=[x for x in eq if x['candidate_id']==k];ist=[x for x in improve if x['candidate_id']==k];body+=f'<h3>{k} {c["result_id"]}</h3><p>涉及{c["n_samples_present"]}个sample、{c["n_datasets_present"]}个dataset。最大sample {html.escape(c["max_sample"])}占{pct(c["max_sample_contribution"])}；最大dataset {c["max_dataset"]}占{pct(c["max_dataset_contribution"])}。来源构成：{c["source_composition"]}。</p>'+tab(es,['high','n_P_valid','n_P_NA','n_R_valid','n_R_NA','n_high_zero_recall','P_sample_equal','R_sample_equal'])+tab(ist,['field','high','n_valid','n_P_up','n_R_up','n_joint_up','n_NA'])
body+='<p>K07在NMF-high的53个有high样本中17个零召回，K31在UCell-high的50个有high样本中7个零召回。不能仅以“出现于多个dataset”称普遍改善。NA分母与零召回均未被marker门槛过滤。'+link('02_coverage/all_candidate_source_vs_C8.tsv')+' · '+link('02_coverage/candidate_source_coverage.tsv')+'</p>'
body+='<h2>4. Marker、功能和QC</h2>'
for c in rec:
 k=c['candidate_id'];mm=[x for x in models if x['candidate_id']==k];gg=[x for x in markers if x['candidate_id']==k];ff=[x for x in fg if x['candidate_id']==k and x['direction']=='positive' and x['control_type'] in ['rest','local'] and float(x['padj'])<.05];ff=sorted(ff,key=lambda x:(x['control_type'],float(x['padj'])))[:18];qq=[x for x in qct if x['candidate_id']==k];otherfun=[x for x in fgall if x['candidate_id']==k and x['direction']=='positive' and x['control_type']=='rest'][:6]
 body+=f'<h3>{k}</h3><p>展示marker：{c["display_markers"] or "无符合展示规则基因"}；非输入marker：{c["non_input_markers"] or "无"}。状态：{c["marker_status"]}。固定最近竞争群为{c["local_cluster"]}，额外carrier为{c["extra_cluster"]}。</p>'+tab(mm,['comparison_id','control_type','control_cluster','n_pairs','n_datasets','status','reason'])+tab(gg,['gene_id','log2FC','FDR','pct_left','pct_right','positive_fraction','local_log2FC','local_FDR','carrier_log2FC','carrier_FDR','in_any_14MP_top50','function_family','selection_reason'])+tab(ff,['control_type','function_family','pathway','NES','padj','leading_edge_input_n','leading_edge_n'])+'<p>不限功能关键词的vs-rest正向前列通路：</p>'+tab(otherfun,['pathway','NES','padj','leading_edge_input_n','leading_edge_n'])+tab(qq,['comparison_id','n_pairs','nCount_positive','nFeature_positive','nCount_median_delta','nFeature_median_delta','mt_median_delta'])
body+='<p>K31的IFRD1、NEDD9、EZH2均不属于14MP输入基因；vs-rest在24/24样本表达更高。最近竞争群为同分区C18：log2FC分别约0.679、0.858、0.774，FDR约0.000313、0.000275、0.00116，方向分别10/11、11/11、9/11。相对额外carrier C8也保持正效应，仍是同批发现RNA而非独立验证。</p><p>K31前列功能是转录/应答（如NGF-stimulated transcription、cellular response to heat），并有sterol biosynthetic process支持（NES约1.601，FDR约0.0448）；脂肪酸相关正向条目未达到FDR0.05。GO中出现免疫/分化等名称可能由共享信号基因贡献，不能据此把冻结上皮细胞改注释成T细胞或其他细胞类型；应结合leading edge审阅。EZH2入选也不等于证明该群是增殖亚群。</p><p>K07没有基因同时满足rest与最近局部对照的展示门槛；其vs-rest前列正向通路包括翻译与氧化磷酸化，不能只强调脂质条目。K07最近竞争群只有3对，局部区分的功效有限；其vs-rest的23对中nCount/nFeature全部更高。K31相对最近竞争群11对中nCount/nFeature全部更低。TMM与配对设计不保证这些与状态相关的检测差异已消失。展示基因必须结合完整效应、检测率差、方向和carrier表审阅，不能仅凭P值或RNA检出判断特异。</p><p>若局部表没有满足门槛的基因，意味着本次规则下缺少足够证据，不等于没有表达差异。通路名须区分脂肪酸、胆固醇/甾醇、其他脂质、应答/应激、增殖和上皮背景；FDR显著类别不能直接升级为机制。'+link('08_review/all_shortlists.tsv')+' · '+link('08_review/function_family_summary.tsv')+' · '+link('08_review/marker_plot_availability.tsv')+' · '+link('08_review/display_marker_counterexamples.tsv')+'</p>'
body+='<h2>5. 成员稳定性与共享sample表达</h2>'
for c in rec:
 k=c['candidate_id'];jj=[x for x in stability if k in [x['candidate_A'],x['candidate_B']]];jj=sorted(jj,key=lambda x:-float(x['Jaccard']));bestj=[]
 for typ in ['adjacent_resolution','RAW_Z']:bestj += [x for x in jj if x['type']==typ][:3]
 body+=f'<h3>{k}</h3>'+tab(bestj)
body+='<p>K07与旧C8 Jaccard约0.286；K31与旧C8约0.436，且K31主要与旧r0.6 C12重叠（Jaccard约0.628），与C16几乎不重叠。群名相同或resolution相同不是成员一致性证据。以上成员变化与不同high冠军分离，更符合MP10跨多个状态分布、边界受尺度/分辨率影响的解释；当前尚不足以确立覆盖充分且稳定的单一亚群。</p><p>'+link('06_stability/adjacent_cluster_flows.tsv')+' · '+link('06_stability/candidate_Jaccard.tsv')+' · '+link('06_stability/candidate_Jaccard_and_profile_concordance.tsv')+' · '+link('06_stability/candidate_pair_shared_sample_marker_directions.tsv')+' · '+link('08_review/shared_eligible_samples.tsv')+' · '+link('08_review/shared_sample_marker_directions.tsv')+'</p>'
shared=rd('08_review/shared_sample_marker_directions.tsv') if (R/'08_review/shared_sample_marker_directions.tsv').exists() else []
if shared:body+='<p>关键候选和旧C8的共同合格sample数：'+shared[0]['n_common_samples']+'。仅在这些真实共同sample上重述表达方向，未重拟合/挑样本提升FDR。</p>'+tab([x for x in shared if x['gene_id'] in set(';'.join(c['display_markers'] for c in rec).split(';'))])
outbytes=sum(p.stat().st_size for p in R.rglob('*') if p.is_file() and p.suffix!='.zip');heaprows=[line.split() for line in (R/'logs/03_aggregate.log').read_text().splitlines() if line.startswith(('Ncells','Vcells'))];heap=sum(max([float(x[-1]) for x in heaprows if x[0]==typ] or [0]) for typ in ['Ncells','Vcells'])
body+=f'<p>实际输出（含必要缓存、打包前）：{outbytes/1024**3:.2f} GiB；RNA聚合日志记录的R分配堆峰值约{heap/1024:.2f} GiB（不是OS RSS）。模型重型步骤顺序执行，BLAS/fgsea各1线程，原RNA只加载一次。</p>'
body+='<h2>6. 敏感性、技术完成与限制</h2>'+tab(rd('08_review/recommended_q80_q95_descriptive.tsv'),['candidate_id','quantile','cutoff','n_high','precision','recall','F1'])+'<p>q80/q95仅描述固定推荐候选，未参与排名或更换候选。原患者映射未独立确认，样本重复不等于独立患者；未新增空间、CNV、预后、机制或外部验证。所有DE、富集及非输入marker均属于同批发现RNA；单比较BH不是整个分辨率搜索的错误率控制。</p>'+tab(notest)+'<p>'+link('04_pseudobulk_markers/model_summary.tsv')+' · '+link('00_audit/old_partition_reproduction.tsv')+' · '+link('FINAL_ACCEPTANCE.tsv')+'</p>'
body+='<h2>7. 推荐审阅图</h2>'
figs=['ALL__NMF__PR_all_clusters','ALL__UCell__PR_all_clusters','ALL__sample_precision_recall','ALL__FULL_RAW__candidate_marker_DotPlot','ALL__FULL_Z__candidate_marker_DotPlot','ALL__adjacent_resolution_member_flows','ALL__RAW_Z_candidate_Jaccard','ALL__RNA__candidate_highlights']
for f in figs:body+=f'<figure><a href="07_figures/{f}.pdf"><img src="07_figures/{f}.png" loading="lazy"></a><figcaption>{f}；完整图注和实际数据见图索引，点击打开PDF。</figcaption></figure>'
body+='<p>'+link('07_figures/figure_index.tsv')+'</p><h2>输出绝对路径</h2><ul>'
for p in ['08_review/FINAL_CANDIDATES.tsv','08_review/RECOMMENDED_CANDIDATES.tsv','08_review/all_shortlists.tsv','04_pseudobulk_markers/model_summary.tsv','05_function/resource_provenance.txt','00_audit/INPUT_MANIFEST.tsv','CONFIG.yaml','METHODS.md','DECISIONS.md','FINAL_ACCEPTANCE.tsv']:body+='<li>'+link(p,str(R/p))+'</li>'
body+='</ul>'
(R/'REPORT.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>MP10 resolution scan v3</title><style>body{max-width:1280px;margin:30px auto;padding:0 22px;font:16px/1.7 system-ui,sans-serif;color:#202332}table{display:block;overflow:auto;border-collapse:collapse;font-size:12px;margin:15px 0}th,td{padding:7px;border:1px solid #ddd;max-width:480px;overflow-wrap:anywhere}th{background:#f4f5f9}h2{margin-top:2em}.key{background:#f4eff8;border-left:4px solid #a43777;padding:12px 20px}a{color:#344a9e}img{width:100%;height:auto;background:white}figcaption,code,li{overflow-wrap:anywhere}</style><body>'+body+'</body></html>')
dec='''# 决策与停止记录

本轮全部候选按固定覆盖指标在DE前选择，marker/富集没有改变冠军。保留RAW r0.40 C4（NMF-high F1冠军）与Z r0.80 C13（UCell-high F1冠军）并列供人工讨论，不把两者合并成“更完整单群”。

当前不支持已经找到高precision且高recall、跨来源稳定的完整MP10单一亚群。RAW改善NMF覆盖但有明显来源集中与深度关联；Z提升UCell纯度但召回未提升，也有来源和局部QC限制。K31的IFRD1/NEDD9/EZH2作为非输入表达组合可优先人工审阅，但UCell-high召回仅31.6%，不得视为完整MP10群。K07无同时通过rest和最近局部门槛的展示marker。K31功能更偏转录/应答伴部分甾醇支持，不能称脂肪酸代谢已验证亚群。具体局部对照与不足以最终表为准；无合格pair是不可检验，不写成无差异。不能按DE基因数量评优。

本轮停止于单细胞审阅。未运行空间、CNV、预后；未增加resolution/seed、修改基因集/high定义、合并群冲榜、删除不利样本。若下一轮希望加密resolution，只能另行讨论，当前不追加。优先讨论high定义是否对应研究问题、来源集中/质量耦合以及marker的独立验证设计，不根据UMAP外观追求漂亮群。

实现修复限于兼容性：YAML默认小数精度将1/15写为0.0666667，已改为17位精度；实际聚类一直使用1/15且旧分区复现。其他修复和状态见logs，未改冻结输入。
'''
(R/'DECISIONS.md').write_text(dec)
readme=f'''# MP10 resolution scan v3

首先打开 [REPORT.html](REPORT.html)。目录：`{R}`。

24套固定FULL_RAW/FULL_Z分区，38个固定候选，108个固定对照（{unique_models}个精确模型键）。状态：{dict(counts)}。旧分区成员复现通过。两种high及候选/局部对照均在DE前冻结。

|候选|细胞|NMF P/R|UCell P/R|展示marker|
|---|---:|---|---|---|
'''
for c in rec:readme+=f"|{c['result_id']}|{c['n_cells']}|{pct(c['precision_NMF'])}/{pct(c['recall_NMF'])}|{pct(c['precision_UCell'])}/{pct(c['recall_UCell'])}|{c['display_markers'] or '无符合规则的基因'}|\n"
readme+='''
没有保证找到“漂亮群”。结果显示high定义依赖、来源集中、零召回和成员边界变化；所有marker/功能均是同批发现RNA，不是外部验证。样本不是已确认独立患者。

主要表：[全部候选](08_review/FINAL_CANDIDATES.tsv)、[24套概览](08_review/24_partition_overview.tsv)、[模型状态](04_pseudobulk_markers/model_summary.tsv)、[完整短名单](08_review/all_shortlists.tsv)、[功能摘要](08_review/function_family_summary.tsv)、[图索引](07_figures/figure_index.tsv)。完整DE和过滤在04_pseudobulk_markers/models，每个比较通过comparison_model_index.tsv追溯；功能全表按同model_key在05_function。

复现：code/01_audit_partitions.R后执行code/run_remaining.sh；出图/报告/验收脚本见code。相同配置/精确成员的缓存可复用。原矩阵/Seurat未打包，轻量审阅ZIP包含完整DE/功能表、短表、图、代码和说明。
'''
(R/'README.md').write_text(readme)
print('Reports generated',dict(counts))
