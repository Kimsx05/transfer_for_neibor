from pathlib import Path
import csv,gzip,html,json,math,datetime,collections
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4')
def read(f):
 p=R/f
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as h:return list(csv.DictReader(h,delimiter='\t'))
def write(f,rows):
 with (R/f).open('w') as h:
  w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def num(x):
 try:return float(x)
 except:return math.nan
def fmt(x):
 v=num(x)
 if not math.isfinite(v):return x
 return f'{v:.4g}'
def table(rows,cols=None):
 if not rows:return '<p>Not available; see status.</p>'
 cols=cols or list(rows[0]);return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in cols)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+html.escape(fmt(r.get(c,'')))+'</td>' for c in cols)+'</tr>' for r in rows)+'</tbody></table></div>'
def link(f,t):return '<a href="'+html.escape(f,quote=True)+'">'+html.escape(t)+'</a>'
def fig(f,caption):return '<figure>'+link(f+'.pdf','PDF')+'<img loading="lazy" src="'+f+'.png"/><figcaption>'+caption+'</figcaption></figure>'
models=read('02_rna_DEG/model_summary.tsv');pairs=read('01_score_comparison/paired_tests.tsv');scores=read('01_score_comparison/cluster_score_summaries.tsv');cover=read('03_sample_robustness/composition_coverage_all_levels.tsv');panel=read('06_supplement/fixed_panel_all_models.tsv');tops=read('06_supplement/all_models_top50.tsv');paths=read('06_supplement/focus_pathways_all_models.tsv');index=read('06_supplement/figure_index.tsv');elig=read('03_sample_robustness/model_sample_eligibility.tsv')
source={}
for m in models:
 q=[x for x in elig if x['branch']==m['branch'] and x['control']==m['control_cluster'] and x['eligible']=='TRUE'];c=collections.Counter(x['dataset'] for x in q);source[(m['branch'],m['control_cluster'])]='; '.join(f'{k}: {v} eligible pairs' for k,v in sorted(c.items()));m['eligible_dataset_pairs']=source[(m['branch'],m['control_cluster'])]
write('07_review/model_execution_and_sources.tsv',models)
# Supplement index with actual source composition and provenance in English captions.
for n,f in enumerate(index,1):
 f['supplement_id']=f'Supplementary Figure S{n}'
 for br in ['full','exclude_P07T']:
  for ctrl in ['rest','C14']:
   if f['figure']==f'{br}_{ctrl}_top20_DotPlot':
    m=next(x for x in models if x['branch']==br and x['control_cluster']==ctrl)
    f['caption']+=' Actual RNA model: '+m['origin']+'; '+m['n_pairs']+' eligible specimen pairs; '+source[(br,ctrl)]+'.'
 if 'fixed_panel_effects' in f['figure']:
  br='exclude_P07T' if f['figure'].startswith('exclude') else 'full'
  f['caption']+=' Model-specific sources: '+' | '.join(c+': '+s for (b,c),s in source.items() if b==br)+'.'
write('06_supplement/figure_index.tsv',index)
(R/'06_supplement/FIGURE_CAPTIONS.txt').write_text('\n\n'.join(f["supplement_id"]+' — '+f['figure']+'\n'+f['caption'] for f in index))
# Explicit full-versus-exclusion score/coverage comparison.
comparison=[]
for br in ['full','exclude_P07T']:
 c=next(x for x in cover if x['branch']==br and x['level']=='ALL');u=next(x for x in scores if x['branch']==br and x['cluster']=='C18' and x['score']=='UCell');n=next(x for x in scores if x['branch']==br and x['cluster']=='C18' and x['score']=='NMF')
 comparison.append(dict(branch=br,n_full=c['n_full'],n_samples=sum(x['branch']==br and x['level']=='sample' for x in cover),n_K32=c['n_K32'],UCell_mean=u['mean'],UCell_rank=u['rank'],usage_mean=n['mean'],usage_rank=n['rank'],NMF_precision=c['NMFprecision'],NMF_recall=c['NMFrecall'],UCell_precision=c['UCellprecision'],UCell_recall=c['UCellrecall']))
write('07_review/full_vs_exclude_score_coverage.tsv',comparison)
# Paired QC descriptive differences, never a covariate-adjusted re-clustering.
qc=read('03_sample_robustness/sample_QC.tsv');g=collections.defaultdict(dict)
for x in qc:g[(x['branch'],x['measure'],x['sample'],x['dataset'])][x['K32']]=x
qrows=[]
for (br,measure,s,ds),v in g.items():
 l=v.get('TRUE',{});r=v.get('FALSE',{});dif=num(l.get('value'))-num(r.get('value'));qrows.append(dict(branch=br,measure=measure,sample=s,dataset=ds,n_K32=l.get('n',0),n_rest=r.get('n',0),mean_K32=l.get('value','NA'),mean_rest=r.get('value','NA'),difference=dif if math.isfinite(dif) else 'NA',eligible_30_each=int(l.get('n',0))>=30 and int(r.get('n',0))>=30))
write('03_sample_robustness/QC_sample_differences.tsv',qrows)
qc_summary=[]
for br in ['full','exclude_P07T']:
 for measure in sorted({x['measure'] for x in qrows}):
  a=[x for x in qrows if x['branch']==br and x['measure']==measure and x['eligible_30_each'] and math.isfinite(num(x['difference']))];v=sorted(num(x['difference']) for x in a);qc_summary.append(dict(branch=br,measure=measure,n_pairs=len(a),n_K32_higher=sum(x>0 for x in v),mean_difference=sum(v)/len(v) if v else 'NA'))
write('03_sample_robustness/QC_direction_summary.tsv',qc_summary)
write('04_14MP_profiles/cluster_sizes_full_exclude.tsv',[{k:x[k] for k in ['branch','cluster','n','rank']} for x in scores if x['score']=='UCell'])
rank_summary=[]
for br in ['full','exclude_P07T']:
 for ctrl in ['rest','C17','C7','C14']:
  for size in [20,50]:
   a=[r for r in tops if r['branch']==br and r['control']==ctrl and int(r['up_rank'])<=size]
   if not a:continue
   rank_summary.append(dict(branch=br,control=ctrl,top_n=size,n_genes=len(a),n_input=sum(x['in_any_14MP_top50']=='TRUE' for x in a),n_noninput=sum(x['in_any_14MP_top50']=='FALSE' for x in a),n_mapped_symbol=sum(x['gene_symbol']!='NA' for x in a),n_with_lipid_name_pathway_annotation=sum(x['n_lipid_pathway_annotations']!='NA' for x in a),n_FDR_lt_05=sum(num(x['FDR'])<.05 for x in a),genes=';'.join(x['gene_id'] for x in a),annotation_caution='Frozen pathway membership only; not a comprehensive functional classification or new enrichment.'))
write('06_supplement/top20_top50_input_annotation_summary.tsv',rank_summary)
# Numerical GSEA/mapping limitations preserved for every model.
limits=[]
for m in models:
 fg=read(m['model_dir']+'/function/enrichment_all.tsv.gz');rk=read(m['model_dir']+'/rank_mapping.tsv.gz');cc=collections.Counter(x['mapping_reason'] for x in rk);limits.append(dict(branch=m['branch'],control=m['control_cluster'],n_pathways=len(fg),n_NA_NES=sum(x['NES']=='NA' for x in fg),n_NA_FDR=sum(x['padj']=='NA' for x in fg),n_at_p_floor=sum(num(x['pval'])<=.0001 for x in fg),rank_mapping_reasons=json.dumps(dict(cc),ensure_ascii=False),action='Fixed parameters retained; NA not called nonsignificant; nonfinite signed sqrt(F) excluded by common mapping rule.'))
write('00_audit/function_numerical_limitations.tsv',limits)
canon=['C5::GOBP_FATTY_ACID_BIOSYNTHETIC_PROCESS','C5::GOBP_FATTY_ACID_METABOLIC_PROCESS','C2::REACTOME_CHOLESTEROL_BIOSYNTHESIS','C5::GOBP_STEROL_BIOSYNTHETIC_PROCESS','C2::REACTOME_REGULATION_OF_CHOLESTEROL_BIOSYNTHESIS_BY_SREBP_SREBF','C5::GOBP_FATTY_ACID_BETA_OXIDATION'];pathsmall=[];seen=set()
for x in paths:
 k=(x['branch'],x['control'],x['pathway'])
 if x['pathway'] in canon and x['control'] in ['rest','C14'] and k not in seen:pathsmall.append(x);seen.add(k)
write('07_review/key_pathway_full_exclude.tsv',pathsmall)
htmlparts=['''<!doctype html><html lang="zh"><meta charset="utf-8"><title>K32 validation v4</title><style>body{font:16px/1.6 system-ui,sans-serif;max-width:1250px;margin:35px auto;padding:0 24px;color:#243142}h1,h2,h3{color:#173B75}table{border-collapse:collapse;font-size:13px}th,td{padding:6px 9px;border:1px solid #d9e0e7;text-align:left;white-space:nowrap}th{background:#eef3fa}.scroll{overflow:auto;max-height:620px}img{max-width:100%;background:white}figure{margin:24px 0}figcaption{font-size:13px;color:#475466}.note{background:#f0f4fa;padding:16px;border-left:4px solid #173B75}a{color:#174c89}</style><h1>K32 补充验证与 Supplement 审阅</h1>''',
'<p>冻结对象：K32=C18，FULL_Z_r1.00；123,699 个原上皮分析细胞、57 个生物样本、5 个 dataset、27 群。K32 3,282 细胞。Normal 22,436 细胞及原纳入 Cycling 均保留。本轮未重新拟合 NMF/UCell/聚类/UMAP。</p>',
'<div class="note"><b>主要结论：</b>K32 在两种 MP10 原始均值中均第一，排除最大贡献样本后仍第一。它是 MP10 与 MP07 共同升高的状态；RNA 支持 SREBP/部分甾醇相关特征，但对 C14 不支持“普遍增强的脂肪酸代谢”。保留 K32 工作标签，尚不能改名为 SREBF2⁺ epi 或代谢通量增强亚群。</div>',
'<h2>1. 分数排名、样本效应与统计</h2><p>总体 UCell：K32 0.359668 vs C14 0.244698；原 usage：28.6331 vs C7 13.7194。下表的差值是合格样本内均值差的等权平均，不是总体细胞加权差。所有比较均无 ties/零差，使用双侧精确 Wilcoxon；Holm 在每个分支的两项预定检验内校正。NMF 检验用原尺度，只有展示用 log1p。</p>',table(pairs),fig('01_score_comparison/full_UCell_27cluster_violin','Single-cell UCell distribution; sample-paired Holm p=0.0625, no significance star.'),fig('01_score_comparison/full_NMF_27cluster_violin','NMF raw usage displayed log1p. Mean diamonds are log1p(mean(raw usage)); Holm p=0.03125.'),
'<h2>2. RNA 模型与主要结果</h2><p>三项 v3 模型精确复用；新增 full C14、exclude_P07T rest/C14 三项均成功拟合。比较方向均为 K32−对照。模型样本数不是已确认独立患者数。</p>',table(models,['branch','control_cluster','origin','seed','n_pairs','n_datasets','n_tested','residual_df','status','eligible_dataset_pairs']),
'<p>对 C14，非输入基因 SREBF2、SREBF1、HSD17B7 的 log2FC 分别为 0.673、0.683、0.735，FDR 为 0.0216、0.0271、0.0205；五个样本的普通 log-normalized RNA 均值方向均为正。FASN 的 log2FC=0.374，FDR=0.185；INSIG1 FDR=0.245。SREBF2、SREBF1、HSD17B7 在 K32 的等权检出率反而低于 C14，说明表达效应不等于独占检出，不能把它们写成特异阳性标签。</p>',
table([x for x in panel if x['control']=='C14'],['branch','gene_requested','log2FC','FDR','pct_left','pct_right','positive_fraction','up_rank','in_any_14MP_top50']),
'<h2>3. 固定 top20/top50 与功能</h2><p>排序仅限定 log2FC&gt;0，然后按 FDR、PValue、gene_id 升序，不按脂质功能筛选。完整上/下调和 top50 含未注释转录本；输入属性、检测率和方向均保留。rest top20 有 13 个输入基因；C14 top20 有 8 个。C14 前列 ANLN、POLQ、CENPP、TPX2 等提示还应考虑细胞周期/分裂背景。完整富集同时出现染色质组织、RHO/CDC42 与形态发生等信号；带神经名称的共享基因通路不证明神经细胞身份。</p>',table(rank_summary,['branch','control','top_n','n_input','n_noninput','n_mapped_symbol','n_with_lipid_name_pathway_annotation','n_FDR_lt_05']),
'<p>“脂质注释”只表示命中冻结 GO/Reactome 中名称含 lipid/sterol/cholesterol/fatty acid/SREBP/SREBF 的集合，是可追溯注释回查，不是全功能分类比例或新富集。</p>',
fig('06_supplement/full_rest_top20_DotPlot','Complete positive-effect top20 vs rest. Equal eligible-sample RNA mean and detection; 12 pairs, four datasets; reused v3 model.'),fig('06_supplement/full_C14_top20_DotPlot','Complete positive-effect top20 vs C14. Five pairs from three datasets; newly fitted RNA model.'),
'<p>vs rest：胆固醇/甾醇合成、SREBP 调控及脂肪酸合成/代谢均正向显著；β氧化未显著。vs C14：SREBP 胆固醇调控 NES=1.895、FDR=0.00799；胆固醇合成 NES=1.638、FDR=0.0998，甾醇合成 NES=1.558、FDR=0.0992。脂肪酸合成/代谢与 β 氧化的 NES 为负且不显著。因此有部分 SREBP/甾醇相关 RNA 支持，不能推广为所有脂质途径都更高。</p>',table(pathsmall,['branch','control','pathway','NES','padj']),
'<h2>4. 排除 P07T 的敏感性</h2><p>精确排除 CNP0000460|CNP0000460_P07T 的 2,912 个细胞，包括 1,461/3,282（44.52%）K32；P07N 保留。剩余 120,787 个细胞、56 个样本、1,821 个 K32。全局 Z 和两种 high 定义完全不变。K32 两种 score 仍第一，但校正后分数比较均未达0.05。</p>',table(comparison),
'<p>重新拟合后，vs rest 的 SREBF2/SREBF1/FASN/HSD17B7 及主要脂质通路仍保留；vs C14 的 SREBF2/SREBF1/HSD17B7 仍正向显著，FASN 仍不显著。SREBP 胆固醇调控 NES=1.762、FDR=0.0378；胆固醇/甾醇合成仍未显著。该分析表明特征不完全由 P07T 决定，但排除后覆盖率降低，且 C14 仅剩4对，不能视为独立复制。</p>',
'<h2>5. 来源、患者映射与 QC</h2>',table([x for x in cover if x['branch']=='full' and x['level'] in ['dataset','source']],['level','unit','n_full','n_K32','contribution','K32_fraction','NMFprecision','NMFrecall','UCellprecision','UCellrecall']),
'<p>K32 在全部五个 dataset 中出现，但 GSE326225 仅4细胞，无合格 RNA 配对。最大 dataset 为 CNP0000460（53.14%），其次 HRA003620（33.30%）。Normal/NMIBC/MIBC 均有 K32；不能据此自动称为恶性亚群。全部57样本的零候选、零召回和 NA 分母均保留，覆盖表不使用30细胞过滤。</p>',table(qc_summary),
'<p>原 patient_recorded 由名称规则 derived/provisional 构建，缺少独立确认的真实映射，因此未进行患者级检验。按现有未核验键，rest 模型中 HRA003620::05 对应 BC-05 和 NC-05 两个合格样本；这提示可能存在重复来源结构，不能将12样本对写成12独立患者。分数比较的未核验患者结构另表列出。QC 差异只作描述，未用于删除或重新分群。</p>',
'<h2>6. 14MP 复合 profile 与冻结 UMAP</h2><p>K32 的 MP10 全局 Z 均值2.606、MP07为1.829，二者均为27群第一；MP03为0.904（第4）、MP11为0.695（第7）。排除 P07T 后 MP10=2.376、MP07=2.058，仍均第一。MP10 是最突出的相对偏离之一，但状态同时携带其他程序，不应描述为单一 MP10 开关。各 MP 原分数不代表不同功能的绝对强弱，Z也不证明机制主导。</p>',fig('04_14MP_profiles/full_frozen_Z_heatmap','Frozen cell-level global Z, then cluster means; no re-standardization of cluster means.'),
'<p>在冻结的14MP分数 UMAP上，K32集中于上方局部区域但存在少量离散细胞；在原 RNA UMAP 中，K32分散于多个区域。此观察支持使用“细胞状态”工作表述，不能把二维形状当作独立生物学验证。</p>',fig('05_umap/ALL/FULL_Z_ALL_triptych','Frozen MP-space coordinates; original raw score and usage, K32 highlight.'),fig('05_umap/ALL/RNA_ALL_triptych','Frozen original RNA coordinates; no optimization for a separate K32 island.'),
'<h2>7. 支持、限制与技术状态</h2><p>继续使用 K32 的依据是双分数高均值、较高 UCell-high purity、部分可复核非输入 RNA 效应，以及排除主要来源后保留的 SREBP 相关方向。限制是总体 recall 较低、来源集中、与 C14 配对数有限、部分检测率差反向、MP07共升高、RNA空间分散和患者映射不可靠。所有结果属于同一发现数据内候选刻画；单模型 BH-FDR 不控制先前候选搜索，也不是独立外部验证。维持 K32 标签，不自动改为 SREBF2⁺ epi，不等同于代谢通量。</p>',
'<p>六项 RNA 模型均完成，无依赖关键输入的阻塞。患者级敏感性未评价；部分 fgsea 结果为数值 NA，完整记录如下。使用固定参数，不以提高显著性为理由增加随机数或改对照。</p>',table(limits),
'<h2>优先审阅与完整 Supplement 索引</h2><p>'+link('06_supplement/fixed_panel_all_models.tsv','固定候选面板')+' · '+link('06_supplement/top20_top50_input_annotation_summary.tsv','完整 top20/top50 基因列表')+' · '+link('06_supplement/focus_pathways_all_models.tsv','所有重点通路（含反向、非显著及NA）')+' · '+link('06_supplement/FIGURE_CAPTIONS.txt','英文图注')+' · '+link('00_audit/INPUT_MANIFEST.tsv','输入清单')+' · '+link('METHODS.md','实际 Methods')+'</p>']
for m in models:
 htmlparts.append('<p>'+m['branch']+' vs '+m['control_cluster']+': '+link(m['model_dir']+'/DE_annotated.tsv.gz','完整 DEG')+' · '+link(m['model_dir']+'/top20.tsv','top20')+' · '+link(m['model_dir']+'/top50.tsv','top50')+' · '+link(m['model_dir']+'/function/enrichment_annotated.tsv.gz','完整富集及 leading-edge')+'</p>')
htmlparts.append('<ul>'+''.join('<li>'+f['supplement_id']+': '+link(f['pdf'],f['figure'])+' · '+link(f['png'],'PNG')+' · '+link(f['plot_data'],'data')+'</li>' for f in index)+'</ul></html>')
(R/'REPORT.html').write_text('\n'.join(htmlparts))
(R/'METHODS.md').write_text('''# Actual methods: K32 validation v4

K32 is the previously selected candidate C18 in frozen FULL_Z_r1.00. The original gene-expression-based NMF, top50 signatures, unweighted UCell scores, original FULL column standardization, 14-dimensional graph partition, and embeddings were reused as distinct historical steps. NMF usage is an untransformed program coefficient, not a pathway matrix or a maximum-coefficient cell label. No NMF, UCell, graph, clustering, PCA, Harmony, balanced sampling or UMAP was fitted in v4.

FULL contains 123699 epithelial-analysis cells, including the previously included Cycling and Normal cells, 57 dataset|orig.ident specimen keys and five datasets. All sources were joined by exact cell_id, with unique-ID, set, score, high and K32 member-hash checks. M_RAW and M_Z are frozen 14-column matrices; Z reconstruction was checked against the frozen mean_FULL/sd_FULL. Patient names are derived/provisional, not verified independent patient identifiers. Scores, labels and embeddings were read-only.

The v3 NMF-high cutoff is 16.820936842046201, raw usage >= cutoff; UCell-high is the v3 FULL type7 90th percentile 0.25970159376059676, >= cutoff with ties. Neither cutoff was recomputed in v4. P07T exclusion uses only CNP0000460|CNP0000460_P07T; all other frozen definitions remain unchanged. P/R denominators are the retained full population. Zero candidate gives precision NA, zero high gives recall NA, and positive high/no candidate gives recall0. Coverage does not use model eligibility filtering.

Single-cell violins retain all cells and zero values, scale='width', with raw-scale mean/median summaries. NMF violins display log1p values and log1p(mean(raw usage)) mean markers. Each score is independently ordered by its raw FULL cluster mean. Paired score tests use per-specimen raw means with >=30 cells on each side: C18 vs C14 for UCell, C18 vs C7 for raw usage. Two-sided exact Wilcoxon applies when absolute differences have no ties or zero differences; otherwise a normal approximation with continuity correction is explicitly recorded. All four realized tests were exact. Holm adjustment is separate across the two prespecified tests in each population. The candidate was previously selected, so these are exploratory discovery-data comparisons.

Existing sample-by-cluster original counts, sum of ordinary LogNormalize expression and detected-cell sums were reused after bin counts and membership checks. The v3 aggregate lacked C14. The original RNA object was read once; its single counts layer was subset sparsely to C14, checked nonnegative integer and aligned by ID; its data layer was checked against LogNormalize(scale.factor=10000). Only sparse C14 specimen aggregates were added. The original whole gene-by-cell matrix was never made dense. Rest aggregates were sample totals minus C18, with count non-negativity and conservation checks.

Three full models (rest, C17, C7) were reused exactly from v3 after candidate and control membership hashes and configuration checks. Full C14 and exclude_P07T rest/C14 were newly fitted. Exclusion models re-evaluated pair eligibility and repeated filtering/TMM/fitting/FDR. RNA eligibility requires >=30 cells per side, >=3 pairs, full-rank ~sample+group design and positive residual degrees of freedom. Right is the reference; groupLeft is K32 minus control. edgeR filterByExpr uses group, min.count=10, min.total.count=15, large.n=10, min.prop=.7; library sizes are recomputed after filtering, TMM normalization, estimateDisp robust=TRUE, glmQLFit robust=TRUE legacy=FALSE, glmQLFTest(groupLeft), BH over all tested genes. Counts are model input; expression and detection summaries are descriptive equal-weight means across those same eligible samples. Positive_fraction describes the direction of ordinary log-normalized sample means, not a count-model coefficient per sample. Dataset-specific directions are retained.

Top20/top50: log2FC>0 followed by ascending FDR, PValue, gene_id. No significance, detection or lipid prefilter. Fixed panel is SREBF2,SREBF1,FASN,HSD17B7,ELOVL5,ADIPOR2,HMGCS1,ABHD3,INSIG1,ACSL3,ACSS2. Input membership is from actual frozen signatures; non-input is not external validation.

Frozen human GO:BP and Reactome sets from v3 (9377 sets) and v3 annotation are used. Exact source hash/release caveat is in 00_audit/resource_provenance.txt. Full tested-gene ranking is sign(log2FC)*sqrt(F). Unmapped/nonfinite values are excluded with reasons; duplicate symbols retain highest logCPM then gene_id. fgseaMultilevel: sampleSize101,nPermSimple1000,minSize10,maxSize500,eps1e-4,nproc1,scoreType='std'. New stochastic work uses seed123; reused v3 fgsea retains seed42. All positive/negative pathways, leading edges, input memberships, NA estimates and p-floor limitations are saved. No ORA. Lipid-name pathway membership and leading-edge/top-list intersections are annotations, not new tests or complete functional classification.

Frozen FULL_Z and original RNA embeddings use identical ID sets across panels. All five dataset views retain FULL axes and score limits. Nonnegative continuous values use plasma; exact zeros black, missing cyan cross (none observed); positive points are ordered ascending score. Grey background is drawn before K32 #173B75 at equal point size. PNG is transparent 300dpi and PDF transparent with rasterized point layers. Heatmaps aggregate cell-level raw UCell/frozen Z, without cluster-mean re-scaling. K32’s MP10 and MP07 ranks are descriptive, not measures of metabolic flux or mechanistic dominance.

All six models, analyses and figures are discovery-data characterization. Within-model BH does not correct the prior resolution/candidate search. No new spatial, CNV, prognosis, mechanism or functional naming analysis was performed. Heavy steps execute sequentially with at most two data/qs threads and one BLAS/fgsea thread. Resume scripts validate configurations and input identity. Session information, input hashes, model design/filter/eligibility and model status are saved.
''')
(R/'README.md').write_text('''# K32 validation v4

Open REPORT.html first. K32=C18 of frozen FULL_Z_r1.00 (3282/123699 cells).

- 01_score_comparison: all27 cluster summaries/violins, raw-scale sample-paired tests.
- 02_rna_DEG: full rest/C17/C7 reused; full C14 and exclude_P07T rest/C14 newly fitted. Complete DEG, ordered top20/top50 and full GSEA in every model directory.
- 03_sample_robustness: all specimens/source coverage, eligibility, QC and fixed-panel sample expression.
- 04_14MP_profiles: raw and frozen-global-Z means, ranks and sizes.
- 05_umap: frozen FULL_Z/RNA, ALL and every dataset; per-cell plot data.
- 06_supplement: figure index, English captions, fixed panel and functional summaries.
- 07_review: concise full/exclusion and model-source tables.

Execute code/run_all.sh to resume with matching inputs/configuration. Do not source v3 scripts. Historical source objects and large aggregate caches are not part of the review ZIP; paths/checksums are recorded in 00_audit/INPUT_MANIFEST.tsv. All per-cell UMAP data and complete DEG/enrichment tables needed for review are included.
''')
(R/'STATUS.md').write_text('''# Execution status

All requested score, six RNA-model, source, exclusion-sensitivity, 14MP and frozen-embedding modules completed. Three historical RNA/fgsea results reused; three new models fitted with seed123. No key-input blocked module.

Not evaluable: patient-unit sensitivity because patient linkage is derived/provisional only. Some fgsea pathways have numerical NA under fixed parameters; see 00_audit/function_numerical_limitations.tsv. These are neither negative tests nor missing whole analyses. No parameter adaptation to gain significance.

Final file/integrity acceptance and package status are in FINAL_ACCEPTANCE.tsv and COMPLETED.json after packaging. No downstream analyses are authorized by this completion record.
''')
print('Report and Supplement captions written:',len(index),'figures')
