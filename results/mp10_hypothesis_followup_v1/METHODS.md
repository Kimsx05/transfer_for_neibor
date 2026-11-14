# 方法与冻结设计

这是已知发现后的定向内部检验，不是前瞻性独立验证。M0时间、原登记及哈希见 config/analysis_plan.yaml、M0_FROZEN.sha256 和 model_registry_v1_frozen.tsv。最终登记增加了预先允许的条件性临床敏感性36项：原Type字段40/40一致，3个dataset内含两类；T-stage缺14/40不插补。条件分支登记晚于主结果，完整记录于 clinical_addendum1.yaml，不按P选字段。原cluster字段实现修复见 technical_revision_cluster.md，原NE表保留。

## 输入与统计单位

最终c2l exact reference成员及原40个合格Tumor标本；116448个Epi。原MP10 Top20/分数不重算。RNA原始counts与原normalized RNA，非integrated/Harmony。原Assay5已核实90011 features ×384612 cells。Epi_all是已有真实全上皮原计数profile，并逐标本核验H50+L50互斥完整相加相等。状态从原冻结分位点与评分源核验；二次文本导出可改变极少数边界比较，因此不从舍入后的表重定义成员。所有细胞barcode/feature/profile唯一且对齐。

T4463209a为主，b仅替换；其余身份沿用derived/provisional，不把记录患者单位写为已核实独立患者。每profile≥50主、≥30敏感性；至少4记录单位、满秩、残差df≥3。不同High有重叠，自己的合格集合与共同集合分别解释。

## 分层与评分

无可追溯Phase，故在全部原Treg原RNA上用Seurat cc.genes.updated.2019 + CellCycleScoring，seed123、nbin24、ctrl为两组覆盖基因数较小者。背景是所有Treg检出基因；S/G2M operational cycling、G1 operational non-cycling，并非G0或实测增殖。分类无效者必须单列。本轮没有无效Phase。MKI67/TOP2A只用于检查，不调阈值。

仅新增必要Treg_NC/CY raw-count pseudobulk，复用原All-Treg和Epi_all。每个谱系表达宇宙固定自原≥30 profile的edgeR filterByExpr；Treg All/NC/CY共用All-Treg宇宙。新分数在保留完整RNA library size的TMM-logCPM(prior.count=2)上计算，固定每谱系main-a参考profile；ssGSEA使用GSVA 2.4.9、alpha=.25、normalize=FALSE。新评分至少5成员且该版本至少50%覆盖，未达者NE。

Epi M20从下游集合和rank/scoring universe移除Top20，不改counts/library/分数/成员。TF主为CollecTRI M20，敏感性DoRothEA A/B；ULM至少5targets，在main-a Epi_all参考逐基因均值中心化后计算。记录Full、unique-target和移除NECTIN2版本；NECTIN2不在当前SREBF靶集中，防循环版本不产生独立证据。旧v1活动使用其原filter/TMM/centering（过滤后library重算），完整保留原标签，新Epi_all定义不伪称数值等同旧活动。

MP10沿用原Q50/Q95固定SD；新X在main-a Epi_all≥50参考定标。Treg所有状态Y使用All-Treg main-a≥50的同一均值/SD。LODO、b、min30、Q95不重定标。不同receiver谱系基因背景不同，跨谱系只作描述比较，不作特异性结论。

## 推断与FDR

新增连续分数：Y~dataset+一个X，OLS HC3；双侧P与95%CI都使用同一HC3方差及残差df的t分布。h=1时HC3不可定义，保留NE而不删除singleton dataset。该限制不阻断满秩、df充分的edgeR计数模型。比例：cbind(cycling,noncycling)~dataset+X，quasi-binomial，实际模型Pearson离散度，t推断及OR/CI；无变异/退化/不收敛等NE。

主Epi10项、桥接18项、NC9项、次要Treg45项固定计划数BH；NE P/q为NA，n仍是计划数。共同样本All9作为单独对照家族。敏感性按登记设置单列家族，不能替换主q；同样β/P可能因家族计划数不同而q不同。模型表保存全输入样本、设计秩、df、覆盖、实际协变量、杠杆和Cook距离；不自动剔离群点。

旧配对DE沿用edgeR TMM robust QL和specimen+state。新增Treg计数分解使用同一registered样本/暴露，edgeR filterByExpr/TMM/robust QL，完整全基因P/FDR及主All/NC/CY联合基因q保留。GSEA排序sign(logFC)*sqrt(max(F,0))，固定原Hallmark/Reactome、minSize15/max500、scoreType std、eps1e-10。修复所有诊断异常模型的全部通路，nPermSimple10000→100000只按数值诊断升级。P在eps下界不是精确尾概率，单列标记；旧q不覆盖。新Treg分解Hallmark每模型完整族BH，主own All/NC/CY联合族另存。GSEA基因排序证据与样本级ssGSEA-HC3不同，不能混P/CI或当独立验证。

四代谢分支、共享上游及文献来源见signature_registry.tsv。小FA/异戊二烯化集合虽可ssGSEA评分，未满足GSEA minSize15时仍NE，不改变门槛。呼吸链结构亚基注释仅用于解释，不把多个基因当独立生物学重复；完整gene-level表保留。

## 资源与可复现性

所有R使用--vanilla，BLAS/OMP/MKL/data.table=1，CPU亲和性限制；本轮计算并发最高不超过4核，实际主要为1–2个单核进程。没有重跑Seurat整合/注释/CNV、CellChat、SCENIC或NicheNet，没有启动或干预HD。只保留小型pseudobulk、评分矩阵、表、标签与图。重画入口只读结果表与绘图配置，不读原始对象、不重算统计。

精确NC质量敏感性对照：匹配All模型使用NC的质量协变量值，另补9项OXPHOS和1项G2M线粒体比例控制；原主模型不变。见matched_control_addendum2.json。结构亚基同一UniProt别名只选择全Treg原counts最多的固定原始feature用于解释性汇总；完整DE不改动，见respiratory_complex_gene_registry.tsv。

历史冻结成员与分位数表报告数存在11个不同边界细胞的浮点精度差异，完整审计留档；未为匹配计数重分组，未改变任何30/50资格。配对结果按历史实际membership解释，不能将重读文本后的计数称为同一精确重建。Phase/G2M共享基因，不作为独立验证。
