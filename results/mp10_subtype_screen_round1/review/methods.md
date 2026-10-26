# MP10 第一轮方法与解释边界

本轮只执行共同样本表、比例、CellChat 与候选汇总。原始对象及既有结果均只读。

## 输入和统计单位

- 注释：`subtype/v2/final_annotation_pre_reference/frozen_annotation_barcode_mapping_v2.csv`，完整母对象为同目录 `final_annotation_v2_pre_reference.qs`。使用 `major_celltype`、`subtype_v2`。沿用已有污染剔除 `annotation_status_v2=excluded_before_round2`，不按 `spatial_reference_include` 过滤；原有 unresolved、dataset_bias、状态型标签全部保留。Cycling 按既有上皮合并政策纳入 Epithelial；其余七大群共同构成 TME。所有 subtype 的父群沿用映射，不把 Pericyte/SMC 重新划出 Fibroblast。
- 样本和临床：`cellchat/sample_by_clinical_group.csv`，标本为 `orig.ident`，UID 加 dataset 前缀。正常/肿瘤来自 `clinical_group`；T 分类仅解析明确 `fine_stage`，不把 NMIBC/MIBC 或 overall stage 改为 T1/T2。当前 T 记录没有 pT/cT 前缀，记为 T_unspecified_basis。未找到已确认的技术文库合并映射，不擅自合并标本。
- `patient_uid_recorded` 复用现有患者规则映射。derived/provisional 不能冒称经临床确认的患者；`patient_uid_verified` 保留 NA，所有“患者数”明确是映射 ID 数。配对和相关模型均以该映射为条件。Inhouse a/b 合并患者但保留两个标本，仍需核实。
- 恶性依据：`Multi2/3.3 recognize malignant/rough_cluster_infercnv_score_cell_metadata.csv` 的既有每细胞 rough_cluster，与 `Multi2/3.4 malignant cell process/NMF1/rough_cluster_malignant_like_status_summary.csv` 的冻结 `malignant_like_main` 映射。仅 Tumor 上皮进入目标；这是原有 malignant-like 判定，不是本轮独立验证的恶性真值。其余肿瘤上皮标为 unresolved，不强称良性。
- MP10：直接复用 `neighborhood_mp10_calibration_v5/01_score_fidelity/cell_level_fidelity.parquet:Top20_score`，基因与任务指定 20 基因集合逐一一致。该既有分数采用平均 ties 的 AUCell-like 排名统计，背景 90,011 features，K=ceil(0.05×features)，归一分母 sum(K−1,…,K−20)。不是本轮改算的官方随机 ties AUCell，不与 Top50 UCell 或 NMF usage 混称。
- 标本分数为上述目标细胞已有分数的算术均值。沿用项目 CellChat 目标细胞覆盖先例，≥30 个有分数的恶性样上皮细胞用于主要 MP10 关联；低覆盖仍保存均值、覆盖数并保留其他临床比较。各 dataset 合格 Tumor 标本中位数为阈值，严格 > 为 high，≤ 为 low。相同分数不随机拆分。

## 比例与模型

四种分母分别为全部有效细胞、全部非上皮 TME、每个 parent、全部 TME。父群存在而亚群缺失为 0；父群不存在为 NA。单一 subtype 的 parent 内恒定比例不检验。所有统计以标本为观测，细胞不作生物学重复。

使用 [speckle 官方 getTransformedProps](https://github.com/phipsonlab/speckle/blob/master/R/getTransformedProps.R) 原始函数（保存在 scripts；环境未安装完整 speckle 包），每个分母独立计算 `(count+0.5)/sum(count+0.5)` 后 natural logit。将转换矩阵传入 limma lmFit/eBayes，trend=FALSE，≥3 可估计群体时 robust=TRUE。参考 [官方教程](https://bioconductor.org/packages/release/bioc/vignettes/speckle/inst/doc/speckle.html)。组别模型调整 dataset；连续 MP10 使用原始分数，分期完整子集的未调整与分期调整使用同一批样本。每个 dataset 同时单独拟合，单队列自然不加 dataset 项。排名不足或残差自由度<2 不检验。

重复标本首先采用 duplicateCorrelation + patient block。仅一组肿瘤重复 ID 时部分子集无法单独估计相关；此时借用同一分母家族、全图谱标本 `~dataset+tissue_status` 拟合的患者内相关作为工作协方差。表中逐项标记 borrowed_full_atlas，并提供每个记录患者按 UID 字典序保留一个合格肿瘤标本的敏感性分析；借用相关依赖跨组织可迁移假设，不把暂定映射当已证实患者独立性。

BH 按比较问题×scope×注释层次×分母口径，在全部实际受检群体内统一校正；subtype/parent 跨所有父群校正。报告模型系数、moderated 95% CI、P/q、有效标本/映射患者/队列数，以及未调整原始比例差（百分点）。连续系数为每 1 个原始分数单位的 natural-logit 变化，不是原始比例倍数。T 总体分析为补充 omnibus F 检验，不作全部两两比较。

## 图 A

每队列柱高先对有效标本比例等权平均；总体柱仅纳入对应 parent 两组均有效的共同 dataset，再对 dataset 等权。各组使用相同 dataset 集合与权重。表保留有效标本/映射 ID 数、队列集合及权重规则。黑线连接左右柱相同累计分界，无 ribbon，不表示细胞迁移。subtype/TME 不重新归一化为 parent 的 100% 图。透明 PDF 和 300 dpi PNG。

## CellChat

使用现有 CellChat 1.6.1 human DB；输入为完整最终注释对象已有非负 RNA:data（LogNormalize），稀疏提取 1,015 个数据库基因。按 sample×communication_type 固定 seed=123 封顶 100 细胞，低于上限全部保留。记录原始与抽样细胞数。每条件网络采用 triMean、raw.use=TRUE、population.size=FALSE、nboot=100、seed=123、min.cells=30、内部 P<0.05。既有 MP10 Top5 细胞划分的旧网络不兼容本轮标本分组，故重新计算。本轮未获得已确认的分选信息，population.size=TRUE 不作为全队列主张或必要补充。

细分网络涵盖全部可评估亚群；大群图由这些网络聚合，标明是“细分 LR 有向边计数之和/强度之和”，不冒称大群独立推断。Normal_epithelial 与 Malignant_like_epithelial 区分；正常组缺少恶性群体为不可评估，不是已测零通讯。

比较使用条件发现 LR 并集。源/靶细胞不足、条件缺失、LR 未通过过表达预筛、计算后未检出分别标记；只有可评估群体的“当前标准下未检出”才在支持强度图中置零。原始 probability 缺失保持 NA。CellChat 内部置换 P 仅说明该合并条件网络，不是患者层面组间差异 P。CellChat 各对象内部按最大表达缩放，因此 pooled 概率差异只作探索，不能解释为患者层面校准效应。

全部标本×亚群的候选数据库基因平均 log-normalized RNA 和检出率先行缓存。shortlist 核查使用原始全细胞样本均值，不用抽样重复作重复，不检验独立 CellChat 对象的 probability。多亚基复合物逐基因保存支持，组合描述用平均表达几何均值与最小亚基检出率；双方每标本≥30 细胞才可评估。表达支持和网络支持单独陈述。

图 B 为真实单 LR 有向 circle：同一对照共用节点、顺序、坐标、颜色、等大节点与概率最大边宽标尺；无可评估群体显示不可评估，检出不足显示当前标准未检出。未使用 pathway 和代替单 LR，未预设 WNT5A–MCAM。绘图来源和共享标尺另存节点/边表。接口按安装版函数签名核对，参照 [CellChat computeCommunProb](https://github.com/jinworks/CellChat/blob/main/man/computeCommunProb.Rd) 和 [netVisual_circle](https://github.com/jinworks/CellChat/blob/main/man/netVisual_circle.Rd)。

计算加速仅将 CellChat `computeCommunProb()` 内的 matrix/group `aggregate(..., triMean)` 替换为 C++ 精确 type-7 四分位数实现，避免逐列 R data.frame 开销，未改变 LR、表达汇总定义、置换次数或随机种子。验证中分组汇总最大浮点误差 4.44e−16，小型真实输入网络 probability 和内部 P 最大差均为 0；见 `triMean_acceleration_validation.tsv`。安装的 CellChat 未被修改。其余步骤仍调用原包。初始未加速计算已在本任务内停止并以同参数加速版本继续，其他任务不受影响。
