# 分析方法与复现

## 项目、分组和证据单位

实际项目根为 `/home/data/t070721/codex_workspace/Bladder Metabolism`。指定的 `/home/jinshengxi/bladder_metabolism` 在本机不存在。用 rg 核查当前工程、父目录及用户目录，未定位到适用于此项目的 AGENTS.md，也未找到指定的 verified_resources JSON。使用本地项目工件及用户给定的 accession 和版本。结果目录在本轮首次执行时新建；中断后续写同一次运行。

从原 30 类注释的标本与条形码映射复用 Epi 身份、MP10_cell 和 pseudobulk。原分位数来自 `mp10_cellchat_high_definition_sensitivity_v1/scripts/01_prepare.R`，聚合来自 `mp10_immune_endothelial_mechanisms_v1/scripts/01_pseudobulk.R`。对每个标本的 Epi 用 type 7 分位数；H05 > Q95，H10 > Q90，H50 > Q50，L50 ≤ Q50。保持等于阈值的并列关系。仅纳入肿瘤来源 Epi，每组 ≥50 个细胞。按原记录保留 T4463209a，多个标本不重复贡献一个已确认患者。其余记录身份部分由数据集内命名派生，独立性尚未全部核实。主比较的有效数为 24 个标本及 24 个记录单位。原评分和 20 个成员保持原版本，见 `scripts/modules.tsv`。ACAT2 采用原基因，CTRP 的 Entrez ID 为 39；本轮 FSP1 指 AIFM2。

## 已完成、复用和本轮新增

| 已有工作 | 本轮复用 | 本轮新增 |
|---|---|---|
| immune_endothelial_mechanisms_v1，Epi 原始 pseudobulk、细胞检出、配对 edgeR、PROGENy | M0031 主 DE，M0016/M0001 敏感性，Epi 的有符号 ULM 足迹和原患者关系 | 三机制模块评分，Full/去重叠，基因检出、驱动、混杂、留一数据集及全 Epi 连续分析 |
| hypothesis_followup_v2，SREBF2 符号修正和 mTORC1 | 原结果显示 H05 的 SREBF2 去重叠差 +0.0276，mTORC1 去重叠 NES +2.744，两项同源于本项目 RNA | 用作脂质调控与一般生长竞争解释，核查可用网络的验证范围 |
| resistance_validation_v2_expanded，真实顺铂反应及耐药验证 | 对照已有实测 GDSC 顺铂，GDSC1 中 MP10–LN(IC50) ρ≈0.05，n=18；核对状态和报告 | 指定四药独立重取 CTRP 原配套表达与 AUC，另取 DepMap 同期表达和基线依赖 |
| dataset_distance_heterogeneity_v1，空间距离与生态关系 | 原空间单位和切片审计提示切片、spot 和患者计数有区别 | 作为组织 IF 的采样背景，未增加 Epi 统计单位或用空间距离代替膜定位 |
| 公开实验资料 | 作者的筛选分数、配套原计数、TPM、DE 与功能图源表 | 定向查询、真实分组及干预/补救对照，单列证据边界 |

## 本项目统计

冻结登记见 `PREREGISTRATION.md` 和 `scripts/modules.tsv`。机制模块为任务指定节点与候选酶构成的等权定向程序，生化参照包括 [GPX4/FSP1](https://www.nature.com/articles/s41586-019-1705-2)、[MUFA/ACSL3](https://pmc.ncbi.nlm.nih.gov/articles/PMC6430697/) 与 [CHMP4B 相关膜修复](https://www.nature.com/articles/s41467-026-74022-y)。易过氧化底物、甾醇、辅酶 Q 和蛋白异戊二烯化分列。每个模块按参与组的 TMM 表达矩阵计算平均 log2 CPM（edgeR prior.count=2）；精确符号匹配，计入至少有一个正计数的成员。逐模块输出覆盖和 High−Low 均值差的成员，另存原配对 edgeR 检验。家族平均可掩盖低丰度或相反方向成员，需结合其单基因检出。

主分析用每个标本配对差的均值、双侧 t 区间及 p。BH 家族为同一对比、版本、模型中的所有可估计模块，含对照与原分数读出。敏感性按预设 H10/H50、逐数据集剔除，以及回归配对差对增殖、应激与 log2 每细胞 counts 比值；调整效应为这些差异取零时的截距。补充梯度用全 Epi 模块对原 Q95 的 1 SD 增量回归，控制数据集。该量描述标本之间的梯度。细胞数量、局部区域、guide 或重抽样次数均未当成独立 n。

MAPK、PI3K 复用已存 Epi PROGENy 带符号 top500 先验及版本化表达背景上的 decoupleR ULM。连续梯度无原 Epi_all 足迹，使用原软件库和同类计算，在全 Epi 样本中筛选足够表达的基因、TMM 归一化、按基因居中，再拟合 ULM。保存实际覆盖。加权成员差 `mor × High−Low 表达差` 是定位驱动的辅助读出，不能直接加总还原 ULM t 值。另依原 30-cell 资格背景和原 TMM 流程重建每个 profile 的 ULM，分解经中心化表达与有符号权重的协方差，按该 profile 的相关系数转换因子缩放，再取标本配对差。实际成员贡献之和逐一还原已存足迹差，见 `tables/PROGENy_H05_exact_contributions.tsv`。RNA 跟踪转录足迹，膜胆固醇分布、膜力学和受体磷酸化单列待测。

## 配套药敏与依赖

CTRP v2.1 从官方归档读取七个指定文件。按 master_cpd_id 及精确名称双重核对。AUC 为拟合活力曲线的面积；各药 16 点，每步二倍。表达用配套 Affymetrix/RMA log2，以 gene_symbol=1 的合法 gene_primary_name 构成完整映射背景，同名探针先求均值，样本内全背景排名（并列平均），再平均 Top20 成员百分位。823 个有表达模型均覆盖 20/20，背景为 16,863 个基因。DepMap 用 24Q4 Public v1、ACH 行索引、精确 GPX4 (2879) 列。对 log2(TPM+1) 分块作全基因排名，合格模型中映射背景 19,193 个基因，Top20 全覆盖。各平台独立评分和统计。

组织学用 DepMap 24Q4 OncotreeCode=BLCA，CTRP 泌尿道模型通过唯一标准名精确匹配。身份注释补全后，重算所有交集与关联；四药始终使用 CTRP 的原表达。模型内重复 AUC 取平均。各药仅在同时有表达和实测终点的模型上，按同一患者同谱系取均值。患者映射缺失的泛癌模型保留模型 ID。BLCA 双侧 Spearman p，2,000 次按记录单位 bootstrap 的百分位 95% CI。逐一留出和最大秩残差用于标注影响，无事后剔除。泛癌保留 n≥3 且两变量有变异的癌种，在癌种内转成百分位秩并居中，合并相关；p 为消耗癌种自由度的近似检验，CI 在癌种内分层重抽 1,000 次并重新排名。癌种内原始 Spearman 另存。总体 BH 按平台与人群；单癌种表在各终点内校正。多药来源重叠，临床治疗效果未在此测量。

## 真实扰动

- **GSE182638**。18 列逐一匹配 GSM、真实列名和 passage。每条细胞系内 TMM log2 CPM，两个方向均用 3 个 passage 配对；记录是否更接近 UT、以及 3 个批次中有几个更近。成员 DE 用 passage+处理的 edgeR QL。模块 BH 在细胞系、版本、对比内。RSL3 5 µM 3 h；Fer-1 2 µM 预处理 2 h。
- **GSE126696**。官方 CP0041 实际取回，核对 construct ID、barcode、序列和 Addgene Brunello 的唯一正向映射，同时统计反向互补匹配。以整库 RPM 归一化，比较 log2(RPM+0.5)，gene 用 guide 中位数；全基因百分位用于识别全局耗竭背景中的相对富集。原研究用药 5 µM，个别 GEO 元字段写 5 mM。每个时点单列，只给描述效应；DMSO 只有一列。参见 [原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC7353921/)。
- **GSE175862/175863**。前者 8 条 GSM 而作者 CSV 只有 gene 与 Screen score 两列，保留排序、覆盖和原分数。后者 Excel 的 TPM Abundance 有 12 个样本列，按 GSM 精确对应，分 3 个模型/干预类型做扰动减对照。两组各 n=2，以 log2(TPM+1) 差与 Welch 区间描述，样本量导致推断精度有限。DESeq2 sheet 的作者合并设计另存；各独立对比及掩码内 BH。参见 [NK 原研究](https://www.nature.com/articles/s41467-022-29469-0)。
- **GSE149933**。用作者单独发布的两套 library/readCounts、GEO family SOFT 的 175 条 GSM 及固定提交的 mValidationAnalysis.Rmd。匹配 guide ID、序列、基因，并以 babelgene 的 ≥3 来源支持、全表双向唯一 Entrez 映射核一对一同源。基因 symbol 无合格映射或无文库覆盖时记 NE。mTKO 额外 18 列未匹配该系列 GSM，列名完整保留。先按模型（Renca 和 Renca-HA 分开）、文库、处理确定中段/终段，再按该阶段的真实 control/treatment 日数对照，明确异日对比。每个文库按样本全库 RPM，取 log2(RPM+0.5)，每个重复取基因的 guide 中位数，再估计三对三 Welch 差/区间，按文库、模型、阶段及处理的覆盖候选作 BH。未处理组相对 T0 的变化作基线生长参照。该候选估计保留原研究设计，效应量和 q 属本轮计算，与原文的 DrugZ NormZ 分列理解。参见 [原研究与代码](https://github.com/NMikolajewicz/Lawson2020)。
- **HMGCS1 功能源表**。索引 96 个 sheet，提取 Fig. 1c、2h、3a、3b、4e，保留行列坐标和原值。每条件三个源值，按模型、处理、对照和补救计算差；Welch 区间及逐 sheet BH 属本轮描述性重算。RNA KD 的参考真值、膜修复速率与患者疗效另列。GSE262114 按指定 HOLD，GSE211072 为可选 SOX13 敲低补充，本轮未用。

## 运行与输入完整性

从当前工程执行 `bash results/mp10_epithelial_function_parallel_v1/scripts/run_analysis.sh`。使用已校验缓存，自动建立带日期的新子目录，单 CPU 串行；每个阶段的共享系统锁覆盖整个进程生命周期。没有新增共用依赖。R 使用项目已有 decoupleR 库；版本见 `logs/R_session.txt`、`logs/python_versions.tsv`。下载代码及重试日志保留，官方矩阵按版本、原始长度和 MD5 校验。磁盘配额不足后将本轮两个 DepMap CSV 转为无损 XZ；对解压字节重新核查官方 MD5。最终资源台账为 `tables/resources_final.tsv`；`tables/resources.tsv` 是首次下载历史。补充总表使用 SOFT 取得完整 GSM 字段。项目和实验统计均未使用预测 IC50、新 CNV 或 FASTQ。

