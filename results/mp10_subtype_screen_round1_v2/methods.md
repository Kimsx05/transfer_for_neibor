# Round1 v2 方法与解释边界

## 输入、单位与保留原则

使用 `subtype/v2/final_annotation_pre_reference` 的冻结注释和 RNA 数据。374,329 个有效细胞，Cycling 按已有规则并入 Epithelial；只沿用 `excluded_before_round2` 排除，不按空间参考 include 标志删细胞。上皮纳入不受旧 malignant-like、rough cluster、CNV 或旧 score barcode 列表限制。Plasma 保持独立 major。所有未解析、队列偏倚等既有标签保留，不重注释。

标本主键是 dataset + orig.ident。临床映射来自 `cellchat/sample_by_clinical_group.csv`；保留已记录患者 ID 的 derived/provisional 来源，不称为核实患者。58 个标本对应 42 个记录 ID；分开记录标本数与 ID 数。T 只用原始明确的 Ta/T1/T2/T3/T4 等字段；当前未注明 pT/cT 的记录标为 T_unspecified_basis，不从 MIBC/NMIBC 或整体 stage 推断 T。Normal 来源、治疗与富集方式未记录；不臆测。

## 全上皮 MP10

固定 20 基因是 HMGCS1、MSMO1、INSIG1、FDFT1、SQLE、FDPS、IDI1、ACSS2、ELOVL5、LDLR、MVD、SCD、HMGCR、ACAT2、SLC20A1、ACSL3、ACSL1、ADIPOR2、MXD1、NPC1。

123,699 个已有细胞按 barcode 复用 `cell_level_fidelity.parquet:Top20_score`；18,182 个新增上皮从最终母对象 RNA counts 补算。冻结背景取原 `official_genes.txt` 并与 counts_cache/genes.txt 逐项核对，共 90,011 features。每细胞按 counts 降序、并列平均秩；K=4,501，分数为 Σmax(0,K−rank_gene)/89,810。零 counts 为 NA；不在信号基因、HVG 或 20 基因子集内重排。抽查五队列共 60 个重叠细胞，检查原 counts 和分数的一致性，见 `score_validation.tsv`；抽查覆盖不同分数及深度。score_source 分 reused / filled_same_definition / unavailable。

标本 MP10 是全部上皮有限分数的算术均值，另给中位数。肿瘤标本须完整可信覆盖且 ≥30 个有分数 epi，才参加标本 MP10 模型。各 dataset 在这批标本上固定均值的中位数：high 严格大于阈值，low 小于等于阈值。后续不随分期、父群覆盖、患者去重重新分组。

标本内 epi 状态另行定义：每个肿瘤标本按全部有效 epi 分数的中位数分 high/low，严格大于为 high，其余为 low；不拆平分 ties。只有该标本两类 epi 均 ≥30 时进入主状态网络。TME 不赋 epi 状态、不复制。两套分组回答不同问题，字段与文件名称分别保存。

## 比例分析

分母包括 major/all_cells、major/TME、subtype/parent、subtype/TME。计数全零补全；分母>0 的零比例保留，分母=0 为 NA。Plasma 单标签的 parent 内比例恒为1，不检验。父群 B<30 不剔除 B/TME、B/all_cells、B-subtype/TME 主分析。

唯一主转换是 asin(sqrt(raw proportion))，不加伪计数。用 limma lmFit/eBayes，trend=FALSE，至少3个可变类别时 robust=TRUE。全部结果保留效应、95% CI、P、q 与不可估计原因；效应是 asin-square-root 差或每1原始 MP10 单位的变化，不叫 logFC。多 T 的 F 检验无方向。Normal、T1、sample low 为参照。ALL 模型调整 dataset；每个队列单独分析。

比较包括 Normal/Tumor、保留配对记录 ID 的 Normal/Tumor 补充、T1/T2、可靠 T 总体检验、sample MP10 二分、连续 MP10，以及同一分期完整子集的连续 MP10 调整前/后模型。二分类每组至少两个不同记录患者 ID、设计满秩、残差自由度≥2。正常/肿瘤配对补充以记录 ID 固定效应保留配对；记录ID嵌套于dataset，故不再加入与患者固定效应冗余的dataset列。其他重复标本优先在当前分析内估计 duplicateCorrelation；仅一个重复 ID 的小块或无法估计时，按 sample_uid 固定保留每个 ID 首个合格标本，标为 one_specimen_fallback。不得借用全 atlas 相关系数。另输出 MP10 每记录 ID 一标本的敏感性。

subtype/parent 的三组预设分析：primary_asin (parent>0)、parent30_asin (parent≥30)、parent30_logit。logit 严格采用官方公式：在该父群全部类别计数上加0.5、归一化后 log(p/(1−p))。两种 parent30 转换使用相同实际拟合标本和协变量；分期调整前后也使用同一实际集合。若其中任一模型需要患者去重，对配对比较统一采用相同固定子集。样本映射见 `model_sample_membership.tsv`。

BH family 固定为 analysis_variant × comparison × scope × annotation_level × denominator_type；subtype/parent 跨全部父谱系受检亚群一起校正，不单独对 B 缩小 family。

图 A 画原比例，队列内标本等权；总体只取该比较、该父群两组均有有效标本的共同 dataset，然后队列等权。黑线连接累计堆叠分界，只是视觉连接，不代表迁移。原始 B 点图保留低 B 标本，B<30 空心标记。marker 复核采用项目冻结 B 身份/Naive/Memory/GC/IFN 面板，不重聚类或改标签；同时标细胞和标本覆盖。

## CellChat 主状态与次级条件

主分析每 dataset 建一个共享对象：同一批合格肿瘤标本的 Epi_MP10_low、Epi_MP10_high 和一次 TME；另建 ALL pooled 对象作为概览。每标本两类 epi 分别等量抽取 min(100,low数,high数)；TME 按标本×节点各封顶100，稀少群体全留。major 与 fine 各自从原始细胞独立抽样和推断，不由 fine 汇总替代 major。

次级 Normal/Tumor、T1/T2、sample_MP10_group_allEpi 各自确定标本集合，使用全部上皮的 Epithelial 节点；各条件分别推断，不强行限制为主状态共同标本。五队列和 ALL、major 和 fine 共68个主/次级网络。另一次主要 B/Plasma 共同覆盖敏感性：仅保留两类 epi、B、Plasma 各≥30的标本，在 major 层重新推断五队列及 ALL 六个主状态对象，不逐 LR 重跑或调阈值。

RNA:data 使用母对象已归一化非负稀疏表达；直接只读引用 v1 的 1,015 个数据库信号基因缓存，其 barcode 覆盖整个原始注释对象，逐 barcode 核对。CellChat 1.6.1 human DB，triMean、raw.use=TRUE、population.size=FALSE、nboot=100、seed.use=123、min.cells=30；保持其表达筛选和置换路径。精确 C++ triMean 加速只替代聚合循环，复用此前对原函数的逐值验证（概率和内部 P 差均为0）；编译缓存独立放 v2。没有重用 v1 已推断网络。

每个对象的完整数值数组保存在独立小型 qs，非零原始 LR 边导出完整表；raw_probability 不因内部 P 不支持而置零。主要 LR 对照只在两端节点及基因可测、同 LR 已计算时求 raw_high−raw_low。未计算、基因无法测量/匹配、节点不足与已计算但无支持分别标记。RNA feature 在母对象存在不保证某队列原始测量清单完备：额外报告逐队列基因非零证据；全队列全零且无法核实测量时记 measurement_unresolved，不能直接解释成测得的生物学零。不同节点内的真实零在有队列测量证据时保留。

仅展示用的 supported_display_probability 在内部 P<0.05 且概率>0时取原概率，否则已计算值显示0；未计算仍为 NA。它不是患者差异 P。次级不同对象各自按最大表达归一化，因此跨对象、跨 dataset 绝对概率和边宽不作可比患者效应。ALL 不是独立重复。

## 标本/记录患者表达复核与候选选择

未抽样的完整表达按标本×epi状态和标本×TME节点缓存，major/fine 分别直接聚合原细胞。候选必须同标本 low epi、high epi、目标 TME 各≥30。epi→TME 比较 epi ligand，目标 receptor 共享；反方向比较 epi receptor，TME ligand 共享。复合体逐组分保存均值及检测率，辅助汇总用几何均值、最小组分检测率。每状态/另一端独立报告≥10%检测支持；允许一类 epi 的实际零值，不要求两状态同时通过。乘积仅为表达辅助量，不是新 CellChat probability。

重复记录 ID 的标本差先平均为患者差，再报告正/负/零计数和 leave-one-patient 方向范围；不制造差异 P。候选按同向队列覆盖、方向稳定队列数、可用记录患者覆盖和配对表达差排序，不做任意加权总分，不按字母序打破无信息并列。至少两个 dataset、各至少两个可评估 ID 的同向证据只代表更完整覆盖，不是显著。没有任何队列达到两个 ID 的 GC B 明确不能排序；保留低覆盖结果。

圈图直接从真实 LR ID 原始概率提取。两个状态面板始终保留两个实际 epi 节点，分别只画对应状态的边；另给同圈并列图。颜色/位置/等大节点固定，两面板共用边宽最大值。所有边保存实际 source、target、raw probability、内部 P、status、展示值及标尺。包括反向、无支持、反向队列和 NE 图。

## 空间承接及来源

本轮只读取既有 `neighborhood_niche_celltype_screen_v6/00_audit/niche_label_map.tsv` 和轻量 unit composition summary，未重新计算空间数据。其 accepted SAMPLE_R12_LEIDEN 编号和粗粒度标签与本轮细胞亚群的对应仍需下一轮显式核对。B-rich niche 的 B 富集是定义的一部分，不是独立验证；共享 TME 网络不证明真实接触或激活/抑制。

实现口径参照本地历史评分源码与冻结 counts_cache、保存的 speckle 原函数、安装版 CellChat API。可追溯资料：[speckle 说明](https://bioconductor.org/packages/release/bioc/vignettes/speckle/inst/doc/speckle.html)、[原转换函数](https://github.com/phipsonlab/speckle/blob/master/R/getTransformedProps.R)、[CellChat 原方法](https://www.nature.com/articles/s41467-021-21246-9)。本项目状态分组、共享 TME 和表达复核是本任务预设设计，不声称为某篇文献已完整验证的同款流程。
