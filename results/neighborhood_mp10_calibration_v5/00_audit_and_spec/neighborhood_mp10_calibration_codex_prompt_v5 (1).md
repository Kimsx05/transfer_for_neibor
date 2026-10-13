# Codex 执行任务：MP10 评分与上皮归属校准，以及固定 niche 的同样本直接比较（v5）

## 0. 用户决定、实际任务和停止范围

请完整阅读后执行实际分析，不要仅预检、复述方案或再次询问是否启动。

本轮承接 `neighborhood_mp10_source_v4`。用户已经看过结果，认可在选择机制研究方向前，先判断：

> niche1、niche6、niche7 等环境中的中央 MP10 差异，究竟能够代表上皮内在状态，还是可能由评分定义、中央混合组成或表达归属步骤造成？

### 0.1 必须完成的工作

1. 复用尽可能完整的上皮单细胞 counts 和原始 MP10 usage，按患者/来源核验冻结 Top20 的代表性；不能仅重复 v4 的 209 个 usage 匹配细胞。
2. 复用 v4 保存的 pseudo-spot **观测 counts 与真实上皮 counts**，实际运行上皮表达归属并比较估计值与真值，补上此前缺少的恢复验证。
3. 在有限、可复现的模拟中，分别检验固定上皮状态、改变中央非上皮组成，以及固定组成、改变上皮状态；重点识别 GC B、SPP1 TAM、Neutrophil 背景下的系统偏差、信号衰减或符号反转。
4. 冻结 niche0–niche8 后，实际完成同切片/同一生物分析单位的三个重点直接对比：
   - `niche1_minus_niche2`
   - `niche1_minus_niche6`
   - `niche7_minus_niche6`
5. 复核全部 Top20 的基因级表达、检测、排名及 AUC 贡献；保留全九群摘要和全部 dataset，不只汇报三个关注群或符合预期的来源。
6. 给出下一步优先研究对象的**证据分级**，而非强行给出一个阳性机制结论。

### 0.2 不变事项

- 正式 niche 固定为 `SAMPLE_R12_LEIDEN`，名称 `niche0`–`niche8`；不改群号、标签、群数、reference anchors 或接受规则。
- 环境继续使用已经冻结的 q05 组成；中央 gate 为 posterior mean Epi / all30 ≥0.25。
- 正式空间评分保留原 Top20、顺序、历史 AUCell-like 定义及已有 v4 数值，不替换基因、不因结果不理想扩大/缩小 signature。
- 原 NMF/固定基底 usage 仅用于评分代表性和状态校准，不自动替换正式空间 Top20 得分。
- 取消的“C3/C6 代表切片人工空间复核”仍不做。
- 病理资料、样本映射沿用 v4 已有内容；未知身份不阻塞本轮，不重启全量病理检索或要求用户先审完 133 ROI。
- 已经查看过空间 MP10，所有后续比较都是**结果知情后的探索/校准**；禁止称为前瞻验证、预注册、独立验证或从未解盲。

### 0.3 本轮允许与禁止

允许：读取现有模型和 counts、有限技术校准、已缓存表达按需重算、最多两个小规模**合成数据** cell2location 拟合任务，以及下述固定对比。

禁止：重新训练真实 63 张切片的 cell2location、重新训练全量 reference、重跑 CopyKAT、NMF、Harmony 或 niche 聚类；不启动 DestVI、CellChat/NicheNet、MISTyR、轨迹、生存、药物筛选或湿实验机制分析。

模拟用的辅助 reference 可以由训练供体的 counts 构建，其名称、用途与原模型严格分开；它不替换任何真实空间模型。

**完成本轮校准、直接比较、基因级核查和证据矩阵后停止。** 来源不确定、效应偏负或无效应均是允许的科学结果，不是需要“修成通过”的技术错误。

---

## 1. 工作目录、只读输入、版本与资源

```text
ROOT=/home/data/t070721/codex_workspace/Bladder Metabolism
OUT=<ROOT>/neighborhood_mp10_calibration_v5
V4=<ROOT>/neighborhood_mp10_source_v4
V3=<ROOT>/neighborhood_ring_harmony_v3
V2=<ROOT>/neighborhood_ring_comparison_v2
C2L=<ROOT>/cell2location_full_run/v2_30type
SCROOT=/home/data/t070721/Metabolistic_scRNA/Multi2
```

以上是定位线索，服务器实际 manifest 与 accepted 版本优先。目录已有同名任务时核对规格与 checkpoint；不得覆盖不同分析。

审阅来源：

```text
repository: Kimsx05/transfer_for_neibor
v4 reviewed commit: 49dee812902252b073d9ad1787491d5d79d7bc0d
v4 repository directory: results/neighborhood_mp10_source_v4
v3 repository directory: results/neighborhood_ring_harmony_v3_upload
v2 repository directory: results/neighborhood_ring_comparison_v2/neighborhood_ring_comparison_v2
```

旧文件、已训练模型、原始 counts、冻结分数、标签、患者映射、历史完成标记全部只读。新结果独立保存。不要直接调用会写回 v4 的 `main()` 或旧 `resume`；旧 `base.py` 把 OUT 写死为 v4，必须先消除写入旧目录的风险。

- 新代码统一 `seed=123`，派生种子用稳定字符串哈希；已有 v4 模拟保持其原始 seed 与 counts 不变，不伪称统一重生成。
- 先查配额、RAM、GPU、现有作业；不终止其他进程，不挤占整卡，不创建全体 spot×基因×30类型×抽样的巨型数组。
- 单细胞和空间表达使用稀疏矩阵/分块；只缓存必要的 score、Top20 计数/贡献、模拟真值和少量诊断。
- 不自动外发邮件、不上传私有患者信息、病理原图或完整 counts；轻量回传包必须去识别化。
- 新的所有函数都显式传 OUT，写入前检查位于本轮新目录。

---

## 2. M00：核对现有资产并冻结分析规格

### 2.1 优先读取的 v4 文件

```text
reports/RUN_SUMMARY.md
reports/METHODS.md
reports/WHAT_IS_STILL_UNCERTAIN.md
00_audit_and_freeze/resolved_inputs.tsv
00_audit_and_freeze/MP10_gene_availability_reaudit.tsv
03_epi_attribution/model_and_gene_universe_manifest.tsv
03_epi_attribution/attribution_technical_validation.tsv
04_niche_comparison/anchor_MP10_all_layers.parquet
04_niche_comparison/independent_unit_effects.tsv
04_niche_comparison/section_niche_score_summary.tsv
05_source_validation/source_validation_manifest.json
05_source_validation/single_cell_reference_scored.tsv.gz
05_source_validation/existing_EpiMP10_usage.tsv
05_source_validation/single_cell_usage_check.tsv
05_source_validation/central_composition_blueprints.tsv
05_source_validation/composition_only_pseudospots.tsv
05_source_validation/composition_donors_and_cellsize.tsv.gz
05_source_validation/state_change_controls.tsv
05_source_validation/state_control_paired_effects.tsv
05_source_validation/counts_cache/<section_id>/axis_manifest.json
05_source_validation/counts_cache/<section_id>/composition_observed.npz
05_source_validation/counts_cache/<section_id>/composition_true_epi.npz
05_source_validation/counts_cache/<section_id>/state_observed.npz
05_source_validation/counts_cache/<section_id>/state_true_epi.npz
code/02_scores.py
code/03_compare.py
code/04_source.py
```

`axis_manifest.json` 包含 genes、composition_ids、state_ids；必须按这些轴做一对一对齐，不能按合并后 DataFrame 的行号猜 counts 对应关系。

预期历史规模仅作审计线索，不作为必须凑齐的通过标准：

- 共同 anchors 129,632；接受标签 126,075；完整 Top20 配对主分析 125,855。
- 目标切片 63，有 anchors 的 61；D1/B13 无共同 anchors。
- GSE238145_B5 的 220 个 anchors 缺 HMGCS1，19 基因结果独立保留，不进入完整20基因的正式层间比较。
- v4 来源模拟为 9,600 单细胞、360 蓝本、1,800 composition-only bags、99 对状态蓝本；usage 实际匹配 209 个上皮细胞。

数量不符先查版本/UID，不删点或造点配平。

### 2.2 精确冻结 Top20

读取 v2 `00_preflight/MP10_top20.tsv` 和原 loading 来源。应核验的顺序为：

```text
HMGCS1, MSMO1, INSIG1, FDFT1, SQLE, FDPS, IDI1, ACSS2, ELOVL5,
LDLR, MVD, SCD, HMGCR, ACAT2, SLC20A1, ACSL3, ACSL1, ADIPOR2, MXD1, NPC1
```

不能以此文本强行覆盖源文件；不符要报告原因。保留 HMGCS1 缺失审计，不替补。基因符号、Ensembl 映射、重复行合并只复用可追溯规则，不能按预期效应任意改别名。

### 2.3 原始上皮对象和 usage 的定位

先从 manifest、脚本和冻结记录查找以下真实资产，不先遍历整个服务器：

- `07_EpiMP_basis_sample_balanced.csv`：原约 3,000 genes×14 MPs 的固定基底。
- `cell_by_14_EpiMP_usage_raw.csv`、`EpiMP_raw_usage_freeze_canonical_v1`。
- 原投影 `00_utils.R`、`08_project_to_cells.R`，以及 `06_final_program_to_EpiMP.csv`。
- 原完整上皮/上皮+cycling counts 对象（历史规模约123,699，只有线索，不是当前必须达到的数量）。
- 先核验历史 `PCA_ALLGENE1/epi_cycling_allgene_harmony_clustered_final_after_small_removal.qs` 是否为可用来源，再沿读取它的脚本找实际 counts/metadata。
- 历史 ST/mapping 下的 basis 和本地 Windows 传输副本只是线索；不假设服务器具有同样的 Windows 路径。

v4 `existing_EpiMP10_usage.tsv` 是可复用索引，但必须核对其定义和生成方式。**同名 K32、lower80、MP10Continuous、LLPS score 或其他状态分群不能冒充原 EpiMP10。**

区分 `original_stored_usage` 与 `fixed_basis_reprojection_usage`；后者并非完全独立生物真值，原程序基因包含 Top20 时存在预期相关性。

### 2.4 冻结分析单位和比较方向

- 复用 v4 的 accepted analysis-unit overlay：patient、section、sample、specimen 的概念分开。
- 已核实同一患者多切片不当作多个患者；确有来源支持的匿名独立单位可保留，身份依据分列。
- 仅有切片名但不清楚独立性的来源，单独描述；不伪造患者 UID。
- 重点对比固定为 `1−2`、`1−6`、`7−6`，正数总是前者较高，不因结果改变方向。
- 优先级 niche1/7/6 是结果知情后的资源分配，不是先验生物结论。

产物：`resolved_inputs.tsv`、`frozen_analysis_spec.yaml`、`input_checksums.tsv`、`identity_and_independence.tsv`、`execution_scope.md`。规格冻结后直接继续 M01–M06，不能只返回预检。

---

## 3. M01：完整上皮数据中的 Top20—原 MP10 一致性

### 3.1 实際扩大核验范围

读取可用的完整上皮 counts 与原 usage，按 canonical cell_id 一一对齐，保留 dataset、patient、reference_batch、normal/malignant 证据。不得用新 CNV 运行或样本 N/T 名称凭空补恶性标签。

原 usage 缺失时，仅对已定位且 preprocessing 可复现的细胞进行原固定14基底 NNLS 投影；不重新发现基底、不调整基因权重、不仅拟合 MP10 一列。不具备准确 preprocessing 时报告 `PROJECTION_NOT_COMPARABLE`，继续已匹配细胞及其他任务。

对完整可用上皮实际计算：

1. 同一既定基因背景下的原冻结 Top20 排名分数。
2. 原 usage 及其类型（raw/relative，使用哪一列，不能事后选相关性最高者）。
3. 细胞库大小、检测基因数、Top20 检测数和各基因 AUC 贡献。
4. 有明确依据时，正常和恶性上皮分别汇总；无依据就只叫 Epithelial。

### 3.2 一致性评估

- 输出总体 Spearman；更重要的是每患者/分析单位的 Spearman、覆盖数及散点/分箱趋势。
- 默认每单位至少30个匹配细胞且分数有足够变化才计算相关；不足不填零。
- 同时查看患者内 usage 排名高低和 Top20 排名的对应，避免仅由 dataset/患者间偏移制造总体相关。
- 记录相关与库大小/检测数的关系。分层或控制技术变量后的结果是辅助诊断，不宣称“真实状态已去噪”。
- 固定选取每个单位均衡细胞用于图；表格保留全量，不能只展示正相关患者。
- 补一个 donor-balanced 汇总；相关矩阵的单位是患者/供体，不把十万细胞当作独立患者。

**不将单个0.475、0.6等阈值作为通用有效/无效标准。** 报告数据能否支持较稳定的状态排序，以及适用的计数深度范围。

### 3.3 评分实现核查（有限技术对照）

历史实现为 `average` ties、K=ceil(0.05×G) 的 AUC-like 分数；不谎称与 R AUCell 的随机 tie 处理逐位相同。

- 在按来源/深度预先抽取的有限细胞与 pseudo-spots 上，对照官方 AUCell 的版本和 tie 策略；固定最多3个随机种子，使用相同 G、Top20 和 cutoff。
- 对 fractional Epi counts 与整数 observed counts 分开记录 ties、零值率、排序变化；不将小数 counts 取整。
- 官方对照只是实现敏感性，不因为其空间方向更“好”而替换主分数；不扫描多套 UCell/AUCell 参数挑赢家。
- 全零归属表达记 NA；有表达但 Top20 排名无贡献可以合法为0，二者必须分开。

产物：

```text
01_score_fidelity/cell_level_fidelity.parquet
01_score_fidelity/usage_provenance.tsv
01_score_fidelity/fidelity_by_patient_dataset.tsv
01_score_fidelity/coverage_and_missing_ids.tsv
01_score_fidelity/tie_and_implementation_audit.tsv
01_score_fidelity/score_fidelity_report.md
```

---

## 4. M02：Top20 基因级证据与统一敏感性

本模块覆盖**全部九群**，分别汇总原始观测表达和模型归属 Epi 表达；必须同切片、同模型基因背景。Epi 表达未缓存则分块复用已验证的 v4 分配函数，不重训模型。

### 4.1 基因级输出

对每个 Top20 基因、niche、section、dataset 保存：

- 原始 counts 检测率，原始 counts/库大小归一化表达摘要。
- 模型归属 counts 的均值、中位数、总量和相对于原始该基因 counts 的归属份额。
- 全基因排序背景中的 rank、是否进入 AUC cutoff、对总分的逐基因贡献。
- 三个直接对比的基因级方向和跨单位一致性。

`E_attributed > 0` 常继承原 spot 的检测，不叫“上皮实测检出率”。同样，稀少非零小数归属并不证明来源是上皮。原始分数和 Epi 分数使用的基因背景须可追溯，不能每个细胞独立删零基因后排序。

绘制两类图：按dataset的 Top20 基因×niche 图，以及三个重点对比的逐基因方向图。检测、表达及 AUC 贡献分开，不仅给一个总score。

### 4.2 Leave-one-gene-out（一次固定的20项诊断）

对每个基因依次移除一次，保留原全基因排序背景，明确重新归一化后的19基因诊断分数与原20基因分数的区别；同时保留固定原分母的贡献删除量。

- 对全部20个基因使用同一规则和三项同样本对比。
- 输出对比符号、效应大小和单基因贡献集中度。
- 不基于结果删除某基因并发布“新版更好signature”；不挑有利的 leave-one-out 作为主结果。
- 不能把 B5 的原有缺失误作一项随机删基因实验。

产物：

```text
02_gene_evidence/gene_by_section_niche.tsv.gz
02_gene_evidence/gene_direct_contrasts.tsv
02_gene_evidence/leave_one_gene_out_all20.tsv
02_gene_evidence/contribution_concentration.tsv
02_gene_evidence/gene_level_report.md
```

---

## 5. M03：实际执行上皮状态恢复校准（本轮核心）

### 5.1 先说清楚“真值”和“oracle”是什么

- `observed_synthetic_counts`：模拟中所有成分的 counts 之和。
- `true_epi_counts`：模拟生成时已知来源为上皮的那部分真实 counts。
- `true_epi_score`：在和估计结果一致的基因背景中，对 `true_epi_counts` 计算的分数。
- `estimated_epi_counts / score`：估计器仅从允许输入得到的结果，严禁读取 `true_epi_counts` 来构造分配权重。
- `true_cell_fraction`、`true_RNA_fraction` 与 c2l 的 signature weight 不同，不可直接互换。

真实上皮 counts 的总和或将其除以 observed 得到的逐基因份额，仅能作为真值/数学单元测试，**不能当作待验证模型的预测**。

实际条件分配公式沿用已核验的 v4 / 官方实现：

```text
mu[s,g] = m[g] * sum_f(w[s,f] * C[g,f]) + ambient[batch(s),g]
E_hat[s,g] = Y[s,g] * (m[g] * w[s,Epi] * C[g,Epi]) / mu[s,g]
```

这里 Y 必须是当前待验证的合成观测 counts；真实数据回查时才使用真实观测 counts。共同的 spot detection 因子在该比值中约去，不能重复乘入；ambient与reference的尺度必须一致。`mu=0,Y>0` 属于模型/轴/尺度异常，不能静默填0。`E_hat` 不取整，不以`Y×Epi_fraction`替代。

将“原模型权重回放”“给定已知混合组成但未知逐基因真值的恢复”“重新估计丰度后的恢复”分层报告，不能把第一种称为端到端验证。

### 5.2 M03A：先复用 v4 的全部缓存做固定权重回放

v4 缓存已保留 1,800 composition bags 和99对状态蓝本。对所有可用缓存：

1. 按 `axis_manifest.json` 校验矩阵行、基因列、bag_id/state_id；校验 counts 非负且 `E_true <= Y`。
2. 对应的真实 section/barcode 模型参数从原只读结果提取：C、w、m_g、ambient、gene顺序；gene映射必须一一对应。
3. 用该真实 anchor 的已存参数，把**合成观测 counts**代入与 v4 同样的条件分配函数，重新得到 `E_hat` 和 score。
4. 保存 `E_hat_score−E_true_score`、Top20计数偏差，以及高减低状态差的恢复情况。

其真实用途是检查“在既有权重下，改变观测 counts 会怎样影响归属”。**这里的 w 不是针对合成数据重新估计的，且旧模拟没有按真实模型加入平台效应/ambient，因此这只是 FROZEN_WEIGHT_REPLAY 压力测试。** 必须单列模型与模拟的失配，不能将偏差全部归咎于真实数据的 c2l 失败，也不能凭这一分支宣布恢复可靠。

不能为通过校验把合成Y中的真实E偷偷注入分子，或用Y自身构建真实参考签名。

### 5.3 M03B：构建有已知组成、可分解误差的控制面板

复用 M01 扩展后的上皮供体和现有30类型真实单细胞；不能全靠同一参考均值生成“太容易恢复”的伪数据。

#### A. 数学/模型匹配基线

用冻结模型的 reference 和已知生成权重、同一套平台系数、显式已知 ambient 生成一个小的正控制面板。进行分配并与期望及实现真值比较。随机生成的真实分子分配不要求逐基因完全恢复。

用途仅是核对参数尺度、维度、概率与数值实现，不能证明对真实供体异质性的泛化。

#### B. 真实供体混合面板：三种重要背景

至少包含 GC B-rich、SPP1 TAM-rich、Neutrophil-rich 三类中央组成背景，以及低目标类型对照。背景来自中央 c2l mean 组成分布，不是 R1/R2 的平均环境组成。

设计必须包括两种受控条件：

1. **状态固定、组成变化**：固定同一上皮供体集合、上皮RNA量、总UMI；改变非上皮部分中目标类型的比例，按比例替换其他非上皮类型，而非增加总细胞/RNA量。同一次重复中优先共用完全相同的 realized Epi count 向量，使真实Epi score严格相同；非上皮 counts 可随组成变化。另保留一套遵循实际中央 Epi 梯度的模拟，用于量化纯度/深度影响。不得把状态固定组内的纯抽样波动误报为归属偏差。
2. **组成固定、状态变化**：同一非上皮 realized counts、同一 Epi UMI、同一总UMI，仅将上皮换成低/高原 MP10 usage 供体。高低尽量在同患者内部定义，避免拿一个患者的高组对另一个患者的低组。

优先用原 usage 定义低/高状态，不能按本次待验证的 Top20 score 挑组后宣称独立状态恢复。某供体没有足够 usage 高低细胞，就记录不足；复用旧 AUC-defined 状态只能标成循环性技术对照，不能放入主要证据。

记录期望上皮状态与实际 `true_epi_score_high−low`。由于有限计数，期望高组的 realized Top20 分数未必更高，不能删掉这种样本；分别报告“原usage对Top20的保真度”和“估计器对已实现真值的保真度”。

#### C. 计数、基因与平台控制

- 保留 cell number/比例、各类型平均每细胞RNA量、RNA fraction、真实总UMI四种信息；严禁用细胞比例直接当转录本比例。
- 用真实单细胞平均 unnormalized profiles 混合后抽样，或一个明确记录的等价生成机制；评分前先混合 counts，不平均单细胞分数代替混合分数。
- gene-specific platform factor、ambient 是否加入、如何作用于各类型，全部显式记录；不能用无ambient模拟承受含大量ambient的估计却不说明。
- 用预定实际深度分位数（如P25/P50/P75）进行有限深度比较；上皮0.25/0.5/0.7等是实验条件，不新增真实数据筛选门槛。
- 不形成全参数笛卡尔积。初始可用360个既有蓝本，每蓝本5次技术重复；额外控制面板总量上限3,000 bags，先以来源覆盖、深度、纯度和目标类型范围选定。
- 保留供体、reference、基因背景、模拟种子和每种条件的原始清单，不挑选结果方向。

#### D. 参考泄漏控制

核对所有模拟供体是否参与过原 c2l reference 和 MP basis 的构建。

- 真正未参与过 reference 的供体：明确作为 held-out donor 技术验证。
- 找不到这种供体：原reference分支标 `IN_REFERENCE_OR_UNVERIFIED`，不得写“独立供体验证”。
- 可以从全量30类型训练供体构建一个辅助的 donor-disjoint reference，测试供体完全排除；按已核实患者划分，身份未知则按 reference_batch 划分并降低独立性表述。
- 辅助reference只用于模拟，不能替换真实63切片结果，也不能叫原模型的独立验证。
- c2l参考每类的尺度要和模拟的 per-cell profile/权重对应。给定真混合比例并不等于已知c2l权重；若使用类型尺度因子，应在训练供体或正控制面板中估计并冻结，不能用测试 `true_epi_counts` 优化。

#### E. 给定组成的恢复

采用冻结 reference 与已知生成权重/有依据的类型尺度，执行实际条件归属，允许 estimator 不知道上皮高低状态。

输出 `KNOWN_MIXTURE_FIXED_SIGNATURE`。这一层是在“丰度已知但细胞状态未知”条件下检查固定reference能否保留上皮状态，不称为完整 c2l end-to-end。

### 5.4 M03C：有限的合成数据丰度估计→归属全链校准

用户已同意补少量端到端校准；本轮允许实际拟合**合成数据**，不需要再次确认。

- 先固定一个具有完整20基因、30类型可映射的模型基因背景，构建不超过2,000 bags 的合成面板，覆盖上述三背景、深度/纯度和上皮状态。
- 为防某一基因背景代表不了其他数据，可以再做一个独立基因背景的小面板；**最多两个新合成拟合任务**。选择仅依赖基因完整度/来源覆盖，不能按 niche 关联方向挑选。
- 对新 synthetic AnnData 从头注册和运行 cell2location，使用已冻结的参考、对应生成条件和可追溯的训练参数；不要把真实切片的 w 当成新的后验结果。
- 不把 niche_id、真实 Epi state、E_true score 或 true_counts 输入模型，不将 high/low 分开定义为技术 batch。
- 默认 checkpoint 5,000 次优化迭代，尚未收敛才按预定步骤续至10,000，再到最多20,000；不无限重跑、不用真实上皮恢复成绩挑 epoch。应按训练目标变化和数值状态报告是否收敛。
- 只能写新目录 `03_attribution_calibration/synthetic_refits/`；不得保存回真实模型路径。
- 如原正式reference与供体尺度未匹配，先正控制查因；不把明显平台失配等同于上皮来源失败。
- 运行后从估计的参数分配Y，再对比真实E、true cell/RNA composition和上皮 score。报告丰度错误及状态恢复错误，区分它们。
- 对未收敛或不可恢复的运行标记 `TECHNICALLY_UNRESOLVED`；不给出“c2l无法恢复”的泛化结论。具体缺失依赖阻断时记录异常，仍完成M03A/B、真实对比和全部报告。

这是**小规模技术benchmark**，不是对真实患者上皮的外部验证。即使成功，也不能声称所有空间spots都已被净化。

### 5.5 所有校准分支统一输出什么？

至少包括：

- 真值与估计分数散点、Spearman、MAE、bias、估计/真值变化斜率；按供体、背景、深度和基因背景分层。
- 真实状态固定时，估计−真值误差是否随 SPP1 TAM、GC B、Neutrophil 组成单调变化。
- 成对状态 `Delta_true=score(E_high)−score(E_low)`、`Delta_hat=score(Ehat_high)−score(Ehat_low)`；保留零变化，不计算除以近零真差的“恢复比例”。
- 符号反转率、信号减弱幅度、估计误差的可重复性；变号需相对于realized truth，而不是仅相对于预期high标签。
- 除了严格符号，按重复抽样噪声与一个最低rank步长给出“接近零”的描述性标记，避免把微小正负抖动当系统反转。任何容差要预先写入规格，不以正向结果调阈值。
- Top20 基因级的 E_true/E_hat 份额，以及是否由固定Epi参考的少数基因主导。
- 背景标签是模拟条件，不是真实独立病例；五次重复不当成五个患者。

不要用一个相关阈值一票通过，也不要把总counts重构误差当来源恢复指标。若只完成回放，`end_to_end_status` 必须是未完成，不能借文件名“oracle”谎称完成。

产物：

```text
03_attribution_calibration/calibration_design.json
03_attribution_calibration/donor_reference_overlap.tsv
03_attribution_calibration/simulation_axis_and_scale_audit.tsv
03_attribution_calibration/frozen_replay_scores.parquet
03_attribution_calibration/known_mixture_scores.parquet
03_attribution_calibration/end_to_end_scores.parquet
03_attribution_calibration/paired_state_recovery.tsv
03_attribution_calibration/composition_dependent_bias.tsv
03_attribution_calibration/gene_source_recovery.tsv
03_attribution_calibration/calibration_summary_by_background.tsv
03_attribution_calibration/calibration_report.md
```

文件不存在时不造空完成记录；缺失分支在 report中说明 status和可执行的续跑命令。

---

## 6. M04：固定 niche 的真实数据同样本直接比较

本模块可在校准同时推进；结果最后一起解释。评分校准不支持来源时，真实数据结果仍报告，但其含义只保留为observed/model-attributed association。

### 6.1 输入和主要分数

主输入为 v4 `anchor_MP10_all_layers.parquet`，使用 accepted labels + FULL20 + 同中心完整配对的行。

主要两列：

```text
spot_MP10_model_universe
epi_MP10_model_universe
```

`spot_MP10_frozen` 作为历史衔接结果另列。不可跨基因背景直接相减；niche在同切片内的两层比较使用同一个基因列表。B5继续单列。

保留所有七个dataset：

```text
CanDis__UNVERIFIED
GSE171351
GSE238145
GSE285715
GSE319536
Inhousedata__UNVERIFIED
NPJ__UNVERIFIED
```

不因为某dataset反向而排除；来源命名未确认继续保留UNVERIFIED。

### 6.2 三个固定直接对比

对于每张切片 s 和一个方向固定的 pair A−B，计算：

```text
Delta_spot(s,A,B) = median(spot_score in A) - median(spot_score in B)
Delta_epi(s,A,B)  = median(epi_score in A)  - median(epi_score in B)
```

- 两群必须出现在同一张切片，且各有≥20个有效anchors；没有覆盖的对比记NA，不以0替代。
- 同一个患者有多张合格切片：先求各切片效应，再在患者内等权平均；其他已获支持的analysis units同理。
- 如果A只在患者的切片1、B只在其切片2，不能伪装成同切片直接比较；可以另列跨切片探索，但不并入主对比。
- 各dataset汇总独立单位效应：均值、中位数、正/负/近零数、范围和有效单位数。跨dataset给pooled-unit与dataset-equal摘要，说明二者估计目标不同；发现同一已核实患者跨dataset时，全局推断按真正analysis-unit归并/重采样，不能重复计作独立病例。
- 单位数量≥5时做固定2,000次analysis-unit bootstrap，区间标探索性、未包含模型不确定性；少于5只报点和方向，不给虚假的精确区间。
- 不用pooled逐spot Wilcoxon/t检验给出独立患者的结论。若报告任意p值，三pair×两主层共6项作为同一预先固定的探索族调整，并明确不等于独立验证；优先报告效应与区间。

所有pair汇总必须附单位清单和未达覆盖的原因。上轮两个“niche vs rest”效应之差不能替代这里的同样本直接对比。

### 6.3 纯度、深度与空间依赖检查

主分析先保留原样本内效应，不以全部30类型作为协变量把环境定义全部消掉。

同时记录两群中央 mean Epi fraction、model-gene UMI、nFeature、完整环比例、病理已知范围和空间坐标分布。

增加一套固定的同切片技术可比性敏感性：

- 根据 Epi fraction 与 log1p(model-gene UMI) 做匹配，不读MP10决定匹配。
- 固定上皮比例caliper=0.05、logUMI caliper=同切片两群合并SD的0.25；匹配不放回、固定种子，报告匹配算法与平衡。
- 两群匹配后各少于20时标不足，不能不断放宽caliper追求可估。
- 精确/适度匹配只针对中心技术/混合可比性，不声称排除了全部混杂；它改变的target population必须报告。
- 若同一群在空间上高度集中，按既有section-specific 8×8 lattice blocks汇总或进行固定空间块扰动作敏感性；跨block/有重叠rings不等于独立患者。
- Epi≥0.5、≥0.7和complete-rings只复用既有敏感性定义，并附各niche/单位留存率；0.7只剩少量点时不得拿它“挽救”某个方向。

作为统一参照保留全九群的原v4 versus-rest效应，明确其rest随切片变化；不重新对36个pair扫描后选择最好看的组合。

### 6.4 一个有限的连续组成诊断

除三个pair外，只对以下三个目标群增加原始环组成的连续检查：

```text
GC B
SPP1 TAM
Neutrophil
```

分别使用R1和独立R2，不将ring2误作累计两圈。中心MP10取中央spot，环内组成来自原始冻结q1/q2，不使用Harmony坐标当细胞丰度。

- 每切片记录目标类型是否有足够变异；常数/极稀少时不硬算。
- 报告连续原始关联，并做预定技术调整诊断：中央Epi fraction、logUMI、中央对应类型mean组成。
- 严重共线性时标为不可分辨，不同时把全部30个总和为1的比例加进模型。
- 这不是接触/距离排斥检验，负向只叫条件关联偏负，不叫“回避”。
- 按analysis unit与dataset汇总，不用空间点数获得伪精度。

产物：

```text
04_fixed_niche_contrasts/paired_section_effects.tsv
04_fixed_niche_contrasts/paired_unit_effects.tsv
04_fixed_niche_contrasts/paired_dataset_effects.tsv
04_fixed_niche_contrasts/pair_coverage_all_sections.tsv
04_fixed_niche_contrasts/matching_balance_and_effects.tsv
04_fixed_niche_contrasts/continuous_composition_diagnostics.tsv
04_fixed_niche_contrasts/all9_reference_summary.tsv
04_fixed_niche_contrasts/by_dataset/<dataset>/...
```

---

## 7. M05：联合证据分级，而非挑选“能讲通故事”的结果

对全部九群和三个重点pair建立一张可追溯的证据表，至少分列：

```text
niche_or_pair
score_fidelity_in_epi
observed_direction
model_epi_direction
paired_dataset_support
known_patient_vs_anonymous_unit_scope
gene_level_coherence
single_gene_dominance
fixed_weight_replay_result
known_mixture_recovery_result
synthetic_refit_result
donor_overlap_status
depth_purity_overlap
pathology_scope
major_failure_modes
interpretation_allowed
next_priority
reason
```

`interpretation_allowed` 可分为：

- `TECHNICAL_RECOVERY_SUPPORTED_WITHIN_TESTED_RANGE`：在所测深度/供体/背景中恢复较好，不等于纯上皮已被实测确认。
- `SOURCE_ALLOCATION_SENSITIVE`：固定真实状态时估计分数随组成系统改变，或状态变化衰减/反转。
- `SCORE_FIDELITY_LIMITED`：Top20无法稳定代表原usage，或少量基因/深度主导。
- `DATASET_HETEROGENEOUS`：技术恢复尚可，但真实效应明显依队列/组织条件变化。
- `INSUFFICIENT_COVERAGE` / `UNRESOLVED`。

这些是说明性分类，允许同时存在多个标记，不把某个任意阈值称作已验证金标准。niche效应CI跨零与来源恢复准确性是两个不同维度，不能用“效应显著”代替“归属有效”。

### 推荐时必须遵循的规则

1. 正向、负向和接近零都照实报告。niche6不需要被“纠正”为正。
2. 不以差预后故事选择评分或dataset；本轮不导入生存结局。
3. 不因为niche1与GC B有关就称TLS或判断预后；不因为niche7富集中性粒就称促肿瘤；不把niche6负向叫空间回避。
4. MP10主要是转录程序，不能自动等同胆固醇浓度、代谢通量或氧固醇产生；相关生物学可作为未来待检假说，不在本轮宣布。
5. 高AUC可以来自该类型自身表达、分配偏差或上皮状态等，须根据模拟与实际数据分别讨论。
6. 如果恢复不可靠：指出具体背景、基因、深度及reference限制，不广泛宣布模型无效，也不自动开始新模型来追求阳性。
7. 如果恢复较好但空间效应仍弱/异质：建议保留为探索或缩小可外推范围，不无限调参。

输出 `evidence_matrix_all9.tsv`、`priority_decision.md` 和每dataset的 `dataset_review.md`。最终优先级不必与讨论时的niche1→7→6相同，但任何变化需基于上述量化证据，不基于结果更符合预期。

---

## 8. M06：图表、审阅包和运行质量

### 8.1 图表

所有结果按dataset输出，少量切片dataset也保留。图表使用`niche0`–`niche8`中性名与固定顺序；图题标明observed、model-attributed Epi、true simulated Epi。

至少：

1. Top20 score vs 原usage：总体与按供体/患者的fidelity图。
2. 模拟真实Epi vs 估计Epi分数；分GC B、SPP1 TAM、Neutrophil背景显示。
3. 成对状态 Delta_true vs Delta_hat；零附近变化单独标记。
4. 状态固定的组成梯度→估计误差图；不同恢复层不混画成同一种方法。
5. 三个固定pair：每dataset的单位效应点图/森林图，两评分层分开对照。
6. Top20逐基因方向与AUC贡献图，以及全部20项leave-one-out敏感性图。
7. 全九群证据矩阵与覆盖热图，区分NA、0、未做、无足够独立单位。

透明PNG+矢量PDF，图例完整、轴单位明确、同panel色阶可比、无拥挤小字。引用已有空间图不新增人工空间复核任务。

### 8.2 验收与测试

- 对齐测试：spot/cell/bag/genes/30type矩阵轴，不能静默inner join丢行。
- 数值测试：非负、分母有效、Ehat<=Y、各类型+ambient重构等；只作为技术检查。
- 隔离测试：v2/v3/v4及真实模型哈希未被修改；合成模型写新路径。
- 冻结测试：niche标签、signature和主评分定义未改变。
- 模拟审计：真值没有进入估计器、供体泄漏被标记、非E realized counts固定条件逐行验证、总UMI/epiUMI匹配。
- 对比测试：A/B同切片真实共同覆盖，正负方向一致；技术重复不计为独立患者。
- 状态区分：脚本执行成功≠生物验证通过；资源/输入失败≠阴性生物结论。
- 所有NA有原因；缺失数据表保留预定schema，不能以空行/空文件假装完成。

### 8.3 输出结构与续跑

```text
neighborhood_mp10_calibration_v5/
  00_audit_and_spec/
  01_score_fidelity/
  02_gene_evidence/
  03_attribution_calibration/
  04_fixed_niche_contrasts/
  05_evidence_synthesis/
  06_figures/
  code/
  logs/
  reports/
    RUN_SUMMARY.md
    REVIEW_INDEX.md
    WHAT_IS_STILL_UNCERTAIN.md
    completion_matrix.tsv
    deviations.tsv
    review.zip
```

创建 `code/resume_mp10_calibration_v5.sh` 和阶段checkpoint，resume必须核验输入与规格hash再复用，不因旧COMPLETED存在就跳过新任务。结束时实际测试轻量检查入口；不能为测试resume又重跑完整合成训练。

`review.zip` 只包括报告、评分摘要、效应表、小图、必要代码和manifest；不含大counts、整套模型或私人身份。分dataset图包可单独生成，给准确路径与实际大小。

资源不足或依赖缺失时，完成其他模块并返回具体失败和续跑命令；不要把待做写成已完成，不承诺异步运行结束后再通知。

---

## 9. 最终回复必须回答的内容

请用中文返回：

1. **实际上做了哪些分支？** M03A/B/C分别运行多少bags/供体/背景/新模型，哪些失败或尚未运行；不以“oracle已做”笼统带过。
2. **Top20 能否稳定对应原MP10？** 总体与患者内是否不同，匹配细胞数是否实际扩展。
3. **是否出现组成导致的假状态差或符号反转？** 分GC B、SPP1 TAM、Neutrophil以及计数深度说明。
4. **niche1−2、niche1−6、niche7−6实际是什么结果？** 两评分层、每dataset效应、方向一致性、可用患者/匿名单位数，并附全部九群摘要。
5. **更值得往哪个niche推进？** 依据技术恢复和真实数据共同判断；不能以期望正向或不良预后来选。
6. **还缺什么？** 缺模型/metadata/原usage/供体真值/病理时说出具体文件和影响范围，不再次要求全量病理审定。
7. 核心文件绝对路径、轻量审阅包大小与可直接执行的续跑命令。

本轮不得报告接触回避、免疫抑制或上皮特异性已证明。正确的负向结果和未能确定来源的结果都需要完整交付。

---

## 10. 方法来源与表述边界

以下用于核对实现，不代表文献验证了本项目全部组合、阈值或九群：

- cell2location 官方条件性细胞类型表达分配实现：
  https://cell2location.readthedocs.io/en/latest/_modules/cell2location/models/_cell2location_module.html
- AUCell 官方 vignette（排名背景、cutoff、ties需要以安装版本为准）：
  https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
- v4 已审阅脚本和模拟资产：
  https://github.com/Kimsx05/transfer_for_neibor/tree/49dee812902252b073d9ad1787491d5d79d7bc0d/results/neighborhood_mp10_source_v4

主统计定义、min20/min30/min5、配对范围、模拟上限、匹配caliper和证据分级都是本轮工作规格，不是上述文献给出的通用金标准。
