# Codex 执行任务：先检索病理证据，再完成 niche0–niche8 的中央 spot / 上皮归属 MP10 分析

## 0. 本轮授权、目标与停止范围

请完整阅读后直接执行实际工作，不要只完成预检、列计划或重新询问是否启动。

用户已确认：

1. 采用已有 `SAMPLE_R12_LEIDEN` 工作分类：sample 层面弱 Harmony、R1/R2 拼接、Leiden 九群。
2. 展示名统一为 `niche0`–`niche8`，与原标签 `0`–`8` 一一对应；不另造 CAF-niche、SPP1-niche 等正式名称。原始组成描述只作为注释。
3. 取消“C3/C6 代表切片的额外空间复核”任务，不把它作为任何阶段的验收前提。
4. **现在授权读取 MP10 数值并进行关联分析与上皮来源评估。** 旧 v3 的“止于 niche / MP10 deferred”仅描述历史运行，不再限制本轮。不要修改历史冻结文件来伪装旧运行已解盲；在新目录中单独记录本次授权。
5. **先查找可复用的病理标注和样本证据，再启动评分；病理暂时找不到或仍需人工判断，不阻塞评分、图表及上皮表达归属工作。**
6. 本轮同时交付：中央 spot 的整体 MP10、模型估计的上皮归属 MP10、九个 niche 的分布及样本内比较，以及当前能够完成的来源评估。
7. 不进行新一轮 niche 聚类、Harmony 选参、CellChat/NicheNet、MISTy、轨迹、预后或实验机制分析。不启动全量 DestVI、重新训练 cell2location、重跑 CopyKAT。

这里用户所说的“纯净 epi 打分”，在常规 Visium 中操作性定义为：**基于已有 cell2location 模型，从中央 spot 的观测 counts 中条件性分配给上皮的表达，再计算 MP10。**

它不是分选/显微切割后实测的纯上皮表达，也不是已证实的恶性上皮专属状态。正式字段和图注使用 `model-attributed epithelial MP10` / “模型估计上皮归属 MP10”。必须实际尝试获得这一层分数，不能只以高 Epi 阈值代替。

---

## 1. 工作目录、只读来源和资源要求

### 1.1 目录

主项目根目录：

```text
/home/data/t070721/codex_workspace/Bladder Metabolism
```

本轮新输出目录：

```text
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_mp10_source_v4
```

优先复用：

```text
<ROOT>/neighborhood_ring_harmony_v3
<ROOT>/neighborhood_ring_comparison_v2
<ROOT>/cell2location_full_run/v2_30type
<ROOT>/neighborhood_v1/00_inputs_and_audit/metadata/frozen_MP_scores/by_section
<ROOT>/spatial_ecology_top20_v2
/home/data/t070721/Metabolistic_scRNA/Multi2
```

以上是已知定位线索，不代表每个路径或子路径均必然存在。先查现有 manifest、README、输入表和脚本，定位真实文件；输出 `resolved_inputs.tsv`。目录存在替代位置时记录对应关系，不重复复制大对象。

GitHub 审阅版本：

```text
repository: Kimsx05/transfer_for_neibor
reviewed commit: a785e2710565ac31c586c5e1d09cfabfd7277b5b
niche results: results/neighborhood_ring_harmony_v3_upload
legacy score audit: results/neighborhood_ring_comparison_v2/neighborhood_ring_comparison_v2
```

服务器 accepted 输入优先于下载副本，但要核对版本与哈希。不要把目录名称类似的 K32、NNLS、LLPS 或其他 MP10 试验替换成本轮冻结评分定义。

### 1.2 不得改动的对象

旧版 models、cell2location reference、q05/mean 结果、niche 标签、anchor gate、环定义、MP10 Top20 signature、冻结原始分数、患者映射证据及旧 `COMPLETED` 文件全部只读。

不要原样运行旧的 `resume.sh`、`05_mp10.py` 或 `resume_niche_only.sh`。可以审阅和复用其中验证过的函数，但新脚本必须只处理本轮九群，并移除“所有 ROI 审定后才能计算”的全局门槛。

### 1.3 运行约束

- `seed=123`；派生随机种子使用稳定的 section/spot 哈希，记录算法，不使用进程随机化的 Python `hash()`。
- 先检查磁盘配额、可用 RAM、CPU/GPU 和现有作业；不终止、不挤占其他任务。
- 默认逐切片执行，表达分解按小批次处理，CPU 线程保守设置；根据实测资源调整。
- 不生成全体 spots × 全基因 × 30 cell types × 多 posterior draws 的巨型张量。
- 可检索并下载公开资料；不自动发送邮件联系作者，不向公共仓库或外部服务上传未公开病理图像、真实患者标识或原始私有数据。
- 保存主要是 score 表、摘要、必要的 Epi 稀疏表达缓存；不复制全部原始 counts，不自动上传原始病理大图或患者识别资料。
- 所有阶段都需 checkpoint、日志、失败原因与可续跑命令。单切片失败不使其余数据全部停止。

---

## 2. 先冻结本轮使用的 niche 和分析单位

### 2.1 固定九群

读取 v3 的：

```text
frozen_niche_manifest.json
03_main/SAMPLE_R12_LEIDEN/working_config.json
03_main/SAMPLE_R12_LEIDEN/labels.parquet
```

若 `03_main` 是符号链接，解析至真实位置；GitHub 中对应结果可位于：

```text
02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/
```

预期配置：

```yaml
config_id: SAMPLE_R12_LEIDEN
representation: SAMPLE_WEAK_R12
harmony_theta: 1
harmony_lambda: 5
n_pc: 22
leiden_resolution: 0.1
n_clusters: 9
n_common_anchors_expected: 129632
n_complete_rings_expected: 108345
n_reference_expected: 20000
n_target_sections_expected: 63
n_sections_with_anchors_expected: 61
```

这些数字用于版本核验，不得通过删点、补点或更改标签强行配平。若不一致，先查 UID 和版本，并清楚报告。

显示映射固定：`0 -> niche0`，……，`8 -> niche8`，不要改成 1–9。

- 主分析：`accepted_label >= 0`。
- `candidate_label` 仅作明确标记的覆盖敏感性结果。
- `accepted_label == -1` 保留 `UNASSIGNED`，不偷偷分配进九群。
- 保留所有共同 anchor 行；缺分数、缺模型、缺病理分别记 NA 和原因，不整体删除记录。
- 本轮对全部 accepted anchors 评分，不限于 20,000 个聚类 reference。
- 原 gate 继续是 **中央 posterior mean Epi / all30 >= 0.25**；环境继续是既有 q05 R1/R2 定义。这不要求上皮表达分解也使用 q05。

生成 `niche_name_map.tsv`、`niche_input_hashes.tsv` 和新的 `analysis_authorization.json`。解盲后不再根据 MP10 结果调整聚类、分辨率、群数、合并规则或 signature。

### 2.2 样本身份

保留 dataset、section、sample、specimen、patient、analysis_unit 的不同含义：

- v3 `sample_uid = dataset_id::section_id` 是空间测量单元，不自动等于患者。
- 已核实同一患者的多张切片应归并处理。
- 来源设计支持独立性的匿名样本可以有证据注明的 analysis unit；不伪造 patient UID。
- 无法判断独立性的切片可用于描述，不冒充独立患者进入患者级推断。
- N/T 配对必须有记录支持，不能只看命名。

---

## M00. 定位全部现有资产并登记

先读取 v2/v3 的审计和输入 manifest，至少定位：

1. 固定九群的标签、共同 anchor 表和空间坐标。
2. 冻结 MP10 Top20 的精确 20 个基因、顺序、来源、hash。
3. `MP10_score_assets.tsv`、`MP10_gene_availability.tsv` 和原始 AUCell 生成脚本。
4. 各切片原始 counts、gene IDs、表达层、保存的过滤基因清单及 normalization/ranking 设置。
5. accepted `v2_30type` 对应的训练模型、reference signatures、训练 h5ad/AnnData registry、posterior 参数及已有 Epi expression layer。
6. H&E、scalefactors、tissue positions、已存在的病理 mask/spot annotations/polygons 及原始样本映射。
7. 可复用单细胞 counts、30 类映射、恶性/非恶性上皮证据、MP10 原 NMF loading/score 及已完成的来源模拟。

先看这些历史来源验证文件（路径以实际定位为准）：

```text
neighborhood_ring_comparison_v2/07_source_validation/single_cell_audit.json
neighborhood_ring_comparison_v2/07_source_validation/general_pseudospot.tsv
neighborhood_ring_comparison_v2/07_source_validation/prospective_association_validation.json
```

旧的 9,600 单细胞/260 pseudo-spots 只能在确认文件、基因范围、来源和方法一致后复用；不能把“已跑过一般模拟”自动判定为当前九群来源验证完成。

`resolved_inputs.tsv` 至少包含角色、绝对路径、来源版本、文件大小、hash/已有hash、是否可读、UID/基因对齐状态、可以复用的内容和缺失项。

完成资产定位后立即推进 M01 以及可执行的评分，不把预检本身作为最终交付。

---

## M01. 第一项实际工作：查找已有病理/组织身份标注

### M01A. 本地优先，再查原始公共来源

覆盖 7 个 datasets，并为 63 张目标切片均保留状态。先检索项目已有 metadata、样本导入记录、原论文/补充材料、下载清单和代码仓库，再按确证的 accession 查询原始公共资料。

检索重点：

- 原作者的 tumor / normal / stroma 等 spot-level 标签；
- H&E 上的 ROI polygons、分区 mask、QuPath/Loupe 导出、人工病理标签；
- 与确切 section 对应的病理说明、治疗/取材信息及样本配对；
- GEO/CNSA/GSA 等补充文件、作者 GitHub/Zenodo 和原论文补充表。

优先处理现有五个大 ROI 所属切片：

```text
CanDis_2026_Sample8
GSE319536_BC11
CanDis_2026_Sample1
CanDis_2026_Sample3
GSE285715_T1
```

同时检查 `tumor_supported` 与已知非肿瘤来源发生冲突的区域。`CanDis__UNVERIFIED`、`NPJ__UNVERIFIED`、`Inhousedata__UNVERIFIED` 必须先追踪本地导入映射，不因搜索到相似标题就自行认定来源。

每个 dataset 完成一轮有记录的本地、原始数据来源、原论文/作者资料检索即可。下载失败或资料找不到时记录具体尝试，继续其他模块；不无限重复检索。

### M01B. 证据必须分层

生成 `pathology_source_inventory.tsv`，至少有：

```text
source_id, dataset_id, section_id, accession, source_title,
source_url_or_local_path, exact_sample_mapping_evidence,
annotation_level, annotation_format, image_resolution,
coordinate_system, registration_evidence,
source_author_or_reviewer, publication_or_file_date,
applicable_region, evidence_status, limitation
```

`annotation_level` 区分：patient、specimen、whole_section、region/polygon、spot。整张肿瘤来源切片的身份，不可直接下放成所有上皮 spots 都恶性。

- N/T 文件名、Epi 高、diploid/aneuploid 任一个都不能独立充当病理审定。
- CNV 是辅助证据；diploid 不自动判正常，aneuploid 不自动判病理确诊肿瘤。
- 论文示意图可以是线索，不能自动转成精确 spot mask；不能从低分辨率截图伪造精确肿瘤轮廓。
- 有真实的原作者标注即可复用，但记录 `author_annotation`，不假装用户或病理医生已经再次审阅。
- 不由 Codex 代替病理人员签发最终诊断；AI/算法推断与作者标注、真实人工审阅分列保存。

若存在可直接映射的 annotations，先验证 barcode、section、坐标系、图像缩放和旋转/翻转，再生成本轮 `spot_pathology_overlay.tsv`。聚类原始标签不改动。

原有大类 `other` 如有资料，应进一步记录 stromal、immune、muscle、necrosis、lumen、artifact、off_tissue；信息不足保留 other，不全部改名 stroma。

### M01C. 只把剩余问题交给用户

复用已有 133 项优先表，不要求逐个审完 8,257 个连通块。输出：

- `pathology_resolved.tsv`：找到明确资料的区域和依据；
- `pathology_conflicts.tsv`：来源/形态/CNV/配准冲突；
- `ROI_review_to_fill.tsv`：仍需人工判断的条目；
- `WHAT_USER_NEEDS_TO_DO.md` 与可离线打开的 `review_index.html`。

用户说明必须明确：

1. 先看哪几个尚未解决的 ROI；每项涉及多少本轮 anchors，以及缺什么证据。
2. 对应哪个实际可打开的 H&E、ROI 叠加图及局部放大图，全部使用**相对链接**；不要保留只能在服务器点击的绝对路径链接。
3. 需要判断的是 `tumor_epithelium`、`non_tumor_epithelium`、`other`、`mixed_needs_split`、`uncertain` 之一。
4. 填写 `reviewed_label`、`review_comment`、`reviewer`、`review_date`；未看是空白，已看仍不能判断是 uncertain。
5. 混合区域可以提供 polygon 或 spot_uid 成员；记录 parent_roi_id、图像文件、像素坐标基准、scale factor/变换，不强迫整块一个标签。
6. 原始图分辨率不够时，要指出具体 section 缺少哪个文件；不笼统要求用户重新提供所有原始数据。
7. 不让用户寻找公共资料中已经存在的内容，也不伪造审阅者、日期或判断。

**M01 完成一轮检索和整理后，即可继续 MP10。病理未知不是评分失败，不得设置 `all_ROIs_reviewed` 全局开关。**

---

## M02. 中央 spot 的整体 MP10：先复用冻结分数

### M02A. 不重新定义 MP10

主评分继续使用原冻结 **MP10 Top20 AUCell**。不要换成 UCell、Core19、Top50、K32、NNLS 或重新挑代谢基因。其他历史评分不是本轮可自动替换的 primary。

从 `MP10_score_assets.tsv` 定位每张切片的文件，并查明：

- 精确 gene-set 和 gene universe；
- 原始表达输入 counts / data / 其他层；
- 基因 ID 映射与重复基因处理；
- AUCell 版本、`aucMaxRank`、`normAUC`、ties/zero 处理、seed、是否做过空间平滑；
- 评分列是 raw score 还是 percentiles/z-scores/residuals。

读取冻结原始 score，按 `section_id + barcode` 或规范 `spot_uid` 一对一合并。验证重复键、barcode 后缀、命名转换、缺失和非有限值；记录转换表。不得只凭 barcode 跨切片 join。

正式字段：

```text
spot_MP10_frozen
spot_score_source
spot_score_definition_id
spot_score_gene_universe_hash
spot_score_status
```

原始评分存在时不因为新库默认设置变化而覆盖。抽取预先固定的小子集检查复现；若 ties/版本等导致无法精确重现，保留冻结值与审计结果，不为了“对齐”重新覆盖所有分数。

确实缺失的冻结评分，可以依据找到的原流程补算，保存为新文件并注明 `RECOMPUTED_EXACT_DEFINITION`。原参数无法恢复则明确标注 `RECOMPUTED_WITH_DOCUMENTED_SETTINGS`，不能与已验证冻结分数悄悄混为一种来源。

### M02B. 严格处理缺基因

重新核对 raw counts、历史评分矩阵、model genes 三层的基因可用性。区分：

- 基因在 feature list 中，但该样本表达全为零；
- ID/alias 映射问题；
- 原始矩阵可用，但建模或评分前被过滤；
- 原始数据确实没有该基因。

已知审计：`GSE238145_B5` 的 Top20 缺 `HMGCS1`。必须实际核验，不靠手填既有数字。若仍缺失：

- 保留原冻结分数作 `PARTIAL_SIGNATURE` 描述结果；
- 不补基因、不按 20/19 放大、不把缺失当检测到的零表达；
- 完整 Top20 primary 比较与部分 signature 敏感性结果分开；
- 报告受影响 niche/anchors/sections，不能因缺基因直接丢掉整类 niche。

输出 `MP10_gene_availability_reaudit.tsv` 和各层的 gene universe manifest。

### M02C. 主对象是中央 anchor

取该 anchor 自己的 MP10，不计算其 ring 内 MP10 均值来冒充中心。

不使用 Harmony 坐标计算 MP10，不对 score 再做 Harmony，不采用邻域 score smoothing。不乘/除 Epi fraction 来制造“校正上皮分数”。

病理尚未审定时，名称为“各 niche 的上皮富集中央 spot MP10”；不能改成“已确诊肿瘤上皮 MP10”。

---

## M03. 实际提取模型估计的上皮归属表达并评分

### M03A. 资产复用优先顺序

1. 优先检查 accepted 模型已保存的 Epi conditional expected expression layer；核对来源、观测 counts、gene/spot 顺序、参数摘要及生成代码。
2. 没有可验证层时，加载对应的 cell2location `v2_30type` 训练模型及训练 AnnData，从已有 posterior/guide 提取必要变量，生成 Epi conditional counts。
3. 必要时只在**已训练模型上**补做 posterior sampling/export，保存到本轮新目录；这不是重新训练。
4. 如果模型 checkpoint 缺失，但观测 counts、reference 和模型所需后验参数完整，可采用与官方实现一致的解析计算；必须通过小批次接口对照或等价性测试。
5. 只有归一化后的 q05 abundance 或 Epi fraction，而缺少 reference、背景等必要参数时，不足以重建可靠的归属表达。明确记录缺项和失败分支；仍完成其他切片和整体 score。不能凭一个矩阵乘法冒充官方输出。

不要误载旧的 32 类 reference、K32 模型或不对应本轮空间拟合的 cell2location 对象。

### M03B. 使用哪一种表达估计

优先使用已安装版本的：

```python
mod.module.model.compute_expected_per_cell_type(
    posterior_summary,
    mod.adata_manager,
    ind_x=spot_indices
)
```

这是接口意图示例，运行前检查本机函数签名、对象层级、registry 和返回结构，不盲抄历史 API。

必须区别以下两种量：

- `conditional_attributed_counts`：由实际观测 counts 按细胞类型贡献比例条件性分配得到；本轮 primary Epi 评分使用这一层。
- `reference_only_expected_mu`：细胞丰度 × 固定参考表达的模型期望；只有单个 Epi signature 时，变化可能主要是统一的幅度缩放，不能作为上皮内部连续状态的主分数。

参考官方实现，条件性分配的示意为：

```text
r[s,g,E] = w[s,E] * reference[g,E] * m[g]
           / (sum_f w[s,f] * reference[g,f] * m[g] + background[s,g])

Y_epi[s,g] = Y_observed[s,g] * r[s,g,E]
```

具体以 accepted 模型和安装版本的公式为准；如 detection 项在分子分母公共因子中消去，要在方法记录中说明。不能混用参考计数尺度、归一化比例和 log expression。

- 主点估计优先用完整的 **posterior means plug-in**，并准确标注，不称为严格的 posterior mean of attributed expression。
- 环境 q05 的选择保持不变；它不要求本模块只能使用 q05 参数。
- 需要 `w_sf`、reference、gene scaling、background、batch mapping 等一致信息，不能将不同运行的参数拼起来。
- 原保存层若用 q05，可留作 `legacy_q05_attribution` 对照；不与 mean-derived 层无标记混用。
- Epi 类型通过 accepted factor_names/type_mapping 精确选择。若原模型真的有多个 Epi factors，按真实映射合并其**表达 counts**；不能只凭名字子串选择。粗粒度 Epi 中无法再拆恶性/正常时如实注明。

### M03C. 先做技术验证，再自动扩展

先按每 dataset 一张、具备有效模型且 gene coverage 最完整的切片做小批次测试。选片不根据 MP10 高低；同等条件按 section_id 稳定排序，每片最多 128 个固定 anchors。

测试通过后，**自动继续全部有 anchors 且模型可用的切片**，不在 pilot 后再请求人工确认。

必须检查：

1. 表达层 spot/genes/factors 顺序与训练 registry 一致。
2. 输出有限、非负；分配给 Epi 的 counts 不超过同基因 observed counts（允许记录的浮点误差）。
3. 在小批次上检查全部 cell types 与 ambient 的贡献之和是否恢复 observed counts；ambient 不可被误计入上皮。
4. denominator 为零但 observed 非零时必须报数值/模型覆盖问题，不静默填零。
5. 若为减少内存只实现 Epi slice，需与官方全类型函数的小批次结果比较，目标 `rtol=1e-5, atol=1e-6`，并报告实际最大误差。
6. `n_epi_attributed_UMI`、Epi 非零基因数、Top20 检出数以及各 niche 的可计算覆盖率。
7. Epi score 是否几乎常数、是否仅由一个固定参考模板/单一幅度主导；异常照实报告，不通过换基因强行制造异质性。

正式输出至少包括：

```text
epi_model_id, epi_attribution_method, posterior_summary_type,
epi_factor_names, epi_total_attributed_counts, epi_n_detected_genes,
epi_top20_detected_genes, epi_score_valid, epi_score_limitation
```

### M03D. 匹配基因背景后计算两个分数

对每张切片，冻结模型真正覆盖、能够与 observed counts 对齐的广泛基因集合 `G_s`。

在**同一 `G_s`、同一 signature、同一 AUCell 配置**下分别计算：

```text
spot_MP10_model_universe   # Y_observed 在 G_s 中的 score
epi_MP10_model_universe    # Y_epi 在 G_s 中的 score
```

这两列构成可以配对比较的结果；原 `spot_MP10_frozen` 继续保留，不能被覆盖。

**禁止只用 Top20 的 20 行构建 AUCell rankings。** 必须使用完整且可追溯的表达背景，Top20 只是被评估的 gene set。

新增 matched-universe AUCell 的设置：

- 采用与原 AUCell 尽量一致的归一化、rank cutoff 比例及 ties 规则，并写入配置。
- 原 rank cutoff 比例确实无法恢复时，新匹配分数显式使用 `aucMaxRank = ceiling(0.05 * |G_s|)`，`normAUC=TRUE`，保存其独立 definition ID；不将其冒充原冻结分数。
- `G_s` 内的 observed 与 Epi 使用相同 gene order、rank cutoff、signature 和 zero/tie handling。保存 seed/分块方案，防止不同 batch 排列改变 ties。
- 不把小数型 attributed counts 随意四舍五入为整数，不为每个 niche 单独筛基因。
- 对全零 Epi 表达返回 score NA，而不是随机零值排名生成的“有效分数”；分母为零也用 NA。
- 对检出基因数低于 `aucMaxRank` 等情况记录低信息警告；不要因真实低 MP10 或 signature 检出少就事后剔除。

如果 `G_s` 缺少 Top20 某基因：

- 完整 Top20 模型分数记为不可用，并列明缺因；
- 可以另出显式标记的 matched partial-signature 对照：observed 与 Epi 都使用完全相同的可用集合，并记录 `signature_used_hash`、缺失列表和 `PARTIAL_SIGNATURE`；
- 不给该基因编造 Epi expression，不拿其他基因替补，不与完整20基因的结果无标记合并；
- 不为了迁就一张缺基因切片，将所有切片主 signature 改成 Core19。

不同切片 `G_s` 不同的情况下，跨切片结论优先使用同切片的对比及方向一致性；不要将原始绝对 AUC 全部混成唯一排名。不把两层 AUC 差值称为“非上皮贡献比例”。

### M03E. 同时提供直观的 counts 来源摘要

以模型覆盖的 signature 为基础，给每个中心计算：

```text
signature_epi_count_share_of_observed
signature_nonepi_count_share_of_observed
signature_ambient_count_share_of_observed
allgene_epi_count_share_of_observed
```

这些是**模型估计的 RNA counts 分配比例**，不是细胞比例，也不是 AUCell 分数百分比。signature observed 总 counts 为零时记 NA。按 niche、section、dataset 汇总，并同时展示 gene coverage。

低 signature counts 的值不稳定，应显示信息量，而不是用一个任意的 50%/80% 阈值宣布上皮专属。

### M03F. 小范围 posterior 不确定性检查

优先复用已保存的联合 posterior draws。条件允许时在上述技术 pilot 的固定 anchors 上使用约 20 个配套 draws，逐 draw 先分配表达再评分，检查 plug-in 分数与 draw 分布。

- 局部和全局变量需来自配套的 posterior sampling/export；不把不同变量的独立 q05 拼成“一个 posterior draw”。
- 不保存全量高维 draw 张量；保留分数摘要。
- draw-level 结果不等于真实生物学的不确定性，也不等于独立验证。
- 无法完成时单列 limitation；不阻塞全部切片的可用点估计。

---

## M04. 统一汇总九个 niche 的中央分数

### M04A. 必须生成完整中心级结果表

`anchor_MP10_all_layers.parquet` 每个共同 anchor 一行，至少含：

```text
spot_uid, barcode, dataset_id, section_id, sample_uid,
patient_uid, analysis_unit_uid, analysis_unit_type,
source_patient_identity_known,
niche_id, niche_name, accepted_label, candidate_label,
assignment_status, candidate_vote_fraction,
epi_fraction_center_mean, complete_rings,
center_nCount, center_nFeature,
spot_MP10_frozen, spot_MP10_model_universe, epi_MP10_model_universe,
score_definition_ids, gene_universe_hashes,
n_signature_available_raw, n_signature_available_model,
signature_coverage_status, epi_score_valid,
signature_epi_count_share_of_observed,
pathology_label, pathology_evidence_level, pathology_source_id,
missing_reason
```

没有的字段记 NA 和原因；不存在的指标不能生成貌似合理的数值。保留患者证据层级，不能将 section_id 自动填进 patient_uid。

### M04B. 分析集合分开

至少包含：

1. **ALL_ACCEPTED**：正式接受标签的上皮 anchors，允许病理未知；按 score/gene availability 标记覆盖。
2. **FULL20_SPOT_PRIMARY**：整体分数的完整 Top20 集合。
3. **MATCHED_FULL20_EPI**：同一批中心且整体/上皮 matched-universe 分数均可用、完整 Top20 的集合。
4. **PARTIAL_SIGNATURE_SENSITIVITY**：缺基因结果单独描述，不能偷偷混入2/3。
5. **PATHOLOGY_SUPPORTED_TUMOR / NON_TUMOR**：仅在已有可追溯证据时分层；不能把 unknown 填正常。

比较整体与 Epi 层时，先用完全相同的 centers 重汇总两层；不能将“整体分析几十张切片、Epi 分析几张切片”的差异归因为纯化效果。

### M04C. 以样本内效应为核心，先看全部 niche

所有 niche0–niche8 均输出均值、中位数、IQR、零值比例、可计算 anchors、切片数、已知患者/独立单位数、中心 Epi 比例和测序深度。固定列顺序 niche0–8，不按 MP10 高低重命名。

主展示优先：

- 每切片 × niche 的原始 score 中位数；
- 同切片内该 niche 与其余 accepted niches 的 score 差异；
- 每 dataset 的方向一致性与效应分布；
- 同一患者多切片的归并结果。

九个 niche-versus-rest 比较是九个预先指定的描述性对比；其“rest”组成随切片变化，必须说明。无需默认展开全部36对 pairwise 显著性检验。

初拟 `MIN_ANCHORS_PER_SECTION_NICHE=20`，对照集合也至少20，用于稳定汇总的工作门槛；更小者仍显示分布和数量，标记为低覆盖。此门槛在查看 MP10 结果前写入配置，不按显著性修改。

整体显示可以另提供 section 内百分位/标准化分数，但必须另列；不要覆盖 raw score，不把所有切片强行调整成相同生物分布。

若计算置信区间或显著性：

- 以已核实患者或有独立性证据的 analysis unit 为重采样单位；同患者多切片捆绑。
- 不进行全体 spots 一次 Wilcoxon 来假装独立重复；不把重叠 rings 当独立病例。
- 至少5个可比较独立单位才提供探索性汇总区间/检验；更少者只描述。此5是本轮工作规则，不是普适统计标准。
- 主推断以患者/单位的切片内差异为基础；需要切片内区间时采用空间块方法，而不是逐 spot iid bootstrap。
- 可预设1000次单位级bootstrap，固定seed。涉及九个planned contrasts的P值按对应score family进行BH，并报告实际检验数量；没有可靠检验就不填P值。
- 不把不同score family中最显著的一个挑作唯一结果。

增加一个预设的中心混合/技术因素调整结果：在有覆盖的切片内比较 niche，考虑中央 mean Epi fraction 与测序深度；可用固定效应 `section + niche + spline(epi_fraction, df=3) + log1p(nCount)`，不同时堆入高度共线的nCount与nFeature。模型区间不得用 iid spot 标准误充当患者证据；样本少、支持区间不重叠或模型退化时只报告不可估计及原因。

### M04D. 有限敏感性检查，不重新聚类

沿用相同 labels 检查：

- 完整 6+12 邻居子集；
- 中央 mean Epi >=0.5；
- Epi >=0.7 仅作补充，同时报告九群覆盖损失；不能让它替代全部 niche 的上皮来源分析；
- accepted 与 all-candidate 对照；
- 去掉单个贡献最大的 dataset 后的描述性结果。

不因某个子集没有 niche6 或效果改变就再改阈值。某群不存在或数量不足是 missing/not estimable，不是零效应。此处删除来源后重汇总不称为“重新训练的外部验证”。

---

## M05. 本轮上皮来源评估：复用旧结果，补当前九群的关键对照

上皮归属分数属于模型结果；本模块是解释可信度检查，不保证一定得出“上皮来源”。先完成 M02/M03 的实际分数，来源模拟与汇总并行，不能用模拟尚未全部完成阻塞交付。

### M05A. 真实单细胞 Epi 是独立的一层参考

复用真实 scRNA counts 和既有细胞类型注释，按原 Top20 评估 Epi 与其他谱系的 score、检出率和参考表达分布。区分正常/恶性上皮的前提是已有可追溯身份依据；缺失时只写 Epithelial。

这是真实单细胞上的上皮评分，但不等于已经定位到某个 Visium 中央 spot。scRNA 与空间样本若未配对，不得把单细胞硬配给具体 niche。

已有 NMF MP10 usage 可用于检验 Epi 内部评分是否保留状态差异；它是复用原MP定义，不是本轮重建NMF。不具备时清楚说明。

### M05B. composition-only 对照：固定上皮状态，改变中央混合组成

复用现有 counts-first pseudo-spot 工具，针对当前 niche0–niche8 补设计，而不是复跑旧六方案。

关键定义：

- 模拟的是 **中央 spot 的组成**，不是 ring 组成。
- 中央组成来自已有 c2l posterior mean / 原 counts-compatible 组成估计，不能用环 q05 比例代替。
- 各 niche 使用相同的上皮供体/状态采样分布，改变实际中央非上皮组成及其占比；这是“无上皮状态差异”的模型条件对照。
- 先在 counts 层混合、匹配实际文库深度和基因范围，再计算 AUCell；不平均单细胞 AUC 代替 mixed counts 评分。
- RNA贡献与细胞数量不同，记录细胞大小/文库量匹配方案，不能把 cell fraction 直接当成RNA fraction。
- 供体信息、随机种子、模拟组成、原始上皮counts和混合counts来源均可追溯。

运行一个有限但覆盖九群的批次：最多500个中央组成蓝本，按 dataset×niche 分层兼顾 Epi 比例，每蓝本5次固定种子重复，目标不超过2500个 composition-only pseudo-spots。小数据集/小niche尽量覆盖，不强制每格等数；报告实际支持。

检验：只靠中心混合组成，是否可产生与实际类似的 niche–spot MP10 趋势？给出 effect/direction 和模拟分布，不把模拟本身当独立患者或实验重复。模拟未复现某效应不自动证明其来自上皮。

### M05C. 状态变化正对照与“纯上皮”真值

在较小的代表蓝本集合（例如最多100个）中固定中央组成、非上皮counts与总测序尺度，替换上皮内部较低/较高MP10状态（优先使用已有NMF usage；没有则说明按纯Epi评分分组用于评分技术校验）。保留已知的纯Epi counts及其同背景score，测试整体分数的状态敏感性。

如可把当前归属公式用于这些模拟，再比较 recovered Epi 与真实 Epi 的score。必须说明：

- 用真实模拟丰度输入归属算法时属于 **oracle composition** 技术校验，不能冒称已验证真实 Visium 丰度估计的准确性。
- 如果模型比例不是在模拟counts上独立估计，就不能宣称这是完整的端到端来源验证。
- 不为这一项重新训练63张空间模型；缺少可行估计时如实保留技术校验层级。

### M05D. 每个 niche 的来源证据汇总

对九群统一输出：

```text
niche_name
spot_effect_and_coverage
epi_effect_same_centers_and_coverage
matched_universe_status
signature_count_attribution_summary
single_cell_support
composition_only_result
state_change_control_result
high_epi_sensitivity
pathology_scope
source_evidence_level
unresolved_reason
```

允许的结论用保守、可核查的措辞，例如：

- `MODEL_ASSISTED_EPI_SUPPORT`：在模型归属Epi层也保留差异，且对照未显示其仅由已检验的组成因素产生；仍然是模型支持。
- `COMPOSITION_CONFOUNDED_OR_ATTENUATED`：整体差异在来源/组成分析后减弱或可被组成变化解释。
- `DISCORDANT`：整体与Epi归属层的方向或结果不一致。
- `UNRESOLVED`：基因覆盖、模型、独立样本或信息量不足。

不要将任何一项自动输出成 `PROVEN_EPITHELIAL_SPECIFIC`，也不要预设niche3或niche6必须有高MP10。rank-based score 的来源不能通过线性分数相减确定。

---

## M06. 图表与输出：按 dataset 一级组织，同时提供全局九群视图

### 6.1 必做图表

所有图固定显示 niche0–niche8；无数据的群保留空位/NA和说明。至少生成：

1. **九群中央整体 MP10 分布**：violin/box + section/患者汇总点；spot展示可平衡下采样，但汇总使用全部合格数据。
2. **九群上皮归属 MP10 分布**：布局与整体一致，图注明确 model-attributed、gene universe和覆盖。
3. **matched-universe 两层对照**：同一批 centers、同一 gene background 的整体 vs Epi，分别画score分布与niche effects；不比较不同背景原始AUC的“绝对纯化增减”。
4. **section × niche 热图**：raw median 与 section-relative 描述分开；旁边标注dataset、患者信息层级、N和缺基因/缺模型。
5. **dataset/患者层面效应图**：所有九群，不只挑3和6；效应、方向和覆盖齐全。
6. **来源 counts 归属摘要**：按niche显示Top20 observed RNA被估计归给Epi、非Epi和ambient的比例与信息量。
7. **覆盖及敏感性表**：完整Top20、matchedEpi、高Epi、完整环、病理子集分别有多少中心和样本。

已有原始组织坐标可用于新增 score 空间图，但这是评分可视化，不恢复被用户取消的代表切片人工复核任务；不以人工逐图验收阻塞本轮。

沿用已有颜色与图例规范、透明背景，连续空间score采用plasma且0黑（NA与0分开）；输出PNG+矢量PDF。报告原图类型、图注和单位。不要把两个不同量纲的score强制共用色阶。

### 6.2 目录建议

```text
neighborhood_mp10_source_v4/
  00_audit_and_freeze/
  01_pathology_search/
    sources/
    user_review/
  02_spot_MP10/
  03_epi_attribution/
    by_section/
  04_niche_comparison/
    by_dataset/
    all_datasets/
  05_source_validation/
  06_figures/
  code/
  logs/
  reports/
```

### 6.3 核心表与报告

至少有：

```text
00_audit_and_freeze/resolved_inputs.tsv
00_audit_and_freeze/niche_name_map.tsv
00_audit_and_freeze/analysis_authorization.json
00_audit_and_freeze/analysis_spec_frozen_before_scoring.yaml
00_audit_and_freeze/MP10_gene_availability_reaudit.tsv
01_pathology_search/pathology_source_inventory.tsv
01_pathology_search/pathology_resolved.tsv
01_pathology_search/pathology_conflicts.tsv
01_pathology_search/spot_pathology_overlay.tsv
01_pathology_search/user_review/WHAT_USER_NEEDS_TO_DO.md
01_pathology_search/user_review/ROI_review_to_fill.tsv
01_pathology_search/user_review/review_index.html
03_epi_attribution/attribution_technical_validation.tsv
03_epi_attribution/model_and_gene_universe_manifest.tsv
04_niche_comparison/anchor_MP10_all_layers.parquet
04_niche_comparison/niche_MP10_summary_all9.tsv
04_niche_comparison/section_niche_score_summary.tsv
04_niche_comparison/independent_unit_effects.tsv
04_niche_comparison/dataset_niche_effects.tsv
04_niche_comparison/coverage_sensitivity.tsv
05_source_validation/source_evidence_all9.tsv
reports/METHODS.md
reports/RUN_SUMMARY.md
reports/REVIEW_INDEX.md
reports/WHAT_IS_STILL_UNCERTAIN.md
reports/review.zip
code/resume_mp10_source_v4.sh
```

大表用parquet，小表TSV。轻量 review 包包含报告、汇总表、关键图、确实未解决的用户操作说明；不打包全部原始counts、模型、庞大Epi表达矩阵和高分辨率H&E。

---

## 7. 本轮验收与遇到问题的处理

必须完成的实际行为：

- 病理本地/公共来源检索有真实日志和证据，不只是给用户一个搜索计划。
- 冻结九群并连接中央spot MP10，所有niche均有覆盖统计和图表。
- 对已有模型实际尝试Epi conditional expression提取，先小批次验证，成功后扩展；不能只写一份未来方案。
- 输出matched-universe的spot/Epi配对评分以及基因可用性。
- 复用或补做当前可行的来源检查，区分技术校验、模型归属和独立生物学证据。
- 若某分支不能完成，逐section列真实原因，提交已经完成的评分/图，而不是只输出失败结论。

下列情形的处理：

| 问题 | 处理 |
|---|---|
| 病理资料没有找到 | unknown + 已查来源 + 具体用户待办；继续评分 |
| 单切片模型缺失或无法恢复 | 完成其spot分数，Epi列NA和具体缺项；继续其他切片 |
| model genes缺Top20 | 完整模型分数不可用；同背景partial明确分开，不替补 |
| 所有模型均缺可恢复参数 | 交付实际spot分析和来源资产诊断，不能伪造纯Epi；列最小补充文件 |
| 模型Epi分数几乎常量/明显模板主导 | 输出真实结果和诊断，不再换基因/重调niche追求漂亮差异 |
| 独立样本或两群共同覆盖不足 | 描述性结果 + not estimable；不制造患者数和P值 |
| 部分模拟暂时不可运行 | 标记具体未验证内容，不阻塞两层评分结果 |
| 磁盘/内存受限 | 降批次和并行度，优先落盘score与小表；不删除旧原始结果 |

续跑脚本必须有明确范围和幂等性。旧结果只读；本轮写入检查实际完成的产物，不提前创建总体成功标记。

---

## 8. 执行结束给用户的回复格式

请按以下五部分作答，提供真实绝对路径和已创建的轻量包：

**A. 实际完成范围**

- 病理检索覆盖、找到的spot/region标注数、仍待审阅数。
- 整体score和Epi归属score分别覆盖多少切片、anchors、niches。
- 是否已完成matched-universe与当前九群来源对照。

**B. 九个 niche 的结果表**

每个niche列出整体MP10、匹配背景EpiMP10、样本内方向、独立单位支持、缺基因/模型限制和来源证据级别。不要只报P值或只展示结果最好的群。

**C. 用户现在必须做什么**

先说明哪些病理问题已经通过资料解决。仅列剩余明确事项：打开哪个相对链接、判断什么、填哪些列、回传哪个文件。没有必要的人工任务就直接写无新增必需任务，不制造额外审核。

**D. 限制与失败项**

明确“中央spot高分”和“模型估计上皮高分”的区别；明确尚未证明的来源/恶性身份/因果。不可把来源验证部分完成写成上皮专属已证明。

**E. 结果入口与续跑**

给出 `RUN_SUMMARY.md`、九群汇总表、两层主图、病理待办、`review.zip` 及真实可运行的续跑命令。本轮结果完成后停止，不自动进入通讯、机制、临床或新聚类。

---

## 9. 方法与来源登记

请核对当前安装版本的官方实现，方法描述区分原方法支持的步骤与本项目的工作阈值。

### 现有项目证据（已审阅 commit）

- 九群配置：
  https://github.com/Kimsx05/transfer_for_neibor/blob/a785e2710565ac31c586c5e1d09cfabfd7277b5b/results/neighborhood_ring_harmony_v3_upload/02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/working_config.json
- 病理复核说明：
  https://github.com/Kimsx05/transfer_for_neibor/blob/a785e2710565ac31c586c5e1d09cfabfd7277b5b/results/neighborhood_ring_harmony_v3_upload/05_user_review/WHAT_USER_NEEDS_TO_DO.md
- 既有缺基因审计：
  https://github.com/Kimsx05/transfer_for_neibor/blob/a785e2710565ac31c586c5e1d09cfabfd7277b5b/results/neighborhood_ring_comparison_v2/neighborhood_ring_comparison_v2/00_preflight/MP10_gene_availability.tsv

### 本轮评分和表达归属的官方依据

- cell2location 官方教程中 “Estimate cell-type specific expression of every gene in the spatial data”：
  https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
- `compute_expected_per_cell_type` 实现（按本机版本再核对）：
  https://cell2location.readthedocs.io/en/latest/_modules/cell2location/models/_cell2location_module.html
- AUCell 官方 vignette（rankings、aucMaxRank、检测基因数和ties注意事项）：
  https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html

本轮不需要重新证明Leiden方法合理或重做整合benchmark；主要任务是上述病理资料检索、中央两层MP10实际计算及来源评估。
