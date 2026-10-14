# Codex 执行规格 v7：获批宽口径候选的独立 R1/R2/R3 分析

## 0. 本轮授权：现在执行阶段 B，不再停在旧 Gate

你正在膀胱癌单细胞与常规 Visium 空间项目中工作。完整阅读本文件及随附清单后，**实际完成计算、统计、出图和报告**。

用户已审阅 v6 的生物学解释，并明确同意宽口径推进：60 个有单位层面正向支持的非上皮组合、39 个补充探索组合、40 个负向对照组合，共 **139 个不同 niche–celltype 组合**。本轮授权包含总体和全部七个 dataset 的分层输出，不按是否显著或是否符合预期选择 dataset。

随附文件：
- `approved_pairs_v7.tsv`：139 个不同组合，保留原总体 candidate_id 作为 master_candidate_id。
- `approved_gate_v7.tsv`：上述组合在总体及七个 dataset 中展开的 **1,112 个 GO 行**；使用 v6 各 scope 的真实 candidate_id。
- `niche_celltype_review_suggestions_270.tsv`：原评审建议，只读留档；其中原 decision 为空不影响本次明确授权。
- `provenance_v7.json`：来源哈希、审批方式、数量和检验族规模。
- `verify_v7_approval.py`：只校验清单，不执行实际分析。

**这些 GO 行是根据用户对话中的整体批准生成的，不声称用户逐行手填或完成了病理审定。** 保存新的 gate receipt，不修改原 v6 审阅表。latest v7 授权覆盖 v6 中“尚未允许执行 B”的旧限制；其他冻结定义保留。

重要：已检查的 v6 `--phase B` 仅校验 Gate 并返回 `GATE_VALIDATED_PHASE_B_PENDING_IMPLEMENTATION`。**不能只运行这个接口就报告完成。必须在新目录实现并执行真正的 R3、三箱图及连续 MP10 关联。**

### 本轮执行范围

1. 原始 q05 环组成的独立 R1/R2/R3。
2. 139 组合的原始比例、相对同切片背景的差值，以及三个环之间的配对比较。
3. 所有获批 celltypes 在“全体合格上皮中心”和“指定 niche 内部”的连续 MP10 关联，使用已有两层匹配评分。
4. 上皮/TME 总份额的结构性解释；一个有限的 TME 内部组成敏感性检查。
5. 全体及每个 dataset 的三箱图、完整数值结果、来源和统计限制。

不开展：新聚类/整合、模型训练、评分优化、上皮来源恢复、C-SIDE、SpottedPy、LIANA+、热点、距离软件、细胞通讯、预后或机制分析。无需等待全量病理审定；不新增“代表切片空间复核”人工任务。

---

## 1. 输出目录、输入定位和只读保护

项目根目录：
```text
/home/data/t070721/codex_workspace/Bladder Metabolism
```
新输出目录：
```text
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring123_broad_v7
```

只读复用：
```text
neighborhood_niche_celltype_screen_v6/
neighborhood_ring_comparison_v2/
neighborhood_ring_harmony_v3/
neighborhood_mp10_source_v4/
neighborhood_mp10_calibration_v5/
```

优先读取 v6 `00_audit/input_manifest.tsv`、`reports/METHODS.md`、`00_audit/analysis_unit_policy.tsv`、`00_audit/spec_frozen.yaml` 和原物理图生产代码，不从名称猜矩阵轴。经 v6 记录的关键来源包括：
```text
neighborhood_ring_comparison_v2/03_features/common_anchors.parquet
neighborhood_ring_comparison_v2/03_features/features.npz
neighborhood_ring_comparison_v2/03_features/feature_columns.json
neighborhood_ring_comparison_v2/03_features/COMPLETE.json
neighborhood_ring_comparison_v2/03_features/anchor_master.parquet
neighborhood_ring_comparison_v2/code/01_prepare.py
neighborhood_ring_comparison_v2/02_graphs/<dataset>/<section>__R1.npz
neighborhood_ring_comparison_v2/02_graphs/<dataset>/<section>__R2.npz
neighborhood_ring_harmony_v3/01_features/metadata.parquet
neighborhood_ring_harmony_v3/02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/labels.parquet
neighborhood_mp10_source_v4/04_niche_comparison/anchor_MP10_all_layers.parquet
neighborhood_mp10_calibration_v5/00_audit_and_spec/identity_and_independence.tsv
```

全组织 q05 组成的具体路径，从原 manifest/生产脚本解析；常见位置为：
```text
neighborhood_v1/00_inputs_and_audit/subtype_q05/<section>__subtype_q05.csv.gz
```

若上传目录名称含 `_upload` 或解压后重复嵌套，先核对内容、哈希及 UID，再定位实际来源。不要把服务器 `/home/...` 与聊天附件 `/mnt/data/...` 当成同一文件系统。只在上述根目录及用户上传目录内有界查找。

先读取运行进程与磁盘/配额状况，避免与正在运行的任务冲突。不终止他人作业，不删除旧成果，不重训模型。只加载必要的元数据、软组成和已算分数，不加载全部基因表达或大模型。

新随机操作主 seed=123；派生子种子用稳定哈希生成并写入日志，不使用 Python 随机化的内置 hash 作为永久种子。

---

## 2. 清单核验与分层意义

### 2.1 数量契约

| review_tier | 组合数 | 使用方式 |
|---|---:|---|
| POSITIVE_CORE | 28 | 优先展示，但所有推断仍为探索性 |
| POSITIVE_BROAD | 32 | 有原单位正向支持，完整运行 |
| EXPLORATORY | 39 | 原支持较弱或异质，完整运行，不写成“已富集” |
| DEPLETION_CONTROL | 40 | 负向对照，完整运行；不预设新结果必须负向 |
| 合计 | **139** | **29 个不同非上皮类型；niche1–niche8** |

所有新检验都双侧。review_tier 只影响图表顺序和证据注释，不决定检验方向，不为每层单独缩小 BH family。

`RESERVE`、`NICHE0_BACKGROUND`、`EPI_CONTROL` 不进入新增的 139 组机制候选。原 270 项证据完整保留。

niche0 不删除：仍属于全体合格上皮中心与背景池，作为结构性参照显示 Epithelial/TME 总份额，但不生成29条“niche0回避某细胞”的候选分支。其余 niche 的 Epithelial 也作为组成解释，而非新增候选。

### 2.2 授权 receipt

核验：139 个 master 组合唯一；展开后1,112个 candidate_id 唯一；每个展开行在原 v6 Gate 中的 `scope/niche/celltype/accepted_label` 全部一致；每个 master 有8个scope；日期有效；GO字段不为空。

写 `00_approval/approval_receipt.json`，记录对话授权、清单及原 Gate SHA256、各层数量、`manual_row_review=false`，并说明没有病理审定含义。可以复制授权清单到新目录，不回写旧 Gate。

若本机缺少原 Gate，但原 v6 冻结证据表和清单键可核对，使用这些表验证并记录。清单与来源内容冲突时只阻塞对应的身份不明组合，不猜 ID；继续其他已核实组合。目标群和别名不得私自替换，例如 SMC、Pericyte 与 CAF 保持区分。

---

## 3. M00：冻结输入与统计口径

保持：
- `SAMPLE_R12_LEIDEN` 的 accepted labels 0–8，显示为 niche0–niche8。
- `accepted_label=-1` 只留在审计，不作为目标、rest 或主要背景。
- q05 按每个 spot 全部30类型归一化，再对环内有效 spots 等权平均。
- 中央 gate 为既有 posterior mean Epithelial fraction≥0.25；不改成 q05 gate。
- 30类列序固定，所有矩阵经 `spot_uid`/barcode 明确连接。
- 多切片按既有 `analysis_unit_uid` 归并；已核实患者、来源支持匿名样本和其他既有单位政策分别报告，不合并改称患者。

已知旧规模仅用于核验：129,632共同anchors、126,075正式接受；旧完整双环108,345指共同anchors，不是新三环中心数。63张目标切片中61张原有共同anchors。R3后的实际保留数必须重新统计，不能填旧数。

输入有版本冲突时报告差异，不静默选更“好”的一版。新配置、阈值、主检验族在读取本轮关联输出之前保存为 `00_audit/spec_frozen.json`。

---

## 4. M01：只扩展独立 R3，重用并校验 R1/R2

### 4.1 物理图必须来自完整组织

每片独立处理：
- 使用原真实六角邻接图，不能用 Harmony、UMAP 或表达近邻代替物理图。
- 保留 tissue=1 但没有组成输出的物理节点在拓扑中。它们在组成统计时是 missing，不是全零细胞。
- 不先删成只含上皮节点的图，不跨切片，不跨组织外/孔洞添加捷径，不用更远邻居凑足数量。
- Rr 是在这个图中与中心最短距离恰为 r 的独立环。中心排除，各环去重且互斥。
- 完整规则格点的 R1/R2/R3 理论数量为6/12/18；边缘缺口保留真实可达结果，异常结构记录审计，不裁剪邻居以伪造标准。

稀疏矩阵/有限深度 BFS 即可；避免构造全切片的稠密所有点对最短路径矩阵。复用 R1/R2 成员表，新增 R3。

### 4.2 主三圈中心集合

对每片 s，定义 S_s：既有接受标签的上皮 anchors，且三圈有效组成覆盖同时满足：
```text
R1: ≥4/6
R2: ≥8/12
R3: ≥12/18
```
三圈主比较使用完全相同的 S_s，不能每个环各用不同的一批中心。组成数据一行若存在未解决类型缺失，不把个别类型补0；依照既有行级有效性定义处理，保证同一环30类的分母一致。

完整6/12/18邻域作为预定敏感性；缺少R3不反过来覆盖v6的A结果。报告每个dataset/section/niche在各步的保留率和缺失原因。

### 4.3 环组成

记邻居j的原始相对组成为p_jc：
```text
q[i,r,c] = sum(p[j,c] for composition-valid j in exact ring r) / valid_count[i,r]
q12[i,c] = (q[i,1,c] + q[i,2,c]) / 2
```

q12用于与v6衔接，不将其改为18邻居的统一平均。R3不要附加到旧60维再聚类。

抽查每个dataset至少一片，并全量验证可比较行：重新汇总q1/q2应在原float32精度允许误差内匹配旧缓存，主容差1e-6；超差明确定位原因，不以重新生成旧矩阵掩盖差异。每行30类之和应约为1。

保存矩阵轴、环成员/稀疏图、有效邻居数与哈希。后续139组合共用这些缓存，不为每个celltype重新找邻居。

---

## 5. M02：同切片背景与单位摘要

### 5.1 背景保持用户同意的定义

每片、每圈、每个celltype：
```text
b[s,r,c] = mean(q[i,r,c] over all i in S_s)
```

背景包含niche0以及该片其余接受niches，不是整片所有组织spots均值，不是跨dataset混合均值，也不依据MP10高低挑背景。

目标T_sk为S_s中niche=k的中心：
```text
m[s,k,r,c] = mean(q[i,r,c] over T_sk)
e[s,k,r,c] = m[s,k,r,c] - b[s,r,c]
rest[s,k,r,c] = mean(q[i,r,c] over S_s excluding T_sk)
d_rest[s,k,r,c] = m[s,k,r,c] - rest[s,k,r,c]
```

e 为主背景差，d_rest 为敏感性，不能混为一个指标。保存 target 在背景中所占比例 lambda=n_target/n_background；校验均值恒等式：
```text
e = (1-lambda) * d_rest
```
背景包含目标自身会收缩效应，不能把它重复复制成若干“独立对照”。

### 5.2 单位聚合

沿用v6的可比较条件：同片目标≥20个中心、rest≥20个中心。三圈使用同一可比较section集合。没有rest/不足20仅作描述，不能从背景差0宣布无效应。

先每片计算m、b、e、d_rest，再在同一analysis_unit内对该组合可用的切片等权平均。总体汇总时同一真实analysis_unit跨dataset出现也只计一次；每dataset内部仍可分别展示其切片，明确重复来源。不得为了方便按 dataset::patient 再拆同一真实患者。

同时保存现有点全部的描述摘要与满足覆盖条件的推断摘要。主箱体/星号只使用同一推断摘要，图上N须与检验N一致。

---

## 6. M03：三个并列箱线图和 Wilcoxon

不画距离曲线、不连单位点、不以折线强制表达单调趋势。

### 6.1 每个获批组合、每个scope必须有的输出

1. `raw`：R1/R2/R3三个并列箱体，y为单位层目标环比例（%）。
2. `excess`：同样三个箱体，y为单位层e×100（百分点差），水平0线为样本背景。
3. `ringpairs`：仍为三个并列箱体，用单独括号标R1−R2、R1−R3、R2−R3的比较；不混用excess的对0星号。

单位点叠加，颜色/形状可区分已核实与匿名单位，但不制造新的分组检验。数据少仍显示实际点和N；无数据使用统一NA占位图/索引状态，不能伪造全零箱体或静默删除scope。

箱体中线=median，箱体Q1/Q3，须至1.5×IQR范围内最远点，离群点只画一次。这个外观借鉴用户提供论文；我们的配对统计不是原论文Fig3d的unpaired t-test复刻。

### 6.2 对背景0的星号

以单位层e为观测，双侧one-sample Wilcoxon signed-rank（等价于配对目标与背景），不是anchor层两独立样本rank-sum。

严格复用v6的数值规则，并做R/Python小例子交叉核验：
- ≥5有效单位且≥5非零差值才报告常规推断。
- ≥5单位、全部差值为0，p=1；不是缺失。
- 少于要求显示 `N/A`，不是 `ns`；实际描述效应仍输出。
- 去零采用Wilcox zero exclusion；n<50且无零值/绝对秩并列时exact；其余采用并列校正的渐近方法及连续性校正。
- 不为避免ties人为加jitter；浮点容差只用于数值相等判断且预先冻结，不能按方向调阈值。
- signed-rank检验针对配对差值的秩位置，不标注为“均值差的t检验”。

每scope主BH family固定为 **139×3=417**。不足计算的项保持NA，但planned n不因缺失缩小。所有tiers共用此family。

### 6.3 三圈之间的星号

固定三个环对：R1−R2、R1−R3、R2−R3。使用同一单位的e作配对差；单位与section集合一致。单独双侧signed-rank，单独BH family=417。原始比例的环对差值列入描述，不与excess环对混标。

星号规则对所有图统一：
```text
q<0.0001 ****; q<0.001 ***; q<0.01 **; q<0.05 *; 其余ns; 不可计算N/A
```

图上必须写明使用哪个q字段。星号不从v6 A图复制，不因前一轮positive/depletion tier改为单侧检验。保存P、q、正/负/零单位数、效应median/mean/IQR和方法。

除每scope主q外，导出“仅七dataset合并”planned n=7×417=2,919，以及“总体+七dataset”planned n=8×417=3,336的辅助q。它们仍是当前筛选后的探索性校正，不宣称消除了聚类后选择效应。

不要求三个环都显著；不要求向0单调收敛；q≥0.05不等于已与背景等价。

---

## 7. M04：与连续 MP10 的关系，两种范围、两层评分

### 7.1 只读取已有分数

使用v4冻结的：
```text
spot_MP10_model_universe
epi_MP10_model_universe
```

同片相同模型基因背景、同一中心的两层比较。为两层共同定义 S_score：S_s中完整Top20、两层均有效的anchors；不分别挑选中心。B5缺HMGCS1的220个旧anchors只影响此评分模块，不能从组成模块删除。

不重算Top20、不换signature、不拟合c2l、不推导“整体score×或÷Epi比例”。Epi层一直标为“模型归属Epi”，不能称实测纯上皮。v5归属校准限制写入报告。

### 7.2 两种范围与去重

- `ALL_ELIGIBLE_EPI`：每片S_score全部中心，含niche0。
- `WITHIN_SELECTED_NICHE`：同片S_score中该获批niche的中心。

同一个celltype被多个niche批准时，全体范围只计算一次。全体共有29个不同celltypes，群内有139个不同niche–celltype组合。不要把相同全体结果复制成独立检验，也不从每个niche使用不同的全体背景。

注意：所有上皮中心的“相对全体平均组成差”平均为0是恒等式，没有新增检验意义。ALL范围做的是MP10–组成关联，不是检验ALL均值是否偏离其自身背景。

### 7.3 切片内关联的固定实现

每片每个范围、每圈、每层评分输出：
1. 原始 Spearman(MP10_center, q_ring_celltype)。
2. **主调整关联系数**：average-tie rank后，分别将rank(MP10)、rank(q)对 `[intercept, rank(center_mean_Epi), rank(log1p(center_nCount))]` 残差化，再算两个残差的Pearson相关。
3. **中央同类型成分敏感性**：若已有合法的中央mean celltype比例，再额外加入其rank，单列结果。不能拿环内组成、niche标签或MP10本身当协变量消除目标关系。

先阅读v5连续诊断代码，将它与以上主/敏感性设计的差异列明。若v5包含“中央对应类型”调整，保留为第3项，不把它偷偷作为所有类型的同一主检验。不得根据本轮显著性选择协变量形式。

每片每范围至少30个完整有效中心，score和目标组成各至少3个不同数值；恒定、残差方差为0、有效自由度不足等记NA。常量协变量可在记录后移除，多重共线用QR审计，不静默换模型。禁止使用spot-level P/SE作为主星号。

同片减去常数背景不会改变Spearman：可用此性质做一致性测试，不把q与q−b两种关联重复当成两条证据。

### 7.4 单位归并与星号

同一单位多片时，各片相关经固定Fisher-z等权平均后变回rho；用于数值转换的相关限制在(-1,1)，保存原值和截断记录。两个评分层使用相同合格section集合，避免层间差异来自样本变更。

各scope内对单位相关系数作双侧signed-rank against 0，规则与M03一致；主调整关联family：
```text
(29个全体celltypes + 139个群内组合) × 3圈 × 2评分层 = 1,008
```

raw关联、中央同类型调整各有独立辅助family，不用其中更小的q替换主q。七dataset辅助planned n=7,056；含总体的八scope辅助planned n=8,064。无法计算的计划项保留NA及原因，不能降低planned n获得星号。

每个范围×celltype（群内还含niche）×score_layer画R1/R2/R3三个并列关联系数箱体，0线为无关联；主图显示调整关联，raw伴随图/表保留。单位点及N必须与检验一致。不匹配MP10-low、不二分MP10、不画趋势曲线。

---

## 8. M05：有限敏感性和组成闭合解释

### 8.1 必做结构性诊断，不另造niche0机制

按同样三圈中心集合、每片/单位/niche，输出：
```text
Epi_ring = q[r,Epithelial]
TME_ring = sum(q[r,nonEpi])
```
校验每个anchor/每圈二者之和约为1。用三个并列箱体或表呈现niche0及其余niches的Epi/TME总份额，不对niche0的29个非上皮降低再次包装29条新假说。跨类型分别取median后未必相加为1，不把中位数的非可加性误报数据错误。

### 8.2 一个固定的TME内部组成敏感性

只作为分母解释，不替代主q，不重聚类：
```text
q_TME[i,r,c] = q[i,r,c] / sum(q[i,r,nonEpi])  # c为非上皮
```
先汇总环，再在环层面去除Epi；不是先对每个邻居做非上皮归一化后平均。只在三圈TME份额都≥0.05的中心上计算，阈值固定、不按结果调。

在这份新敏感性中心集合上，重新计算同片背景并同时重算raw q主口径的对照，保证“原口径→TME内部”差异不只是中心集合变化。目标/rest和MP10的覆盖规则保留。

数值结果覆盖所有139组合和两种MP10范围，主要导出effect、方向、覆盖；若报告检验，另立family，不替换主星号。不得用两个pooling总体均值相除后宣称新的样本级显著性，也不把TME内部比例视为真实绝对细胞数量。

### 8.3 另外两类预定敏感性

- 完整6/12/18三环：同标签、同算法重算，不改变主集合；不完整环影响多少明确报告。
- 已核实患者单位子集：保留原身份，不把匿名样本伪装成患者；使用相同planned family、单独q与覆盖。

连续关联另外做一个有限的空间重复抽样诊断：优先复用已保存且有物理坐标依据的v5/v3块定义；若不存在，固定为每片原array坐标8×8块，每块按seed123稳定抽取最多一个合格anchor，三圈、两层和各目标用一致抽取。只检查效应方向/覆盖，空间块不作为患者，不以块数增加样本量；样本不足如实NA。这不是完整空间置换检验，不宣称彻底解决空间相关或归属偏差。

不扩展成多种块大小、阈值或评分方法的大规模寻优。主关系仍靠独立单位层面、分dataset的可重复方向与效应解释。

---

## 9. M06：报告必须把三类证据分开

每个组合保留前一轮 `review_tier`、原始效果量、v6 q和来源警告，再加入：
1. **环境组成**：目标类型是否高于/低于同片背景，三个环各是多少。
2. **距离尺度差异**：三个环是否不同；是否只在近圈、广泛分布或近远异向。
3. **中央MP10关系**：全体及群内连续关联、两层评分、各dataset方向。

不能把v6热图上的富集、M03的环境差异、M04的MP10关联混写成同一件事。

允许的解释标签包括：
```text
COMPOSITION_CHARACTERISTIC
EXCESS_ABOVE_SAMPLE_BACKGROUND
DEPLETION_BELOW_SAMPLE_BACKGROUND
MP10_POSITIVE_ASSOCIATION_EXPLORATORY
MP10_NEGATIVE_ASSOCIATION_EXPLORATORY
TWO_LAYER_DIRECTION_CONCORDANT
ONE_LAYER_STATISTICAL_SUPPORT_ONLY
OPPOSITE_SCORE_LAYER_DIRECTIONS
TME_DENOMINATOR_SENSITIVE
DATASET_HETEROGENEOUS
INSUFFICIENT_COVERAGE
```

“两个分数点估计同号”和“两层都通过主q”分开记录；一层ns不能自动当作方向相反。既往tier不强迫本轮方向。关系小、区间跨零时不要用标签掩盖不确定性。

病例数、空间块数、anchors数分别列；患者/匿名单位分别列。多环同向不是三个独立验证，R3也不是独立队列。全体有关、群内不明显可能是环境间关系，不能仅因群内不显著就全部否定。群内属于按组成筛选后的条件分析，保持选择偏差/范围受限说明。

结果可以提出下一轮应讨论哪些组合，但不自动选“最符合预期”的dataset，不删反向结果，也不自动进入通讯或机制研究。未发现MP10关联也要完整交付环境组成结果。

---

## 10. 输出结构与文件组织

```text
neighborhood_ring123_broad_v7/
  00_approval/                 # 输入清单、receipt、哈希
  00_audit/                    # 输入版本/UID/单位/参数/软件/检验族
  01_ring123/                 # 三环缓存、coverage、轴与成员表
  02_background_composition/  # section/unit m,b,e,target-rest
  03_ring_statistics/         # 对0与环间检验；不混用q
  04_MP10_association/        # section/unit/scope原始+主调整+敏感性
  05_sensitivity/             # Epi/TME、TME内部、完整环、已核实患者、稀疏块
  06_figures/<scope>/         # raw/excess/ringpairs/MP10
  07_synthesis/              # 139组合与全部dataset综合证据
  code/
  logs/
  reports/
```

核心数值文件至少包括：
```text
coverage_all_steps.tsv
background_by_section_ring.tsv.gz
section_ring123_composition.tsv.gz
unit_ring123_composition.tsv.gz
ring123_background_tests.tsv
ring123_pairwise_tests.tsv
MP10_section_associations.tsv.gz
MP10_unit_associations.tsv.gz
MP10_scope_associations.tsv
composition_closure_diagnostics.tsv
TME_conditional_sensitivity.tsv.gz
candidate_evidence_139.tsv
candidate_dataset_evidence.tsv
family_manifest.tsv
```

所有值含原分数单位、比例或百分点单位说明，明确分析中心集合/背景/score_layer/ring/analysis_unit_policy。

绘图/回传：
- 总体及七个dataset都有结果和入口；缺失值用NA，切片少不省略。
- 所有可计算的139组合均有三圈主图，不只给显著项。全体MP10×celltype图去重一次，群内按组合输出。
- PNG≥300dpi，主要图另存PDF或SVG；透明底同时让0线、字和星号在白背景上清楚。图例、N、scope、ring和检验字段完整。
- 不需要重画全部旧A热图；引用只读旧成果即可。颜色/label固定，不重新给niche起生物学名字。
- 长类型名（含斜线）不直接作为路径。使用如 `n1_c08_excess.png` 的短文件名，并保存 `display_name_map.tsv`；标题仍显示完整原名。
- 数值覆盖完整，raw辅助图可进入每dataset的完整图包；避免将数千张图片全嵌入单个巨大文档。按预定tier而非显著性排序索引。
- `review_B.zip` 以报告、全部数值摘要、总体图和全scope缩略总览为主；各dataset完整图另存独立包。所有链接是可移植相对路径。包中不含字体、原始模型、全基因表达或旧zip的递归嵌套。

---

## 11. 必须实际实现的入口、核验与停止位置

建议实现并实际测试：
```bash
bash code/run_ring123_broad_v7.sh --stage all --approval /absolute/path/approved_gate_v7.tsv
bash code/run_ring123_broad_v7.sh --verify-only
```

阶段检查点建议：
```text
M00_APPROVAL_AND_INPUTS_VERIFIED
M01_RING123_READY
M02_BACKGROUND_READY
M03_COMPOSITION_BOXES_DONE
M04_MP10_ASSOCIATIONS_DONE
M05_SENSITIVITY_DONE
M06_DELIVERY_VERIFIED
PHASE_B_COMPLETE_AWAITING_RESULTS_REVIEW
```

验收至少包括：
- 清单139/28+32+39+40、29非上皮类型、8scopes、1112真实ID映射通过。
- 原文件/模型/标签/MP10分数无修改；输入hash与轴明确。
- 打乱输入行/类型列顺序后仍能按显式键恢复；重复UID、错列可被检测。
- R1/R2旧值复现，R3独立，中心排除，missing不当0。
- 三圈中心集合一致；背景始终同片；e与target−rest的均值恒等式通过。
- B5的组成保留、完整Top20关联排除；NA不足与全零结果区分。
- 同一真实单位多片/跨dataset不被总体重复计数。
- 主检验family为417/417/1008；未计算项不缩family；raw与敏感性不替代主q。
- 相同全体celltype关联去重，群内与全体、两层评分均明确分开。
- 图上N与实际检验N一致，箱体对0星号与环间括号星号不混。
- TME内部只是敏感性，不篡改原组成；niche0保持结构性参照。
- 抽查全部139总体excess图的存在及范围，并目视检查总体、两个主要dataset及一个稀疏scope；文字、长标签、NA、星号不遮挡。
- 从review包解压到短临时目录检查所有相对链接有效，无重复套目录。

缺少某片或某一模块数据时，记录具体missing path/字段/受影响范围，继续所有不依赖它的任务。不要为了全部成功给NA补0。磁盘或代码错误需修复实际可修问题，但不改变科学阈值；无法修复时返回真实partial状态和可运行续跑命令。

采用切片流式或少量并行、共享小型缓存。Checkpoint原子写入并绑定输入hash；复跑跳过一致的已完成阶段。先完成少量真实切片技术pilot验证后自动扩到全量，不能只完成pilot再次等待泛泛批准。

最后回传：
- `reports/RUN_SUMMARY.md`：实际完成状态、各dataset覆盖、各tier完成数、富集/耗竭与MP10结果分开。
- `reports/REVIEW_INDEX.md`：先看哪些全局/各dataset箱线图和表。
- `reports/METHODS.md`：公式、观察单位、planned families、阈值和来源限制。
- `reports/WHAT_IS_STILL_UNCERTAIN.md`：包含选择效应、组成闭合、空间依赖、未审定病理和Epi归属限制。
- `reports/validation_results.tsv`、`reports/execution_status.json`。
- `reports/review_B.zip` 及各dataset图包。

**本轮在全部已授权阶段B实际完成并交付后停止。不要再次停在旧Gate，也不要进入SpottedPy/LIANA+/通讯。**

---

## 12. 文献与本项目适配的边界

用户提供的 Khaliq et al., Nature Genetics 56, 2455–2465 (2024)，DOI `10.1038/s41588-024-01914-4`：
- PDF第4页Fig.2b图注和第13页左栏说明单侧Wilcoxon rank-sum与BH，支持前一轮A的参考式组成筛选。
- PDF第13页右栏说明箱线图median/Q1/Q3/1.5×IQR及离群点；这一外观用于本轮三个并列箱体。
- 该栏的原文箱线图检验另为two-sided unpaired t-tests；**本轮单位配对signed-rank不是该文原样统计流程。**

本轮独立环、中心范围、背景、unit聚合、139项选择、TME门槛和主检验族是本项目已讨论的适配/工作设置，不称通用文献标准。无需为本轮再开展新一轮文献检索。

以下附表是授权清单的可读版本；准确ID和scope由随附TSV与原v6表核对。

| niche | 核心 POSITIVE_CORE | 宽口径 POSITIVE_BROAD | EXPLORATORY | DEPLETION_CONTROL |
|---|---|---|---|---|
| niche1 | GC B；cDC1 | CD4 Naive/Memory-like；Lymphatic EC；Mast；pDC | CD8 GZMK/Effector-memory-like；LAMP3 DC；NK；Naive/Memory B；Tfh/CXCL13-like；Treg；cDC2 | Arterial EC；C1QC TAM；Matrix_CAF；Plasma；SMC；iCAF；vCAF |
| niche2 | C1QC TAM；CD8 Exhausted；CD8 GZMK/Effector-memory-like；NK；Naive/Memory B；Plasma；Tfh/CXCL13-like；Treg | Arterial EC；CD4 Helper/Conventional；CD4 Naive/Memory-like；LAMP3 DC；Venous EC；cDC2；iCAF；rCAF；vCAF | Capillary EC；Classical Mono；Matrix_CAF；cDC1；pDC | GC B；SPP1 TAM |
| niche3 | Capillary EC；Matrix_CAF；Pericyte | Arterial EC；C1QC TAM；CD4 Helper/Conventional；CD8 Exhausted；CD8 GZMK/Effector-memory-like；NK；Tfh/CXCL13-like；Treg；cDC2；rCAF；vCAF | LAMP3 DC；Lymphatic EC；SPP1 TAM；iCAF | Neutrophil；pDC |
| niche4 | Lymphatic EC；vCAF | — | Arterial EC；Capillary EC；LAMP3 DC；Mast；Venous EC；iCAF；pDC | Classical Mono；Matrix_CAF；SPP1 TAM；cDC2；rCAF |
| niche5 | Matrix_CAF；SMC；iCAF；rCAF | Arterial EC；Venous EC | C1QC TAM；Lymphatic EC；Pericyte；vCAF | LAMP3 DC；Neutrophil；SPP1 TAM；pDC |
| niche6 | Classical Mono；Neutrophil；SPP1 TAM | LAMP3 DC；Pericyte；cDC1；cDC2 | C1QC TAM；Lymphatic EC；Matrix_CAF；NK；Treg；iCAF | CD8 Exhausted；CD8 GZMK/Effector-memory-like；Naive/Memory B；Plasma；Tfh/CXCL13-like；vCAF |
| niche7 | Classical Mono；Neutrophil；SPP1 TAM | LAMP3 DC | CD8 Exhausted；NK；Treg | C1QC TAM；Capillary EC；GC B；Matrix_CAF；Pericyte；Plasma；SMC；cDC1 |
| niche8 | Capillary EC；Pericyte；vCAF | Arterial EC | C1QC TAM；Matrix_CAF；Naive/Memory B | LAMP3 DC；Neutrophil；SMC；cDC1；pDC；rCAF |

附表不是新增生物学结论。若只提供本MD，也可将附表按“niche+celltype”展开，再从本地v6 Gate取真实ID生成同等清单；必须校验139项及各tier数量。不得从描述性文字猜测其他39项或40项。
