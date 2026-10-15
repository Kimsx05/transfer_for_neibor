# Codex 执行规格 v8：MP10 热点与 niche／细胞富集区域的重叠和距离

## 0. 授权、科学问题与停止位置

请完整阅读后实际执行。用户已认可用 SpottedPy 补充此前没有直接检验的“局部高值区域的几何关系”。本文件是交给服务器 Codex 的执行规格，不表示聊天助手已经运行分析。

**本轮问题：MP10 高值聚集区域是否相对于同切片一般上皮位置，更常落在某个 niche 内、更靠近某个 niche 或某种细胞类型的富集区域？**

本轮完成：
1. 版本、坐标、热点统计量和软件默认行为的有限技术核验。
2. 两层冻结 MP10 的热点／冷点识别。
3. 全部 niche0–niche8 的重叠、含重叠距离、区域外邻近距离。
4. 全部29种非上皮 celltype 的原始 spot 组成热点，与 MP10 热点的距离和重叠。
5. 总体及全部七个 dataset 的分析单位层结果、空间图、箱线图、完整表及预定敏感性。

本轮不做：重新打分／筛基因、表达来源恢复或校准、C-SIDE、DestVI、模型训练、Harmony／Leiden／k-means、niche 重命名、CopyKAT、全新病理审定、inner–outer correlation、重新运行139组 ring123、LIANA+、细胞通讯、预后、机制实验或新队列下载。

无需等待病理审定，不新增“代表切片需要用户手工确认”门槛。不以出现显著结果作为运行成功标准。完成技术 pilot 后，若实现通过，继续全部已授权数据；不要只返回预检、安装记录或 pilot 结果询问是否继续。

保留 v7 的非显著连续关联结论。热点分析是另一种空间问题，不能宣称它自动推翻 v7，或者修复上皮来源不确定性。

**停止状态：数值结果、全目标覆盖、所有 dataset 可浏览图包及限制报告完成后，置 `V8_COMPLETE_AWAITING_REVIEW`。此后不自动开展通讯。**

---

## 1. 项目位置、只读输入和版本保护

项目根目录：
```text
/home/data/t070721/codex_workspace/Bladder Metabolism
```
新输出目录：
```text
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_spottedpy_v8
```

只读复用：
```text
neighborhood_ring123_broad_v7/
neighborhood_niche_celltype_screen_v6/
neighborhood_mp10_calibration_v5/
neighborhood_mp10_source_v4/
neighborhood_ring_harmony_v3/
neighborhood_ring_comparison_v2/
```

GitHub 上传别名 `results/review_B` 对应的正式服务器工作目录预期是 `neighborhood_ring123_broad_v7`。不得因为目录名不同重复生成分析；先查 manifest、代码和 UID/哈希确认。如有重复嵌套，以实际内容定位。复用现有输入/日志即可，不为本任务重新拉取整个大型结果仓库。

已核实的 review_B 归档 commit：
```text
Kimsx05/transfer_for_neibor
bbd226fc09934a5496513a71edadc7dbf4450dbc
results/review_B/
```

优先读取：
```text
neighborhood_ring123_broad_v7/reports/METHODS.md
neighborhood_ring123_broad_v7/reports/RUN_SUMMARY.md
neighborhood_ring123_broad_v7/00_audit/input_manifest.tsv
neighborhood_ring123_broad_v7/00_audit/spec_frozen.json
neighborhood_ring123_broad_v7/01_ring123/anchor_metadata.parquet
neighborhood_ring123_broad_v7/01_ring123/graphs/<dataset>/<section>__node_axis.parquet
neighborhood_ring123_broad_v7/01_ring123/graphs/<dataset>/<section>__R1.npz
neighborhood_ring_comparison_v2/03_features/anchor_master.parquet
neighborhood_ring_comparison_v2/03_features/feature_columns.json
neighborhood_ring_comparison_v2/00_preflight/coverage63.tsv
neighborhood_ring_harmony_v3/03_main/SAMPLE_R12_LEIDEN/labels.parquet
neighborhood_ring_harmony_v3/02_R12_pilot/SAMPLE_WEAK_R12/LEIDEN/labels.parquet
neighborhood_mp10_source_v4/04_niche_comparison/anchor_MP10_all_layers.parquet
neighborhood_mp10_calibration_v5/00_audit_and_spec/identity_and_independence.tsv
neighborhood_niche_celltype_screen_v6/00_audit/analysis_unit_policy.tsv
```
q05 全组织路径和中央 mean abundance 路径从 `coverage63.tsv` 的 `q05_path/mean_path` 读取，不凭文件名猜列序。图像、坐标和比例尺优先复用 v3 及原 Seurat 图像提取代码。

只在项目根目录、既定数据目录和当前上传目录有界查找；不要全盘扫描，不把 `/mnt/data` 视作服务器目录。未找到原文件时，明确列缺失项，能从受信任的冻结 manifest 重建的仅重建轻量元数据，不重训模型。

记录输入绝对路径、大小、SHA256、读取用途和 UID 数量。不要覆盖旧输出，不终止其他作业，不删除旧结果。先查看磁盘、配额和已有进程。

---

## 2. 方法来源与适配：必须如实区分

### 2.1 原方法

Withnell E, Secrier M. *SpottedPy quantifies relationships between spatial transcriptomic hotspots and uncovers environmental cues of epithelial-mesenchymal plasticity in breast cancer*. Genome Biology 25, 289 (2024). DOI: `10.1186/s13059-024-03428-y`。

原文支持热点／冷点、多切片距离分析、一般肿瘤位置背景、热点尺度敏感性。它不是我们的上皮 MP10 验证论文，也没有规定本项目的 gate、统计单位或检验家族。

公开代码锁定：
```text
https://github.com/secrierlab/SpottedPy
commit: 9dc364f203b7617b6701aa8fb579cbbb35357b26
spottedpy/hotspot_helper.py
  blob: c34129f78e91ba9a9ebdef4344813079913d3ea4
spottedpy/main.py
  blob: 8d7c8078c853ddbbed45505a503545a36279a3ac
```
该源码调用 `G_Local` 时没有显式设置 `star=True`，热点依据 `Zs` 的方向及原生 `p_sim<0.05`；局部 FDR 代码未启用。距离按 `array_row/array_col` 的欧氏最短距离计算，目标是相应列非缺失的位置；默认可对缺失目标填最大距离。

### 2.2 本轮固定选择

为避免在看到结果后选择 G 或 G*：
- **正式主统计量固定为公开源码路线的局部 G：`star=False, transform='R'`。** 名称使用 `SOURCE_G`，不要误称为 G*。
- 论文写 G*，故在技术 pilot 的两个 MP10 层上额外运行 `star=True`，记录局部统计量、热点比例和 Jaccard 的区别。仅作为实现审计，不据阳性多少更换正式主统计量，也不扩成另一套全项目主分析。
- 使用真实、等距度量正确的组织坐标，而不是直接把未经几何变换的 Visium array_row/array_col 当欧氏坐标。
- 连通热点由**原完整组织一阶图的诱导子图**定义，不在热点子集重新找近邻，把远处热点强接在一起。
- 缺失目标不填最大距离；不同组织物理连通分量之间不计算主距离。
- 跨样本关系使用预先定义的分析单位差值推断，不使用绘图函数默认的逐 spot t-test／GEE 输出替代。

以上是显式项目适配。允许写一个调用 PySAL `G_Local`、复用 SpottedPy 热点逻辑的轻量 adapter；必须与锁定源码在相同输入、相同坐标、相同权重下做统计量一致性核验。不要把自定义了坐标和连通规则的版本称为逐项原样复现。

安装版本的 `esda.G_Local` 默认值及 `p_sim` 尾概率约定必须读取实际签名和源码确认。新版可能出现 `alternative` 参数，不能只凭旧版本印象。参考：
```text
https://pysal.org/esda/v2.9.0/generated/esda.G_Local.html
https://pysal.org/esda/stable/generated/esda.G_Local.html
```

Khaliq et al. 2024（DOI `10.1038/s41588-024-01914-4`）的组成热图已在 v6 完成。本轮不重做该热图，也不声称其单侧 rank-sum 就是本轮热点距离检验。箱线图外观继续使用 median、Q1/Q3、1.5×IQR 规范。

---

## 3. M00：冻结来源池、目标池与不变量

### 3.1 主 MP10 来源池 E0

优先直接复用 v7 `score_eligible=True` 的 UID：
- 当前正式 accepted niche 标签0–8；niche0 保留。
- 原 mean Epithelial/all30 gate≥0.25 不变。
- 原三圈共同覆盖与完整 Top20 条件不变，以便与 v7 对齐。
- `spot_MP10_model_universe` 和 `epi_MP10_model_universe` 均有限，`epi_score_valid=True`。

既有核验规模：63张目标切片，61张有主来源点；v7主组成中心120,630，主评分中心120,439。核验实际输入，不硬填这些数。B5可保留组成数据，但不回填不完整Top20分数；两张无共同中心的切片仍保留状态行。

两层评分用相同E0，热点标签可以不同。禁止为了让热点更清楚而改Top20、排名背景、归属结果、标准化/平滑分数或 Epi 门槛。

本轮不默认启用“扩展边缘 anchors”。R3不是热点方法的必需条件，但现在冻结E0是为了保持比较基础一致。未来如需扩展，只能单列来源池敏感性，不能替换本轮主结果。

### 3.2 组织物理池 T

保留完整 tissue=1 的物理图及其真实连通分量，包含没有组成输出的物理节点。这些节点可用于判断物理连通，但不为其生成MP10或细胞比例；不填0。

### 3.3 niche目标池 N_k

每片使用冻结 v3 的所有 accepted common anchors 中 `accepted_label=k` 的位置作为目标，k=0…8。目标不要求自己进入v7三圈评分池，防止把原本存在的niche边缘人为删去。

重叠比例的计算自然在E0内取交集；距离可到同片同组织连通分量内的全部目标成员。明确输出“全niche目标数”和“来源池内niche数”。不把无标签间质spots涂成某个niche。

### 3.4 细胞热点目标池 C

使用整片组成有效的 tissue spots，而不只是E0或上皮spots：
```text
p[j,c] = q05_abundance[j,c] / sum_over_all30(q05_abundance[j,*])
```
列名、列序与冻结字典完全一致，Epithelial以外恰为29类。保留SMC、Pericyte、各CAF亚型的区别。

不得输入q1/q2/q3环均值、Harmony坐标、dominant label、聚类centroid或TME内重归一化比例来生成主细胞热点。目标称“该类型相对组成热点”，不是已知绝对细胞数量热点。

### 3.5 稀疏、常量与缺失

常量、全零、可用位置不足或无法生成热点的变量记录原因；不加噪声、不补零、不改阈值。热点识别以有限有效值为输入；缺测值保持缺测。

使用少量明确门槛：全切片参与该变量热点计算的位置至少30、至少3个不同数值；KNN所在真实组织连通分量至少k+1个该变量有效位置。其余位置标 `LOCAL_SUPPORT_INSUFFICIENT`。两层MP10共享几何支持域，不因为某层分数高低改变几何域。

---

## 4. M01：坐标、物理连通和KNN几何

### 4.1 正确的度量坐标

优先从已验证全分辨率图像像素坐标获得XY，使用完整组织一阶真实邻接边估计该片的中位spot间距 `pitch_px`：
```text
xy_pitch = xy_fullres_px / pitch_px
```
主距离单位为 `spot pitches`，1表示一阶spot中心间距。仅在有可靠物理标定时额外输出μm；不能把像素或array坐标直接标成μm。

若仅有array格点，先验证六角偏移编码，再转换成欧氏等距坐标。标准 `(Δrow,Δcol)=(0,±2),(±1,±1)` 编码可用 `x=col/2, y=sqrt(3)*row/2`；不是所有数据都先验使用这一编码，必须验证。

保存原坐标、转换后坐标、尺度来源、六个一阶邻居距离一致性、图像叠加核验。图像旋转/镜像只用于复用原显示方向，不为对齐不同患者做几何注册。

### 4.2 主热点权重

对每片、每个输入域、每个k：在同一个**真实组织连通分量内**按欧氏距离取k个近邻，不含self；使用二值有向KNN输入，交由G_Local做行标准化。不同组织块不相连，源域和目标域分别构图。

多个有效组织分量可组成block-diagonal权重，在该片的有效变量值背景上计算局部G；不跨切片、不以dataset为局部G背景。

先按固定canonical spot_uid顺序排列每个输入域，W、值向量和输出均沿用该轴；冻结稳定的距离tie规则（距离、spot_uid），不得让坐标排列顺序改变邻居身份或置换抽样的对应关系。相同域的两层MP10共享权重；全组织29类型有效行一致时共享权重。仅构图一次并缓存。

### 4.3 必须保留真实空间范围审计

k不是物理环数，上皮筛选后的kNN可能较远。逐位置保存：
- 到第k近邻距离（pitch及可用μm）；
- 邻居是否同一组织连通分量；
- 近邻中有多少落在真实R≤3；
- 是否存在较长边或绕过孔洞的风险。

主结果不把k10叫“10个固定物理邻居”或“第二圈”。增加固定 `GEOMETRY_LOCAL_SUPPORT` 描述性敏感性：来源中心及其k近邻均在真实R≤3、且第k邻居距离≤3pitch时保留；在这个来源子集重算相同目标下的H和背景效应。此敏感性不重新选k、不把缺少局部支持的位置冒称为回避。

热点组件使用完整组织R1图在hot成员上的诱导子图。不要用热点点集的kNN代替原组织图；不得跨未测空隙合并团块。

---

## 5. M02：有限技术pilot，然后实际全量计算热点

### 5.1 选择pilot

每个有数据的dataset选择1张切片，共最多7张。选择规则为满足输入完整性的切片中，E0大小最接近该dataset中位数；并列按section_id。不得依据MP10效果或某celltype热点多少选pilot。

用pilot核对：
1. SOURCE_G adapter与锁定源码/直接PySAL调用在同一W下的G、Z和p值。
2. star=False与star=True在两层MP10上的差异，只记录，不换主方法。
3. 相同seed重复运行、置换/UID对齐、坐标距离正确性。
4. 几个人工小图：无目标、全零、单点、两个不相连组织块、共享来源/目标点距离0、已知平移后的距离变化。
5. 保留模拟矩阵与否的p值一致性；确认节省内存不会让p或Z缺失。
6. 每任务耗时、内存和预计全量时长；记录环境与依赖，而非为了加速降低正式置换数。

技术通过即继续全量；不因pilot没有热点而判技术失败。实现/身份不明只阻塞相应切片，继续其余独立模块。若两层/两类方法比较无法保持定义一致，明确停止该比较，不自创替代分数。

### 5.2 固定热点参数

```text
statistic = SOURCE_G
star = False
transform = R
primary_k = 10
sensitivity_k = [6, 18]
permutations = 9999
seed = 123 (必要的稳定派生种子记录)
relative_background = section
component_min_spots = 5
min_variable_positions = 30
min_distinct_values = 3
```

对两层MP10以及29种非上皮spot组成实际执行三个k，共31个变量域/每片；相同目标热点不要因来源层不同重复运行。

### 5.3 两个热点阈值轨道

保留两个**预先固定、绝不互相替换**的结果轨道：

**A. `SOURCE_P05`（原代码风格探索轨道）**
- 保存原生 `p_sim` 和其明确的尾概率约定。
- 按锁定源码 `Zs>0 & p_sim_native<0.05` 标原始hot，`Zs<0 & p_sim_native<0.05` 标原始cold。
- 这是未做局部多重校正的探索图。即使后面单位层的关系检验过BH，也不能把其组成点叫“FDR热点”。

**B. `LOCAL_BH05`（较严格的热点敏感性轨道）**
- 主程序需提供明确的双侧局部置换P：若所用原生p_sim是官方旧实现的小尾/单侧约定，使用 `min(1, 2*p_sim_native)`；若实际库已返回双侧P，直接使用，不能再次乘2。
- 必须以实际源码和小规模保留置换样本核验尾概率，保存 `p_sim_native, p_two_sided, pvalue_convention`。
- 每片、每变量、每k，对全部局部可检验位置做一次BH，然后按Z方向分hot/cold。不能先挑正Z再缩小BH家族。
- LOCAL_BH05使用 `q_local<0.05`；分别在该轨道上计算来源MP10热点和细胞目标热点，不能把一种轨道的来源混接另一种轨道的目标。

两轨道共同保存原始点级hot/cold及真实R1组件ID。正式区域摘要使用≥5个连通成员的组件；1–4点组件不删除原始记录，单列小区域清单。组件数量在过滤后重新计算，不直接信任上游过滤前n_components。

任一轨道没有hot区域，保留 `NO_SOURCE_HOTSPOT` 或 `NO_TARGET_HOTSPOT`，不放宽P、改k或把top分位点补成热点。LOCAL_BH05没有结果，不自动用SOURCE_P05替代它。

### 5.4 热点质量信息

逐片／变量／k／轨道记录：有效位置、hot/cold点数、组件数、各组件大小、原始值中位数/IQR、来源位置自身MP10分位数、热点比例、KNN空间范围、局部P最小可分辨率、孤立点与缺失原因。

SOURCE_G不包含中心自身，hot标签表示局部周围高值，不保证每个成员中心自己都高分。必须同时展示成员自身MP10分布，不暗中额外加top20%门槛。

细胞热点附原始组成值和全组织背景。若热点只是很低比例上的相对高值，标 `LOW_ABSOLUTE_RELATIVE_COMPOSITION` 描述警告，不把它称为大量该类细胞积聚，也不自动删除稀有细胞。

---

## 6. M03：九个niche的空间重叠与距离（主问题A）

以下按每片、评分层、k、热点轨道分别计算。E为源池中具有有效局部统计支持的位置；H为E中属于正式hot组件的位置。N_k为冻结niche目标集。

### 6.1 O：相对同片上皮背景的重叠富集

```text
overlap_hot = |H ∩ N_k| / |H|
overlap_background = |E ∩ N_k| / |E|
DeltaO = overlap_hot - overlap_background
DeltaO_pp = 100 * DeltaO
```

输出分子分母、hot占E比例、niche目标的全组织数/源池内数、rest数。正值=在热点中占比更高；负值=占比更低。不把这个量直接叫OR或风险比。

同片没有该niche时标 `TARGET_NICHE_ABSENT`，不把必然的0差当成支持“无关系”。几乎整片都是该niche、没有足够rest时只作描述。

### 6.2 D：含重叠的距离

对每个来源位置i，计算到最近目标成员j的欧氏距离，采用xy_pitch：
```text
d(i,N_k) = min_j ||xy_pitch[i] - xy_pitch[j]||
DeltaD = mean(d over H_valid) - mean(d over E_valid)
```

主方向：DeltaD<0表示更近，>0表示更远。**距离源到目标是有方向的，不默认等于目标到源。** 不额外跑反向36/29种分支。

H_valid、E_valid只能在含目标的同一组织物理连通分量内产生可达目标。目标只在另一块独立组织中时，不跨空白计算主距离；对应分量的H和E都从该距离比较中移除，并报告移除数和可达覆盖。重叠结果仍使用它自身定义的域，二者分母不同必须明确。

主距离包含重叠产生的0；禁止只报告所有点的最小值。另保存median distance和热点组件等权摘要作为描述性敏感性，不以结果最显著的聚合替换mean。

### 6.3 D_out：区域外接近

```text
H_out = H_valid \ N_k
E_out = E_valid \ N_k
DeltaD_out = mean(d over H_out) - mean(d over E_out)
```
对照也必须排除目标区域，不能拿H_out与完整E比较。空集或太小标NA。该量区分“位于niche内”与“从外部靠近”。

### 6.4 推断覆盖门槛（项目工作参数）

先为所有情况保存描述量。只有满足下列条件才纳入相应单位层正式检验：
- O：E≥30，H≥20，E\H≥20；E内目标≥20、非目标≥20。
- D：实际可达E≥30，H_valid≥20，E_valid\H_valid≥20；同片有效目标至少5个位置。
- D_out：实际E_out≥30，H_out≥20，E_out\H_out≥20；目标至少5个。

允许组件最小5但切片正式hot总数不足20时只作描述。不要降低阈值凑N。所有门槛写入frozen spec，不能看到星号后再变动。

### 6.5 包含自身的背景是明确选择，不是假造独立对照

E包含H，符合用户“样本自身一般上皮环境”要求。保存nonhot背景的补充差值，并验证：
```text
Delta = (1 - |H|/|E|) * (mean_H - mean_nonhot)
```
对O将每点的y设为niche成员指示；对D使用相同可达域。该恒等式用于审计与解释效应收缩，不把重叠的两个集合当独立样本做rank-sum/t-test。

---

## 7. M04：29类细胞组成热点的空间关系（主问题B）

对每个celltype c，将相应阈值轨道、相应k的正式组成热点T_c作为目标：
- 来源仍为两层MP10热点，绝不重新定义源分数。
- T_c从完整组成有效组织得到，不局限于上皮或特定niche。
- `DeltaD`按第6节可达域和相同背景计算，是本模块主终点。
- `DeltaO`和`DeltaD_out`作为解释性辅助，分开输出与校正。
- 没有目标热点=NA，绝不能填最大距离；T_c存在但不在源所在组织连通分量=不可达，并报告覆盖。

主B回答MP10区域与细胞富集区域的关系，不是某个单细胞发生接触。不把niche7等价于Neutrophil热点；两个目标独立计算，后面交叉解释。

禁止调用距离函数的危险默认值。若调用SpottedPy API，必须显式：
```text
empty_hotspot_default_to_max_distance = False
```
变量成员按非缺失值识别时，非成员填NaN，而非0或False。距离阵列主计算用cKDTree或分块最近点查询，不建立全组织N×N稠密矩阵。用小样本与SpottedPy `calculateDistances` 比较验证结果。

---

## 8. M05：分析单位统计及BH家族

### 8.1 先切片，再分析单位

每片分别计算效应，不混合跨dataset原始坐标或分数。按既有 `analysis_unit_uid` 对可用切片的效应等权平均。不能按anchor数或热点数让大切片获得更多生物学重复。

同一真实单位跨dataset出现，总体只算一个单位；dataset结果仍展示各自切片，标记共享来源。已核实患者、匿名独立样本、其他既有政策单位分开计数，不全部重命名为患者。

主箱线图和星号使用完全相同的单位行。辅助提供全部描述切片、组件级原始表，但组件或spots不作为独立患者。

### 8.2 检验

对单位效应vs0做双侧Wilcoxon signed-rank；不是两个重叠spot集合的rank-sum。延续v7：至少5个单位及5个非零差值才正式检验；至少5个单位且全部差值为0时P=1；其余NA附理由。无零/绝对ties且n<50时exact，否则显式记录tie/continuity处理。

统计检验针对差值分布位置，不宣称专门检验均值；同时输出mean/median/IQR、正负单位数、N及效应点。近零浮点误差按固定数值容差处理并审计，不能加jitter打破ties。

保留源hot比例过大、只有一个大团块、target少等警告。单位层P表达的是同片背景差是否跨单位一致，不是已经完成保留空间自相关的严格随机位置零模型。

### 8.3 预先固定的检验家族

在正式 `SOURCE_G, k10` 中，每个scope、每个热点轨道各自：

| family | 成员 | planned n |
|---|---|---:|
| NICHE_PRIMARY | 9 niche × 2评分层 × 2终点(O,D) | 36 |
| CELL_PRIMARY_DISTANCE | 29细胞类型 × 2评分层 × D | 58 |
| NICHE_OUTSIDE_SECONDARY | 9 niche × 2评分层 × D_out | 18 |
| CELL_OVERLAP_SECONDARY | 29细胞类型 × 2评分层 × O | 58 |
| CELL_OUTSIDE_SECONDARY | 29细胞类型 × 2评分层 × D_out | 58 |

全部可计算P使用BH，但保留原planned n，不因缺失目标或单位不足缩小候选家族。主图星号基于各自q_scope；另给七dataset联合q以及包括总体的八scope联合q，明确不是独立复制队列。

SOURCE_P05和LOCAL_BH05的关系结果分表，不从两者中择优。SOURCE_P05过单位BH不意味着其所有点都过局部FDR；LOCAL_BH05热图没有hot时显示NA。

k6、k18的效应与覆盖作为预定尺度敏感性，默认不提供可与主k10竞争的“最佳P”。star=True仅技术pilot，不参与自动候选筛选。已核实患者子集另表同家族，不能用更好的子集替换全体。

### 8.4 两评分层的公平比较

两层共享E0和几何源域，但热点成员和是否存在热点可以不同。单层结果各自完整保留；“两层一致”必须另外在同一可用切片/单位交集上重算并报告，不把不同N的正负号直接称为配对复现。

不把数值接近0的异号自动标成机制冲突。保存实际差值、效应区间或IQR及支持数。上皮归属层仍保留v5的来源恢复限制。

---

## 9. M06：有限且固定的敏感性，不扩成方法竞赛

必须完成的敏感性：

1. **k6/k18**：对两层MP10及29目标都完成热点，采用匹配k的目标和来源计算同样效应；显示覆盖变化与配对同单位的方向，不按最优k选结论。
2. **LOCAL_BH05**：同一统计量、同一k10的严格热点版本；全量保留无hot的NA。
3. **几何局部支持子集**：按第4节固定标准检查长KNN边对主结果的影响，不重新调半径。
4. **已核实患者子集**：使用原患者身份，不把匿名样本改名。
5. **组件等权描述**：与主“hot spots等权”的切片摘要并列，判断是否由一个巨大组件主导；不可把组件当独立n。
6. **技术条件分层背景**：无需匹配MP10-low；固定如下定义，对所有主k10目标统一运行，而非只给显著候选临时调参。

### 技术条件分层背景

在每片有效E上，不看niche/目标，按中央Epi比例和log1p(nCount)各切成三分位，形成最多3×3层。ties合并边界；若任一非空层少于20个E点，按固定顺序整体退到2×2；仍不足则标 `STRATIFIED_BACKGROUND_UNAVAILABLE`，不要改成单层后声称已控制。

用hot的层权重标准化背景：
```text
B_strat(y) = sum_b (n_hot,b / n_hot) * mean(y over eligible background in stratum b)
Delta_strat = mean(y over hot) - B_strat(y)
```
其中y为目标成员指示或最近距离，使用对应metric实际有效域；outside版本也在outside域分层。保存热区/背景的Epi和深度分布。该敏感性仍是描述性调整，不证明消除了所有混杂，不生成匹配low组。

本轮不启动高成本空间平移、Moran谱随机化、全套点过程拟合或新的归属校准。若有后续希望升级为严格空间排斥结论的候选，在限制文件中注明尚需相应验证，不把本轮有限检查说成已经完成。

---

## 10. M07：出图、dataset组织和说明

所有七个dataset均输出；无可估计结果仍提供状态页。niche按0–8固定顺序，celltype按冻结30类中剔除Epi后的顺序。显著性只控制标记，不决定是否出现在表中。

### 必需图形

A. 每切片MP10空间页面：
- 原冻结spot分数和model-Epi分数，分开色标并写清名称；
- source hot/cold/其余有效anchors，缺测位置灰色；
- 固定niche图以及hot叠加；
- k10 SOURCE_P05和LOCAL_BH05并列，注明实际N和组件数；
- 原物理tissue轮廓/H&E能复用则叠加，缺图不阻塞数值。

B. 九niche总体及dataset热图：
- overlap excess(pp)、distance excess(pitch)和outside distance独立图例；
- 每评分层分开，不能将距离负值与重叠正值都不解释地涂成“上调”；
- 数值、q、有效单位数在表中对应，NA不是白色0效应。

C. 29细胞目标距离热图及支持表：
- 主DeltaD之外，附目标hot数量、原始目标相对组成值与可达覆盖；
- 相同target图只生成一次供两层来源复用。

D. 单位层箱线图：
- y为单位DeltaO或DeltaD，0为同片一般上皮背景；
- 箱体median/Q1/Q3，须1.5IQR；每个单位点画一次，不重复画fliers；
- 已核实患者与其他单位形状区分，透明底PNG与矢量PDF，清楚图例；
- 不画ring距离曲线。k敏感性用并列箱体/矩阵，不连接点伪造单调趋势。

E. 候选空间叠加图：
- 全9 niche的来源图均可查看；
- 全29细胞目标热点图均可查看；
- 可额外优先展示已有关注的GC B、SPP1 TAM、Neutrophil、Matrix_CAF及血管相关类型，但不省略其他结果。

连续MP10地图沿用已定plasma/0黑的风格；不要为hotspot检测使用展示截断/显示z-score。niche颜色沿用v3，hot/cold需色盲友好并有图例。所有图保持真实比例和原图方向，不能把UMAP用于几何叠加。

HTML/Markdown索引全部使用包内相对路径；图上的标星可用 * / ** / *** / **** 或数值q，但必须说明是q_scope关系检验、local热点P还是localBH，禁止混用。

---

## 11. M08：完整证据表与候选分级

必须交付：

```text
00_audit/input_manifest.tsv
00_audit/spec_frozen.json
00_audit/source_versions.json
00_audit/source_vs_paper_parity.tsv
00_audit/analysis_unit_policy.tsv
01_geometry/coordinate_and_graph_audit.tsv
01_geometry/knn_extent_summary.tsv
02_hotspots/hotspot_coverage_all.tsv
02_hotspots/hotspot_components.tsv.gz
02_hotspots/point_hotspot_statistics.parquet
03_niche_geometry/niche_section_effects.tsv
03_niche_geometry/niche_unit_effects.tsv
03_niche_geometry/niche_scope_tests.tsv
04_cell_hotspot_geometry/cell_section_effects.tsv
04_cell_hotspot_geometry/cell_unit_effects.tsv
04_cell_hotspot_geometry/cell_scope_tests.tsv
05_sensitivity/k_sensitivity.tsv
05_sensitivity/hotspot_threshold_sensitivity.tsv
05_sensitivity/two_score_layers_matched_units.tsv
05_sensitivity/geometry_local_support.tsv
05_sensitivity/stratified_background.tsv
05_sensitivity/component_weighting.tsv
07_synthesis/niche_geometry_evidence_all9.tsv
07_synthesis/cell_geometry_evidence_all29.tsv
07_synthesis/dataset_evidence_all.tsv
07_synthesis/unavailable_results.tsv
reports/RUN_SUMMARY.md
reports/METHODS.md
reports/WHAT_IS_STILL_UNCERTAIN.md
reports/REVIEW_INDEX.md
reports/INDEX_ALL_DATASETS.html
reports/validation_results.tsv
reports/RESUME.md
```

可以按更合理的存储布局分文件，但最终索引必须映射到所有逻辑产物。避免重复生成多份巨大表达矩阵或距离矩阵。

证据表每行明确：dataset/scope、源层、statistic、k、阈值轨道、目标类别与ID、终点、分析域、效应方向和大小、P/q、N、已核实/其他单位数、正负支持、局部hot覆盖、可达性、稀有信号、两层配对一致性和各敏感性。

建议的解释标签可以包括：
- `EXPLORATORY_OVERLAP_ENRICHMENT`
- `EXPLORATORY_PROXIMITY`
- `EXPLORATORY_EXTERNAL_PROXIMITY`
- `EXPLORATORY_SPATIAL_SEPARATION`
- `SOURCE_P05_ONLY`
- `LOCAL_BH_AND_UNIT_TEST_SUPPORTED`
- `SCORE_LAYER_SENSITIVE`
- `TECHNICAL_BACKGROUND_SENSITIVE`
- `SCALE_SENSITIVE`
- `DATASET_HETEROGENEOUS`
- `NO_SOURCE_HOTSPOT` / `NO_TARGET_HOTSPOT`
- `TARGET_UNREACHABLE` / `INSUFFICIENT_UNITS`
- `NO_SUPPORTED_ASSOCIATION`

这些标签只是证据分类，不能自动写成上皮特异、恶性细胞、免疫排斥、主动招募或促癌机制。方向一致、单层q过线、双层配对q过线分别存字段。

表末保留“本次是v7之后的探索；与既有连续环结果可能不同；不以热点阳性否定既有非显著结果”。不能删除阴性dataset后称为跨队列稳定。

---

## 12. 最小技术验收

至少自动验证：
1. UID唯一、错误UID失败、随机打乱输入行后显式连接恢复同一结果。
2. 标签与两层分数只读不变；full20/有效来源与v7对应；无score位置不补0。
3. 两评分层几何支持域一致，29目标原始组成有效行与冻结字典对应。
4. 六角一阶边距离近似相等、坐标缩放前后pitch距离一致。
5. source/target/component不跨切片；物理组件来源为完整tissue图。
6. adapter/source在相同W下数值一致；实际G/G*及P尾约定写明。
7. 源code_P05与localBH标签不混；BH用全部可检验局部位置；P不无故为0。
8. 真实R1组件诱导子图：两个被空白隔开的热点不能被重新KNN接起来。
9. 无hot、无target、小target、全零、常量、目标在别的组织块各有明确NA状态。
10. 最近点距离小样本与暴力欧氏/SpottedPy一致；自重叠为0；外部距离两集合均排除目标。
11. 禁止最大距离填充、跨组件最近点、全局最小点对替代均值距离。
12. 全体背景包含hot的收缩恒等式在O和D分别通过。
13. 单位等切片归并、已知患者不重复、图N与检验N一致。
14. 36/58及辅助planned families、两阈值轨道、八scope输出齐全，NA不缩小家族。
15. 两层“配对一致性”使用相同实际单位集；不拿不同N作伪配对。
16. k10为主，不按敏感性P值重选；没有热点也是有效计算结果。
17. 索引相对链接和各dataset zip均能离线打开，抽查透明图在白底可读。
18. 全部旧输入哈希不变，无新模型/评分/通讯任务。

---

## 13. 资源预算、checkpoint与异常处理

SpottedPy所需依赖可在独立环境或项目隔离目录安装；不升级破坏已有R/Python环境。保存Python、SpottedPy commit、esda、libpysal、numpy、numba、scipy等版本及实际import路径。

9999置换×31变量×3尺度计算量不小。按section/variable/k缓存并完成以下顺序：
1. pilot和版本核验；
2. 全部两层MP10的主k10；
3. 九niche主关系；
4. 全29细胞目标主k10及其关系；
5. k6/k18预定敏感性；
6. 汇总、验证与打包。

这只是执行优先级，不是只完成前几步就声称全部完成。若会话中断，保存真实状态和续跑入口；不能用缩小候选、减少正式置换或只跑显著目标悄悄替代。

并行默认最多2个切片worker，各 `G_Local n_jobs=1`，BLAS线程1；根据实测资源可降低并行，不突破服务器限制。生产模式尽量 `keep_simulations=False`，前提是pilot证明所需统计量/P与保留模式一致；不跨全部切片保存N×9999矩阵。需要验证尾概率时仅小pilot保存或分块累计。

距离用KDTree；全量点P表用parquet/zstd；中间权重稀疏存储。checkpoint键包含输入哈希、spec哈希、源码/依赖版本、section、变量、k、统计量、种子。

如果真正出现资源不足，暂停未完成分支并报告实际完成数量、原因、剩余工作和可直接续跑命令，不把计划量写成完成量，不自动运行nohup隐藏任务。

---

## 14. 结束回传格式

按以下顺序报告：

**A. 实际完成与执行身份**
数据/分数/标签版本、主statistic、坐标单位、实际切片和变量数、pilot结果、是否所有尺度完成。

**B. 热点是否存在**
每dataset两层source hot/cold数量、有效切片、源码P与局部BH覆盖区别；没有hot与检测失败分开。

**C. 九niche几何结果**
分别列重叠、含重叠距离、区域外距离，给效应和N，不只列q值。不强迫给出“最佳niche”。

**D. 29细胞目标结果**
所有目标的方向/覆盖/统计支持，包括稀有目标、缺失和相反方向。

**E. 稳健性与限制**
两层配对一致性、k敏感性、局部BH、技术背景、几何及已核实患者子集；仍未解决的来源归属和空间依赖。

**F. 文件及续跑**
给实际存在的服务器绝对路径、轻量 `review_spottedpy_v8.zip`、七个dataset图包、完整数值表和运行命令。

实现：
```text
code/run_spottedpy_v8.sh
code/run_spottedpy_v8.sh --resume
code/run_spottedpy_v8.sh --verify-only
```
verify-only只能核验，不重跑全部置换。

轻量包包含本prompt、frozen spec、方法/状态、关键数值表、总体图与每dataset至少一张技术选定的overview、相对索引；大规模各细胞空间图放单独dataset图包。不要上传基因表达大矩阵、posterior draws或依赖环境。

最终停在 `V8_COMPLETE_AWAITING_REVIEW`；若未完成如实用 `V8_PARTIAL_*`。本轮不自动推进细胞通讯、预后或机制。
