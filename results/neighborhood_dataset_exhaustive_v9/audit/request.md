# 膀胱癌 MP10 邻域数据集穷举分析 Codex 执行 Prompt

版本：v9；设计日期：2026-10-18。

请在下述服务器项目中实现并执行本任务。用户已经授权：从当前七个数据集中穷举所有至少包含三个数据集的组合；完整重复现有邻域降维、选择 niche、环分析及 SpottedPy 流程；优先寻找含 inhouse 的、有 BH 支持且可解释的空间关系；严格限制新增磁盘占用。请直接完成已授权工作，不再询问是否运行下一阶段，也不要沿用旧任务“等待审阅后继续”的停机状态。

本轮是已有数据上的穷举探索。允许比较全部结果和选择后续研究对象；如实区分组合内显著、整轮搜索校正后显著、稳健性支持及独立验证。不能承诺必然获得阳性，也不能因已经获得阳性而提前结束穷举。

## 1. 工作目录、输入与授权边界

服务器项目根目录：

`/home/data/t070721/codex_workspace/Bladder Metabolism`

新输出目录建议为：

`/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_dataset_exhaustive_v9`

先检查该目录是否存在。若是同一配置的未完成任务，校验后续跑；若属于不同任务或已完成版本，创建独立版本子目录，不覆盖既有结果。代码、配置、轻量结果及本次任务专属临时目录均置于新输出目录内。

优先使用服务器现成输入和环境，使用 rg / rg --files 定位。读取适用的 AGENTS.md、当前 COMPLETE/manifest、实际代码及 METHODS；将历史任务内的 stop_after、awaiting_review 当作历史运行状态，不当作本轮用户的新审批要求。

定位并核对这些阶段的实际服务器路径：

- neighborhood_ring_comparison_v2：真实环、原始组成特征、共同 anchors、原始参考抽样方法。
- neighborhood_ring_harmony_v3：当前 SAMPLE_R12_LEIDEN、输入契约、Harmony/聚类/标签转移实现。
- neighborhood_mp10_source_v4、neighborhood_mp10_calibration_v5：现有两层 MP10 评分与有效性记录。
- 生成 review / review_B 的后续 niche 解释和 v7 R1/R2/R3、连续 MP10 关联任务。
- neighborhood_spottedpy_v8：现有项目 adapter、源域规则、热点、几何指标、推断和缓存。

路径可能带版本子目录，必须根据文件和 manifest 定位；不能仅凭相似名称把更早的 Wang K12 方案当作本轮基础。代码及方法的只读参照见本文件末尾链接。服务器有完整输入时，不下载 GitHub 的大 ZIP，不复制整个旧项目。

已授权本轮自动创建结果、顺序运行组合、生成必要图表，以及删除本轮任务自己新建且可再生的临时文件。旧项目、原始表达、cell2location、评分、旧标签、旧图包、旧缓存和其他任务文件全部只读。不得为释放空间删除、移动或覆盖这些旧文件，也不得终止其他人的进程。

## 2. 精确的穷举范围

冻结以下七个 dataset_id，保留实际大小写与完整标识：

1. CanDis__UNVERIFIED
2. GSE171351
3. GSE238145
4. GSE285715
5. GSE319536
6. Inhousedata__UNVERIFIED
7. NPJ__UNVERIFIED

用户对 inhouse 是优先偏好，因此本轮不把是否包含 inhouse 设为数据排除条件。

枚举全部大小为 3、4、5、6、7 的子集：

| 数据集个数 | 全部组合 | 含 inhouse | 不含 inhouse |
|---|---:|---:|---:|
| 3 | 35 | 15 | 20 |
| 4 | 35 | 20 | 15 |
| 5 | 21 | 15 | 6 |
| 6 | 7 | 6 | 1 |
| 7 | 1 | 1 | 0 |
| 合计 | 99 | 57 | 42 |

生成不可变 combinations.tsv，记录 combo_id、dataset 列表、n_datasets、contains_inhouse、组合哈希、执行顺序及状态。组合编号按固定七位成员 bitmask 或排序后的 dataset_id 集合生成，不用运行完成顺序编号。检验上述 99/57/42 及逐层数量。

先用全七数据集组合做资源和口径 pilot，并将其结果计入 99 组，不另外复制一次正式运行。之后采用预先确定的顺序运行其余组合，例如数据集个数升序、同大小按 combo_id 排序。排序不得依赖 P/q。

全部可计算组合都要运行到主 SpottedPy 结果，不能因聚类图一般、没有 MP10 富集、样本少、无显著结果或不含 inhouse 而提前裁掉。无法估计的组合或关系保留状态、原因和实际数量；“无热点”“目标缺失”“单位数不足”是有效结果状态。

若实际发现缺少某个数据集输入或有确证的跨来源重复，保留原 99 组合登记和问题证据，按已知患者/标本身份正确去重或标注。不得把缺失输入悄悄当成空数据继续，也不得把不同文件夹名自动解释为七个已证实独立的队列。完成未受影响的组合，报告受影响范围与恢复命令。

## 3. 先做轻量资源和身份核对

### 3.1 磁盘、配额与内存

用户特别强调磁盘空间不足，磁盘容量和 RAM 要分别检查。

- 只读检查输出盘和临时盘的 df、inode；检查可读取的用户配额或项目配额。effective_available_bytes 取可获得限制中最小剩余额度。
- quota 不可读时明确记为 UNKNOWN，不能把全盘空闲量等同于用户配额，也不要因此直接放弃；采用小批量写入和保守预留，并捕获 quota/ENOSPC 错误。
- 检查实际可用 RAM、cgroup/调度器限制和已有任务负载；尊重本机调度规则。默认一个组合 worker、一个重计算 section worker；数值库默认单线程，必要时依据 pilot 设少量线程。
- 不同时启动 99 组，也不在每个组合里再展开多进程热点运算。不生成全组织 N×N 距离矩阵或稠密邻接矩阵。

默认新增磁盘预算写入 config：

- 本轮新增持久结果、缓存、临时和图文件合计预算：min(3 GiB, 启动时 effective_available_bytes 的 40%)。
- 缓存预算不超过总预算的 25%；每组合临时文件预算不超过总预算的 15%；图包预算不超过总预算的 20%。这些是上限，不是必须占满的配额。
- 安全预留：min(2 GiB, max(256 MiB, 启动时 effective_available_bytes 的 20%))。
- 对原子写入、压缩及文件替换预留同时存在新旧两份的峰值空间。
- 每个重步骤和大文件写入前重新检查当前可用空间、已用新增预算与预测写入峰值。不要只在任务开始时检查一次。

先根据 pilot 测得的真实文件大小、运行时间和内存峰值更新估算，不宣称固定预算一定足够。若必须保存的最小结果也无法容纳，完成可完成的小结果和断点后标记 RESOURCE_PAUSED；列出缺口与续跑命令，不通过删旧数据、增大磁盘预算、减少组合或放宽方法来绕过。

### 3.2 身份与可用性

使用完整切片表而非仅有评分中心的子表。现有参考总数是 63 张切片；有 v8 score_eligible 中心的是 60 张。这些数量是对照，不是必须强行凑齐的新结果。

核对 patient_uid、analysis_unit_uid、specimen、section、dataset、组织状态、治疗或时间点等已有字段。inhouse 当前七个样本的 patient_uid 尚未完整登记，沿用现有独立样本政策并标注；如项目已有新增可靠映射则应用真实映射，不能凭样本名猜患者。GSE238145 B17/B18 的来源冲突应明确记录并采用既有可追溯政策，不自动拆出新的患者。

多个切片或跨来源记录属于同一患者时，保持同一 analysis_unit_uid；组合总体统计中不能因 dataset 不同将其计两次。对治疗前后或不同组织状态，保留现有主分析范围和既定归并政策；需要改变科学比较时记录为新分支，不暗中改变纳入范围。

来源映射未完全核实并不自动阻断本次已授权的探索。全部单位与已核实患者子集明确分开，后者只能作为预设辅助分析。

## 4. 两条分析分支

### 4.1 REFIT：99 组真正重拟合，作为主搜索分支

每组使用自己的数据集成员，重新生成参考集、拟合 PCA/Harmony、建图、聚类、转移标签并定义新 accepted anchors，完整运行环分析和主 SpottedPy。

### 4.2 FROZEN_LABEL：99 组旧标签子集对照

使用当前全七数据集学到的 SAMPLE_R12_LEIDEN 冻结标签和 v8 的切片效应，按每个组合重新汇总分析单位与总体结果。只要定义完全一致，可直接读取共享旧切片效应，不复制 99 套逐点数据。

这条分支用于区分“样本组成变化”和“重新定义 niche 的变化”。其源域保留旧版定义，不伪称与 REFIT 的源域完全相同。两分支比较时报告实际源域、有效切片和单位覆盖差异。

FROZEN_LABEL 是对照分支，不能用其最小 q 替代 REFIT 主搜索结果。它不是独立验证，旧标签学习已经包含了各队列。其全局统计单独成 family。

## 5. 建立一次性的共享输入

共享输入只保留路径引用、必要列或一次转换后的紧凑数组。不要复制表达矩阵、Seurat/AnnData、cell2location 模型、整个 v2–v8 文件夹或原图。

建议共享内容：

- 唯一 master spot_uid 字典和整数行索引，dataset、section、analysis_unit 映射。
- 完整组织的物理坐标、真实 R1 稀疏图、连通组件以及必要的 R2/R3 邻接或轻量统计。
- 冻结的 30 类 subtype q05/all30 组成与中心 posterior mean Epithelial/all30。
- 现有 R1/R2/R12 原始组成和 Hellinger 特征，原始评分所需的轻量有效性列。
- spot_MP10_model_universe、epi_MP10_model_universe 两层现有评分。
- 可安全共享的完整组织细胞目标热点成员及支持域。

沿用现有 float32 特征时保留其原始身份和精度来源；PCA/需要精确统计的步骤按当前实现升为 float64。效应汇总和所有 P/q 保存为 float64，不为省空间降低推断精度。

## 6. 每个 REFIT 组合的邻域学习

### 6.1 保留实际已认可的输入定义

- 使用 accepted v2_30type cell2location subtype q05，保留冻结的 30 因子顺序，主组成包含 Epithelial。
- 中心 gate：posterior mean Epithelial/all30 >= 0.25；不将其声称为恶性细胞鉴定。
- 真实 Visium 六角网格，中心排除；R1 理论 6 个外围点、纯 R2 理论 12 个；至少 4/8 个有效邻居。
- 缺组成但有物理位置的组织点保留在拓扑图中，不将缺组成填 0，不补不存在的邻居。
- 沿用实际 v2/v3/v7 口径：每 spot q05 按 all30 归一化，圈内有效邻居等权平均，R12 为 R1 与 R2 的等权平均。不要换成邻居原始丰度求和后再归一化，也不要将 18 个点直接等权混合。
- 以全部合格上皮 anchors 学环境。MP10、其高低分组、末端 P/q 和挑选的生物学候选不进入参考抽样、PCA、Harmony、聚类或选解。

### 6.2 参考集与表示

- 按原平衡参考抽样函数，只在本组合内部重新分配各 dataset/section 的抽样额度；优先沿用原 20,000 个参考 anchors 上限，不足时取全部可用值，不重复抽样补足。
- seed = 123。记录实际参考 UID 和配额，不能从全七数据集拟合好的表示中简单删行后称为本组合重新拟合。
- 在本组合参考集重新拟合 R12 Hellinger PCA，使用最小累计 95% 解释方差规则；不额外标准化或 whitening。原全数据集 22 PCs 只是参照，不强制复用。
- 主 Harmony 固定为 SAMPLE/WEAK：batch_key = dataset_id::section_id，theta = 1，lambda = 5。batch_key 代表切片，不能与 patient_uid 混用。
- 其余 Harmony 参数沿用已审计 v3：sigma=.1、内部 nclust=100、seed123、单线程、迭代和收敛参数使用实际版本及有限继续规则。内部 nclust 不等于最终 niche K。
- 在本组合全部共同 anchors 的 PCA 坐标上联合拟合 Harmony，明确属于联合表示学习。使用实际安装版本支持的正确导出 API，不盲目切换软件版本。
- 本轮不在 NONE/DATASET/SAMPLE、WEAK/STRONG 等多套表示中再按末端显著性选优；主搜索轴是 dataset 组合。已有原表示比较结果作为背景保留。

### 6.3 SNN、Leiden 与解的选择

- 沿用原实现：Euclidean 30NN（按原实现包含 self）、shared-neighbor Jaccard、原 1/15 边阈值、无向图、无 self-loop；权重和 tie 处理不另改。
- Leiden RBConfigurationVertexPartition，weighted，n_iterations=10，seed123。
- 沿用原有限 resolution grid：0.005、0.01、0.02、0.05、0.1；扩展仅限原先确定的有限规则，不无限调参以找阳性。
- 按原组成分离、非退化性、适当群数及稳定性规则选择工作解；原 resolution=.1、K=9 是参考，不能强迫各组合均为9群。
- 候选表示与聚类选择只能依赖环境及结构指标。所有候选的简短指标和选择理由保留，不保存各候选整套大对象。
- 若候选都有质量限制，记录实际最优可运行解和限制并继续探索；不能在看过 MP10 后改变质量门槛。不能仅因没有生物学上偏爱的 CAF/B 等群就重分。
- 保留原 15NN 标签转移、support>=.6、距离不超过 reference leave-self-out Q99 等 accepted 规则；参考 native override、tie、UNASSIGNED 规则沿用原代码。
- 主 seed123，稳定性 seeds124–128，沿用参考 80% 子集和 R12 空间块扰动逻辑。此时只重拟合聚类，表示固定；明确不是患者留出验证。稳定性重采样只保存指标，不保存每次完整模型与 embedding。
- UMAP 不参与选解。99 组扫描阶段不做高分辨率 UMAP 或逐切片大图，最终重点组合再生成。

### 6.4 冻结、解释与新旧对应

在接入 MP10 之前冻结每组候选解、最终 K、标签、accepted mask、参数、结构质量指标和哈希。

用原始 q1/q2 组成解释每个 niche，输出每 dataset/section 的覆盖、患者贡献、组成 profile、原始空间分布摘要。保留 sample/dataset 主导、很小群、低支持等原有 flags；这些 flags 不是看结果后删群的工具。

以原始组成距离及共同 anchors 对应关系建立新旧 niche 映射，允许一对多、多对一、未匹配；不能仅用相同数字 ID 匹配。不同组合间的 niche 同样按组成匹配后才能谈组合稳定性。

## 7. MP10 有效域与 v7 环分析

### 7.1 评分及分析域

不重跑 cell2location，不重新训练 MP10 来源模型，不更换基因集、评分分母或 cutoff。按 spot_uid 接入两层已有模型评分和完整性记录。

按实际 v7 定义重建本组合的新 accepted anchors 对应的 S_s、S_score/E0：R1/R2/R3 同中心支持、完整 Top20、epi_score_valid、两层评分有限等条件继续适用。R3 为纯第三圈，理论18个，主要共同中心支持沿用4/8/12。若某位置没有可靠现有评分，标记缺失；不得为扩大样本或改善阳性补值、重训或从不兼容层替代。

必须显式拆出：基础几何和评分 eligibility、当前组合 accepted mask、实际 E0、每一指标实际 E。新旧聚类的 accepted mask 可能不同，不能把旧 E0 无条件视为新 E0。

### 7.2 完整重复当前环分析，控制文件量

- 重复 R1/R2/R3 组成、同切片背景差、相邻圈配对差与连续 MP10 关联。
- 同一比较使用相同实际中心和有效单位；保持 raw/excess/ringpair 的分母和定义一致。
- MP10 主关联沿用 v7 的秩残差化 partial correlation：rank(center mean Epi)、rank(log1p(center_nCount))；原 raw 和中心对应细胞类型调整仅作既定辅助。
- 评分层比较先匹配实际有效切片，再按 analysis_unit 汇总；Fisher z 合并相关系数及主检验沿用实际 v7。v7/v8 不同的数值零规则不得混用。
- 新 niche–celltype 列表根据新环境组成、沿用原 v6/v7 选择规则生成并在 MP10 结果前冻结；不把旧139对按数字标签硬套。若原手工选择无法程序化恢复，统一使用本组合所有 K×29 个非 Epi 配对，明确记录扩展与新 planned family，不自行挑少量显著者。
- 每组保存完整紧凑数值表，不生成成百上千张 boxplot。此前的环关联阳性或阴性都保留；其 q 不用于替代下述 SpottedPy 主搜索命中。
- 为所有组合按各自预定环分析 family 保留组合内结果，按各类同口径环终点另外计算跨组合辅助校正，不与主 SpottedPy family 混接。

## 8. SpottedPy 主计算与正确缓存

### 8.1 固定项目 adapter

以已审计的 v8 实际实现为准：锁定 SpottedPy commit 9dc364f203b7617b6701aa8fb579cbbb35357b26 及 esda/libpysal 实际版本；本项目是明确的 adapter，不声称未修改地复现官网所有默认统计。

- SOURCE_G：esda.G_Local(star=False, transform='R', permutations=9999, seed=123, n_jobs=1, keep_simulations=False)。
- 扫描主尺度 k=10。所有99组均跑此尺度，两评分层、SOURCE_P05和LOCAL_BH05均输出。
- SOURCE_P05：原生较小尾 p_sim<.05 加 Z 方向，未做局部 FDR。
- LOCAL_BH05：min(1,2*p_sim_native)，在每 section/variable/k 所有可计算位置一起 BH，再按 q<.05 与 Z 方向定义。不能单独给 positive-Z 位置做 family。
- 每变量最少30个有效位置、至少3个不同值；支持连通域满足 k+1 位置。原始热点成员和正式区域区分，正式 R1 连通区域至少5个成员。
- 使用完整组织连通组件限制邻居/最近目标，不跨不相连的组织片连接；v8 的置换背景和权重定义不随意改为另一种零模型。
- 坐标单位仍为 full-tissue R1 边中位长度归一化后的 spot pitch；未校准时不得写作微米。
- kNN tie、顺序、坐标 jitter 和 SOURCE_G/G* 区别沿用已审计规则；不把R12的18外围点与热点k10混为一谈。

### 8.2 主终点和支持门槛

源为本组实际 E0 支持域中的 MP10 正式热点；niche 目标为本组全部 accepted common anchors 中对应 label，不限于 E0。细胞目标为完整组织 q05/all30 的29个非Epi相对组成热点。

沿用 v8：

- O = H 中目标指示变量均值 − E 中目标指示变量均值；正为重叠增加。
- D = H 到最近目标的平均距离 − E 到最近目标的平均距离；负为更近。
- D_out：在其实际域中排除已经属于目标的位置后重算；不重定义源热点。
- 最近目标限制同一物理组织组件；没有目标返回 NA，不能填最大距离。target niche 缺失时 O 也记录不可估计，不以必然0充当支持。
- 点权重均值为主；组件等权、median等按原方法作描述性敏感性。
- O 推断要求 E>=30、H>=20、nonhot>=20、target-in-E>=20、rest>=20。
- D/D_out 要求其实际域 E>=30、H>=20、nonhot>=20、target-full>=5。
- 每切片先算效应，同一 analysis_unit 内等权合并切片；总体不能把spot、热点组件或技术切片当独立患者。
- 两侧 Wilcoxon signed-rank vs0，至少5个单位和5个非零效应；全零且单位数足够按原规则 P=1。v8 零容差1e-12；exact/approx、ties和continuity处理沿用原代码。
- 同时报告mean/median、N、非零N、正负单位数；不要把signed-rank检验说成专门检验均值。保留其分布和独立性假设。

### 8.3 缓存不得改变统计对象

可共享：原始物理图、坐标、分母及完整组织细胞目标热点。缓存键必须包含输入内容哈希、排序后的UID、domain/support mask、变量、k、权重/组件政策、版本、seed、permutations、P值约定、轨道及组件大小规则。

MP10源热点缓存尤其必须包含实际源域UID或mask哈希及其形成政策。禁止只按 section/layer/k 命中缓存：新 accepted mask 可能改变 E0、kNN、支持域、P值背景和局部BH数量。

标签目标变化后重算 niche 几何。只有源域、目标域、权重、评分和参数完全相同，才能重用效应。旧两层公平比较也必须按实际匹配域重做。

逐 section/variable 计算所需统计后，保存 compact formal-hot/support masks、区域摘要和最终切片效应，及时释放点级Z/P/q及临时数组。不得保存9999次置换的逐点模拟分布。

如能安全复用 v8 的置换索引 memoization，先做小规模 bitwise/parity 对照，并限制该缓存大小；不能在RAM或磁盘无限保留每一个n的置换索引。

## 9. 穷举后的 BH/BY 口径

### 9.1 明确主要搜索 family

主要搜索只使用：

`REFIT × SOURCE_P05 × 主k10 × 全部既定分析单位 × 两层MP10 × 所有冻结niche × O/D`

第c个组合最终有 K_c 个 niche，则主 planned family 大小：

`M_primary = sum_over_combinations(4 * K_c)`。

该family覆盖本轮99组，不仅覆盖57个含inhouse的组合。局部BH热点轨道、细胞目标、D_out、已核实患者子集、旧标签对照和环关联均为预设辅助，不以其中最小P值补救主搜索。

### 9.2 保留三种校正值

对每条主关系保存：

- p_raw：原始分析单位层关系P。
- q_within_BH：本组合4K_c个计划关系的BH。
- q_global_BH：本轮所有组合的全部主关系原始P一次性BH。
- q_global_BY：相同全局family的Benjamini–Yekutieli校正。

全局池必须包含所有预定主关系；不得先用组合内q或P筛选、不得只收集每组合最小P、不得对q_within再做BH。不得假定q_global_BH数值必然大于q_within_BH。

只有必需计算已经成功完成、关系因预定的目标缺失或支持数量不足而不可检验时，才保存原p/q为NA和明确原因，按现有planned-n方法用P=1占位参与校正分母，不把占位1冒充观测P。技术失败、计算尚未完成或结果文件损坏不能用P=1补齐，必须与科学上不可估计的状态分开。中途的全局q仅为PROVISIONAL。全部组合完成或按事先规定的、与MP10结果无关的结构性不可估计规则结算后，才输出FINAL全局q。

任何组合的必需计算因技术问题未完成，无论最终K是否已经确定，都标记SEARCH_INCOMPLETE，继续未受影响任务并保留恢复信息；此时全局结果只能是PROVISIONAL。例如，标签已冻结但SpottedPy尚未成功计算，不能按普通NA结算。最终K未确定时也不能随意指定K=9补分母或不留痕地删除该组合。对结构性无任何可定义niche的情况保留明确状态及无检验原因，不制造伪检验。

BY针对检验间任意依赖的保证仍以每个原始P有效、family完整为前提。它不能修复患者重复、检验假设不合适、空间零模型不足、共享数据学习及未处理选择导致的P值失真。

flat全局校正对应的是整个候选关系拒绝集合。它不保证最后挑中的一个组合内部FDR，不消除赢家效应偏差，也不能把重点短名单当成每条关系95%可信或独立验证。

### 9.3 辅助family

为LOCAL_BH05的niche O/D、各轨道niche D_out、cell D、cell O、cell D_out、verified-only、FROZEN_LABEL、环统计分别建立命名清楚的family；保留各自组合内校正和同family的跨组合辅助BH/BY。REFIT家族分母随对应最终K/固定29类型确定，缺失按planned-n规则处理。

FROZEN_LABEL始终保留旧版niche0–8的九类计划目标；SOURCE_P05的O/D对照每组合固定36项，全部99组合固定3564项。某子集没有某个旧label时，仍保留该计划目标并按成功评估后的目标缺失规则处理，不按子集中实际出现的label数缩小family。其LOCAL_BH05和其他辅助family同样按固定旧9类或固定29细胞类型计算计划数量。

这些辅助family用于解释稳健性与其他问题，不能在主family未命中时把它们改名为主命中。各dataset结果用于分层展示；不要把组合总体、单dataset和多个辅助scope混成可任意挑星号的主榜。

## 10. 全部结果榜、inhouse优先榜和重点对象

先输出完整结果，然后再输出便于研究讨论的短名单。

必须有：

1. 99组合的完整结构与运行汇总。
2. 全部主关系结果表，包含未显著与不可估计关系。
3. 全部q_global_BH<.05的主关系，不因不含inhouse或方向不符合预期而隐藏。
4. contains_inhouse=TRUE的57组合优先榜，以及42个不含inhouse的对照榜；保留从全部99组计算的global q，不为优先榜重新缩小主family。
5. 仅q_within_BH<.05的组合内探索候选表，与全局命中明确区分。
6. q_global_BY<.05支持标记。

“统计显著”与效应方向分开。O>0为重叠增加；O<0为重叠减少；D<0为更近；D>0为更远。两种方向都记录。不得把GC B历史的负距离差说成回避，也不得把重叠减少直接写成主动免疫排斥。

重点展示最多3个不同组合，每组最多3个重点关系。按以下规则事先固定选择：优先全局主命中；同类中优先含inhouse且inhouse可估计/同向、跨dataset和患者支持较好、结构质量较好的对象；再参考效应和q及组成差异避免三份完全重复的例子。BY通过是支持标记，不要求为了过BY不断改参数。

若没有全局命中，最多选择3个组合内探索候选明确标记；若连组合内命中都没有，展示最有信息的非显著结果及其N/效应，不将趋势写成阳性。这种展示排序不改变全部统计表和全局拒绝集合。

inhouse“被纳入组合”与“支持该关系”是两个字段。分别报告该niche在inhouse的覆盖、可估计切片/单位、效应方向、分布及是否为NA。不能用含inhouse这个标签掩盖该关系在inhouse完全不可估计。

不要要求每个dataset都单独BH显著，也不要求留一dataset每次q<.05。99组合共享大量患者，很多组合同向只能称组合稳定性，不能当成99次独立复现。

## 11. 重点对象的后续检查与图包

99组主扫描结束后，对最多3个选中组合自动完成以下细查，不必再等待用户逐组审批：

- 各dataset、inhouse和全部分析单位的效应分布与可估计N；两评分层按相同实际切片/单位公平比较。
- 冻结该组合标签的leave-one-dataset-out影响诊断：从保留的切片效应重新按剩余analysis_unit归并，正确处理跨dataset共享单位；不是重新学习的holdout，也不是独立验证。
- 既有k6/k18、技术分层、几何支持、组件等权、完整环等敏感性。k6/k18主用效应和支持变化描述，不制造又一套可择优的主P值。
- 必要时对少数关键技术问题复算，不展开新的统计方法、阈值、模型或无界参数搜索。
- 显示全部niche组成热图、dataset/患者贡献、reference UMAP、目标niche与MP10空间图、O/D/D_out分开的森林图或箱点图，以及REFIT/FROZEN_LABEL差异。
- inhouse七张切片保留紧凑overview（无有效结果也写明）；其他每dataset默认选择最多两张技术覆盖合格、具有代表性的切片，选择依据记录。完整切片/单位数值表必须保留。
- 只对重点panel保存透明300dpi PNG和矢量PDF；其余为紧凑HTML或数值汇总。不要为99组合×全部切片×全部细胞目标导出双格式图包，不自动打包大型ZIP。

普通置信区间应按分析单位而非spot重采样，必要时在dataset内分层；只用于重点对象的描述，标注其未经赢家选择校正。不要用新增bootstrap或另一检验的P来替换原主P。

## 12. 必须执行的节省空间策略

### 12.1 永久保留的最小复现信息

所有99组合都保留，而不只保留胜出组合：

- dataset成员、参考整数行索引、输入/配置/软件和代码哈希、seed、选解记录、最终K与质量指标。
- 本组合anchor整数索引、candidate/accepted label（带-1的int16或安全的更紧凑编码）、accepted/support mask；必要的support和距离诊断列压缩存储。
- 原始组成profiles、新旧标签映射、覆盖及状态计数。
- 最终切片效应和analysis_unit效应、P/q、planned family清单及不可估计原因。
- 运行时间、内存/磁盘峰值、有效域哈希、缓存命中/失效原因、断点状态及清理日志。

使用共享UID字典，避免在每个大型表反复写长字符串。Parquet启用字典编码和适度ZSTD压缩；小表用TSV/JSON。一个数值结果只保留一个权威版本，报告按需读取，不同时保留大CSV/TSV/Parquet三份。

标签保存成本很小：按129,632 anchors全长int16保留99份，理论约24.5MiB，另加少量mask。因此不能以省磁盘为由删掉所有失败/阴性组合的标签。

### 12.2 不保存或及时删除的可再生大对象

- 每组合完整AnnData/Seurat、表达、c2l对象副本。
- 全量PCA/Harmony embedding、SNN图、全部候选模型、每次稳定性拟合对象。
- 全部逐点Z/P/q、全量point-to-target距离、每次置换模拟矩阵。
- 99套UMAP/HE图、所有celltype图、缓存的图片副本、重复ZIP/压缩包。

可保留很小的PCA模型、参数或参考摘要帮助复现；大表示和图在当前组合使用完后释放。只为重点组合按预算保留必要模型，其余通过冻结输入、UID、代码及环境重新生成。提供重新生成指定组合图或模型的脚本，不在报告里写不存在的命令。

### 12.3 流式执行、缓存淘汰与安全清理

- 一次只处理一个组合和当前需要的section；先写最终轻量效应，再释放点统计与临时数组。
- 共享cell目标缓存优先保留；源域依赖缓存采用有大小上限的LRU，预算不足就重算，不能跨错误domain复用。
- 临时目录须在本轮输出目录内，带run_id和ownership marker；禁止临时目录默认落到空间更紧张的用户home或不受监控的/tmp。
- 每阶段成功提交、验证必要结果后，立即清理本阶段本轮自产且可再生临时大文件。
- 清理必须先检查realpath、ownership marker、当前运行引用及保留清单；只删除本轮明确登记的文件，不跟随指向旧数据的symlink，不使用宽泛rm -rf或按文件名通配跨目录删除。
- 删除前记录相对路径、字节数、所属阶段、可再生依据；保留小日志。失败阶段先保存参数和错误摘要，再清理其无用大temp，不能删唯一有效结果。
- 单日志滚动截断或轮转，避免海量逐点print；保留完整关键错误与版本信息。

若按最小方案仍超预算，先省去普通图、重复报告和可再生缓存，保留完整数值、标签和复现信息。仍不够时安全断点暂停，不能悄悄只跑一部分组合却标记完成。

## 13. 断点与自动运行规则

每个组合采用明确阶段状态：REGISTERED、REPRESENTATION_DONE、NICHE_FROZEN、RING_DONE、SPOTTED_MAIN_DONE、COMPLETE，另有INELIGIBLE、TECHNICAL_FAILED、RESOURCE_PAUSED。

使用原子写入：先写本轮temp，再核对schema、行数、hash、必要数值，最后提交阶段状态。状态不能先写COMPLETE再写文件。组合级完成标志必须依赖必需结果，而非仅依赖进程退出码0。

续跑需同时核对输入哈希、config hash、代码版本、master UID顺序、domain hash与产物状态；已有轻量结果一致时不重复计算。已完成但大embedding因策略删除，属于正常完成状态，不因为缺少可再生大对象自动重跑。

中途发生单组合技术错误，保留证据并继续其他可运行组合；有限修复和重试需记录，不改变科学参数。资源暂停需保存当前最小断点并退出重计算，不能忙等或自动删除旧项目。

不要在每组niche冻结、每阶段或每个辅助结果后停下等待用户确认。无阳性、来源映射尚不完整、某些目标NA都不是审批停机理由。只有真实输入缺失、无法安全写入、无法执行的技术错误或会改变已授权科学范围的重大矛盾需要明确报告未完成项。

## 14. 必要验证，不做重复大规模测试

先用小规模真实切片验证，再执行扫描；测试聚焦可能改变结果的风险：

1. 组合登记准确99/57/42，无重复或遗漏。
2. 新组合只用成员dataset重新拟合；无MP10进入聚类或选解的代码路径。
3. 原始ring计算口径、UID join、mean gate和q05组成与旧代码一致。
4. 同一热点输入直接计算与缓存读取一致；改变source mask后必然缓存失效。
5. 冻结标签子集从旧section表重新汇总，正确归并重复患者。
6. 无目标、空热点、跨组件、重叠零距离、数值零和exact/approx检验规则正确。
7. 全局BH/BY直接使用所有原始P，planned缺失占位、不同K和不完整搜索状态正确；核对一种已安装统计库的结果。
8. 模拟一次小规模断点恢复和本轮temp清理，确认不会删除旧文件、必要结果或在引用缓存。

已有v8稳定通过的算法无需对所有组合重复执行整个验证套件；复用其证据，只验证这次组合重拟合、域缓存和全局汇总引入的风险。科学数值比对的容差应基于原算法/精度规定，不为通过测试任意放宽。

## 15. 最终交付与状态

至少提供：

- README.md：本任务范围、实际路径、运行/续跑/单组合重建/出图/仅清理本轮temp的真实命令。
- RUN_SUMMARY.md：99组合状态，57/42覆盖，主全局BH命中、BY支持、组合内探索命中，重点对象和限制。
- combinations.tsv、combination_summary.tsv。
- primary_relations.parquet：全部主关系与p/q，完整family_id、combo_id、K、track、layer、endpoint、target、effect、N和状态。
- global_primary_hits.tsv、inhouse_priority_candidates.tsv、without_inhouse_controls.tsv、within_combo_exploratory_hits.tsv。
- 独立的辅助family结果与REFIT/FROZEN_LABEL对照。
- 所有组合的compact labels/indices、profiles、section/unit效果及manifest。
- 最多三个重点组合的review图包和全局overview。
- resource_usage.tsv、cleanup_manifest.tsv、validation_summary、RESUME.json。

结果分类：

- GLOBAL_BH_HIT：REFIT主family的q_global_BH<.05。
- GLOBAL_BY_SUPPORTED：同一主关系q_global_BY<.05的支持标记。
- WITHIN_COMBO_EXPLORATORY：只有组合内BH达到阈值，明确其搜索限制。
- NO_GLOBAL_HIT：本轮没有整轮搜索校正后的主命中，但仍交付完整结果。
- SEARCH_INCOMPLETE：有技术性未结算组合，不能宣称完成全部穷举。

执行任务完成与科学目标达到分开记录。全部99组完成且没有全局阳性，也是执行完成；不要为获得阳性额外新增阈值、模型或继续无界搜索。

最后向用户直接回答：跑完多少组合；哪些含inhouse组合值得关注；其阳性属于哪一层校正；inhouse是否实际支持；哪些结果由单一dataset/患者驱动；新旧标签结论是否一致；新增磁盘占用及已清理多少本轮临时文件；仍有哪些未完成问题。

不要称“BH保证至少一个真阳性”“已经消除所有选择偏差”“99次独立验证”。不要把MP10高分直接等同于胆固醇浓度、恶性上皮特异机制、细胞招募或主动排斥。

## 16. 只读方法与既有结果参照

- 当前niche方法：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/neighborhood_ring_harmony_v3_upload/reports/METHODS.md
- 当前冻结niche清单：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/neighborhood_ring_harmony_v3_upload/frozen_niche_manifest.json
- 完整切片身份：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/neighborhood_ring_harmony_v3_upload/00_preflight/batch_key_mapping.tsv
- 当前v7环分析：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/review_B/reports/METHODS.md
- 当前v8方法：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/neighborhood_spottedpy_v8_upload/reports/METHODS.md
- 当前v8结果：https://github.com/Kimsx05/transfer_for_neibor/blob/main/results/neighborhood_spottedpy_v8_upload/reports/RUN_SUMMARY.md
- SpottedPy原文：https://pmc.ncbi.nlm.nih.gov/articles/PMC11552145/
- SpottedPy锁定源码：https://github.com/secrierlab/SpottedPy/tree/9dc364f203b7617b6701aa8fb579cbbb35357b26
- Benjamini与Yekutieli，依赖条件下FDR方法：https://doi.org/10.1214/aos/1013699998
- Benjamini与Bogomolov，选择后多family推断：https://doi.org/10.1111/rssb.12028

这些链接用于核对出处；实际执行以服务器中可追溯的已审计输入和实现为准。发现差异时记录并据实处理，不根据一段方法摘要另写不兼容算法。
