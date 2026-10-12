# BLCA Visium 双环 niche Harmony 比较执行 Prompt v3

版本日期：2026-10-12  
执行目标：比较不整合、按 dataset 整合、按 sample 整合，完成上皮 anchor 的 niche 发现、全量标签应用、解释与复核材料。  
本轮终点：niche。MP10 评分、关联、解盲和来源验证全部后置。

请完整阅读本文件后，在服务器上实际执行。不要只回复计划、可行性或预检结果。例行适配、有限选参、可逆修复和兼容 checkpoint 续跑均已包含在授权内；不需要每完成一个阶段就询问是否继续。

## 0 本轮决定与停止位置

用户已经确认：

1. 并列比较 NONE、DATASET、SAMPLE 三种整合模式。SAMPLE 必须实际运行，不能只在讨论中列为备选。
2. 分别分析 R1、R2、R12，每种输入均比较 K-means 和 Leiden。
3. 最终主比较为 3 种整合模式 × 3 种环输入 × 2 种算法，共 18 个配置。
4. 先用 R12 比较两档 Harmony 参数，再在 DATASET 内选一档、SAMPLE 内选一档；两个模式均扩展到 R1 和 R2。
5. 保留 v2 原始结果。v3 修改整合与选参，不覆盖旧标签、不改原始细胞组成或物理环。
6. 区域复核材料要告诉用户“看哪张图、判断什么、填写什么、如何返回”；区域未审完不阻塞本轮 niche。
7. 本轮不读取或分析逐 spot MP10，不运行 MP10 来源验证、pseudo-spots、DESTVI 或旧流程的 MP10 模块。
8. niche 结果和复核材料完成后停止，交付供用户审阅。

建立独立的 v3 入口和 scope manifest。不要直接运行旧版可能继续进入 MP10 的 resume.sh。新入口正常退出时必须仍满足：

~~~json
{
  "stop_after": "niche",
  "mp10_authorized_this_run": false,
  "mp10_values_loaded": false,
  "mp10_analysis_executed": false,
  "source_validation_status": "DEFERRED_BY_USER",
  "region_review_required_for_niche_completion": false
}
~~~

可以读取现有数据的字段名和非 MP10 元数据；加载表时使用字段白名单，避免把 MP10 数值一并带入分析对象。不要为了制作 README 或完成状态，调用任何旧的 MP10/来源验证脚本。将来源验证保留在后续事项中即可，不在这轮补做。

## 1 工作目录与已有成果

主项目：

~~~
/home/data/t070721/codex_workspace/Bladder Metabolism
~~~

已完成 v2：

~~~
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2
~~~

建议新目录：

~~~
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3
~~~

其他定位线索：

- 已认可 c2l：主项目下 cell2location_full_run/v2_30type。
- 现有 CNV：主项目下 spatial_ecology_top20_v2/03_copykat_geography。
- v2 的 00_preflight、01_regions、03_features、04_clustering、05_comparison、reports、code。
- GitHub 回传目录：
  https://github.com/Kimsx05/transfer_for_neibor/tree/main/results/neighborhood_ring_comparison_v2
- GitHub 当前回传存在同名双层目录；在服务器上根据 manifest 和实际文件定位，不照搬重复路径。

先读适用的 AGENTS.md、README、v2 input_contract、特征字典、冻结清单和实际脚本。本轮用户决定取代旧流程中与本轮范围冲突的要求，例如“全 ROI 审完才允许 niche 完成”或“niche 后自动进入 MP10”。

现有审阅基线如下，运行时核对实际数据，不能为了凑数伪造：

| 项目 | v2 审阅值 |
|---|---:|
| 目标 datasets | 7 |
| 目标合格切片 | 63 |
| 共同 anchors | 129,632 |
| 完整 R1=6、R2=12 的共同 anchors | 108,345 |
| 固定 reference anchors | 20,000 |
| reference 实际覆盖切片 | 61 |
| 共同 anchors 内有来源证据的患者 UID | 27 |

区分目标切片、具有共同 anchors 的切片和进入 reference 的切片。若数量与此不同，先定位版本、UID 或过滤规则差异；不要静默重抽 reference 或改变共同集合。

资源执行要求：

- 先检查磁盘、内存、兼容环境及当前其他任务。不要升级整个环境。
- R1/R2/R12 按需顺序处理；大型任务控制并发，不干扰正在运行的其他任务。
- 复用旧 feature/cache，必要时只复制小型清单；只读引用的文件不得就地修改。
- 不重跑 63 张 c2l、CopyKAT，不恢复或导出历史大 CNV 文本矩阵。
- 不构造全量 anchors × anchors 稠密距离矩阵。SNN 用稀疏矩阵，silhouette 用固定小型评估集。
- 使用 seed=123；稳定性种子另行固定。保存实际命令、环境、参数、输入哈希和 checkpoint。
- 若 v3 已有兼容 checkpoint，复用；规格不兼容时新建版本子目录，保留此前输出。

## 2 冻结输入与样本身份

### 2.1 唯一主输入合同

严格沿用用户已确认的：

- 环境：accepted v2_30type 产生的 frozen subtype q05 相对组成。
- 中心上皮 gate：accepted posterior mean 的 Epithelial/all30 >= 0.25。
- 不把 q05 环境和 mean gate 的组合视为错误混搭；这是已批准的口径。
- 不把环境主输入替换成 mean，不改变 gate，不把 Epi 当作恶性细胞比例。
- 30 类名称、顺序、映射、共同 anchors、完整环标志和 reference UID 均与 v2 一致。

读取并保存 v2 input_contract 和 feature_columns 的副本与哈希。不得按矩阵当前列位置盲拼，不得把缺失输出填成 0。

### 2.2 SAMPLE 的操作定义

本轮 SAMPLE 指“独立空间测量的切片单位”。构造：

~~~
sample_uid = dataset_id + "::" + section_id
~~~

其中 section_id 必须对应一个实际切片/独立空间测量单元。同一患者的不同切片可有不同 sample_uid；这不代表它们是不同患者。

建立 batch_key_mapping.tsv，至少包括：

~~~
dataset_id
section_id
sample_uid
original_sample_id
specimen_uid_if_known
patient_uid
analysis_unit_uid
analysis_unit_type
source_patient_identity_known
n_common_anchors
n_reference_anchors
sample_key_definition
~~~

original_sample_id、specimen_uid、patient_uid 可缺失。使用来源明确的现有分析单位政策，不凭文件名、数字顺序或 N/T 名称编造患者和配对关系。匿名样本的批次整合不依赖公开患者别名是否恢复。

必须明确：

- DATASET：仅传 dataset_id。
- SAMPLE：仅传 sample_uid。
- 不把 dataset_id 和 sample_uid 同时传入一次 Harmony。
- 不先做 DATASET，再在其结果上连续做 SAMPLE。
- 不用 patient_uid 替代 sample_uid。
- 不在每个 dataset 内独立运行 SAMPLE 后拼接结果。
- 每个模式均在全体合格 anchors 上联合运行；无 anchors 的类别不创建占位零向量。
- 若现有 specimen 与 section 存在多对一关系，如实记录；本轮 sample 校正仍是测量切片层级，患者/标本身份单独保留用于解释。

患者和临床元数据核对仅做到支持本轮来源描述与稳定性分块所需的程度，不将其扩展成新的单细胞来源验证任务。

## 3 复用环特征和同一套初始 PCA

### 3.1 物理环

复用已经修复缺失组成节点后的全组织物理图及独立环：

- R1：精确一跳邻居。
- R2：精确二跳邻居，排除中心和 R1。
- 中心 anchor 不加入环境向量。
- 有效环组成最低覆盖要求沿用 R1 >= 4、R2 >= 8。
- 完整环敏感性集合为 R1=6、R2=12。
- 物理邻接按切片建立；切片之间没有物理连边。
- 环仍可跨越候选肿瘤边界；不先裁肿瘤 mask 再计算邻域。
- 不新增方向、角度分箱、旋转、极性或其他空间特征。

### 3.2 组成与 Hellinger

沿用 v2 的定义。设 a_jc 为 spot j 的 q05 丰度：

~~~
p_jc = a_jc / sum_c(a_jc)
q1_ic = mean(p_jc for valid j in R1(i))
q2_ic = mean(p_jc for valid j in R2(i))

X1_i  = sqrt(q1_i)
X2_i  = sqrt(q2_i)
X12_i = concatenate(sqrt(q1_i), sqrt(q2_i)) / sqrt(2)
~~~

先对每个 spot 归一化，再在环内等权平均。不要改成先把丰度相加后归一化，不要把 R1+R2 混成一个 30 维均值，不要把双环权重改成邻居数 6:12。

q05 派生的比例是相对组成摘要，不是观测细胞数、细胞存在概率或恶性细胞比例。

### 3.3 PCA

每个输入分别拟合或复用一套 PCA：

- 在固定 20,000 reference 上拟合。
- 中心化、不按列标准化、不白化。
- 保留累计解释方差达到 95% 的最少 PCs。
- “中心化”是减去特征均值，不是使用中心 anchor 的数据。
- 用同一组均值、载荷投影该输入的全部共同 anchors。
- 同一输入的 NONE、DATASET、SAMPLE 共用这套初始 PCA。
- R1、R2 各自做 PCA 和 Harmony，不能从 R12 的校正坐标截取列冒充单环。
- 不往特征中添加中心 Epi、中心表达、CNV、患者标签、坐标或 MP10。

v2 的 PC 数 18/17/22 可用于核查，不硬编码。若输入完全一致而出现明显变化，定位原因。保存载荷、均值、解释方差、行 UID、列名与精度；关键参考坐标和建图输入优先保留 float64，避免后续只能近似重建图。

## 4 比较矩阵和计算顺序

### 4.1 R12 pilot 共五种表示

“pilot”表示先在 R12 上选取整合参数。聚类及选参使用同一批固定 reference；为保证全量标签应用的坐标一致，Harmony 本身在全部共同 anchors 的 R12 PCA 上联合拟合。

| 表示 ID | Harmony 输入批次 | theta | lambda |
|---|---|---:|---:|
| NONE | 完全跳过 Harmony | — | — |
| DATASET_WEAK | dataset_id | 1 | 5 |
| DATASET_STRONG | dataset_id | 2 | 1 |
| SAMPLE_WEAK | sample_uid | 1 | 5 |
| SAMPLE_STRONG | sample_uid | 2 | 1 |

WEAK/STRONG 是本项目预定的两组相对参数名称，不是官方最佳预设，也不保证所有指标单调。相同参数在不同批次数量下不代表相同校正幅度。

四个 Harmony 运行都从相同、未校正的 R12 PCA 开始。NONE 必须实际跳过算法，theta=0 不能用来冒充 NONE。

### 4.2 保留三模式并扩展到 18 个主配置

在 R12 pilot 中：

1. DATASET_WEAK 与 DATASET_STRONG 比较，选定 DATASET 的一套工作参数。
2. SAMPLE_WEAK 与 SAMPLE_STRONG 比较，选定 SAMPLE 的一套工作参数。
3. 保留全部 pilot 结果、参数和理由。
4. 两个模式均继续扩展到 R1、R2。不要只留下全局“赢家”。

每个整合模式的工作参数统一用于三种环输入和两种聚类算法；不得每个 ring、每种算法重新挑一档使其看起来更好。

最终主配置 ID：

| 整合模式 | R1 | R2 | R12 |
|---|---|---|---|
| NONE | NONE_R1_KM、NONE_R1_LEIDEN | NONE_R2_KM、NONE_R2_LEIDEN | NONE_R12_KM、NONE_R12_LEIDEN |
| DATASET | DATASET_R1_KM、DATASET_R1_LEIDEN | DATASET_R2_KM、DATASET_R2_LEIDEN | DATASET_R12_KM、DATASET_R12_LEIDEN |
| SAMPLE | SAMPLE_R1_KM、SAMPLE_R1_LEIDEN | SAMPLE_R2_KM、SAMPLE_R2_LEIDEN | SAMPLE_R12_KM、SAMPLE_R12_LEIDEN |

这 18 个配置里，同一“模式＋输入”的 KM 和 Leiden 必须使用同一套 PCA/Harmony 坐标。

v2 的 full-Hellinger K-means 另列 LEGACY_RAW_FULL_KM 基线。可复用 K=2 和已有 K=6 候选结果作桥接比较，但不能用它替代 NONE_PCA_KM，也不能把它计入 18 个主配置。

## 5 Harmony 实现与参数

优先使用服务器已安装、兼容的 R harmony 矩阵接口，避免构造包含表达和 MP10 的大型 Seurat 对象。使用其他实现时必须说明版本、参数对应关系及数值一致性范围，不能静默替换。

概念调用如下，实际参数名按已安装版本的正式接口适配：

~~~r
set.seed(123)

# Z: anchor × PC；metadata: 相同行顺序的 anchor 元数据
stopifnot(identical(rownames(Z), rownames(metadata)))

Z_harmony <- harmony::RunHarmony(
  data_mat = Z,
  meta_data = metadata,
  vars_use = batch_column,
  theta = selected_theta,
  lambda = selected_lambda,
  sigma = 0.1,
  nclust = 100,
  max_iter = 20,
  early_stop = TRUE,
  ncores = 1
)
~~~

必须做到：

- 入参与返回坐标都明确为 anchor × PC；检验维数、行 UID 顺序、有限值和批次类别。
- Harmony 不负责自动按 UID 对齐 metadata，调用前必须显式重排和断言。
- 显式设置 lambda=5 或 1，不能用 lambda=NULL 自动估计代替本轮比较。
- nclust=100 是 Harmony 内部拟合分量数，和最终 niche K 无关。
- 对安装版本展开记录全部实际参数。当前接口常用 max_iter；旧版可能使用 max.iter.harmony，不能把不支持的参数悄悄忽略。
- 如果正式接口支持，显式固定高级收敛选项 max.iter.cluster=20、epsilon.cluster=1e-3、epsilon.harmony=1e-4；通过该版本的 harmony_options/.options 传入，并记录其余实际值。
- 保存 objective/convergence 信息。达到迭代上限不等于收敛；确有未收敛时，仅对该配置做一次上限 50 轮的延长，参数强度不变。
- 不对校正后的坐标额外逐 sample 标准化，不再次自动运行 PCA 改变主比较表示。
- 不将 Harmony 坐标平方、反变换后称为“校正后的细胞比例”。

资源允许时可以明确限制少量 BLAS 线程，但所有模式设置一致，避免与外部并行过度叠加。记录硬件、线程和软件差异；seed 不等于跨所有平台逐位相同。

因为每个 Harmony 表示都在全量 anchors 上拟合，本轮属于联合整合后在平衡 reference 上聚类并传播标签。剩余 anchors 不能称为未参与表示学习的独立 holdout。

## 6 聚类和粒度选择

### 6.1 K-means

在当前模式的 reference PCA/Harmony 坐标上拟合。

- 工作候选 K = 6、8、10、12。
- 保留 K=2、4 作为粗粒度诊断；它们不能通过旧的宽容差规则自动胜出。
- seed=123，n_init=50，max_iter=300，Euclidean/Lloyd；记录实际收敛。未收敛时可对该候选有限延长，不能当成完成。
- 保存每个候选的 reference 原生标签、质心、组成画像和支持表。

明确废止 v2 的“距最佳 silhouette 不超过 0.02，则选最小 K”规则。v2 三种 full-Hellinger 输入的最高 silhouette 均为 K=6；这个事实用于安排优先检查，不代表 Harmony 后必须仍选 K=6。

### 6.2 Leiden

主比较沿用 v2 的图定义以控制额外变量：

- reference 在当前特征表示中联合建图。
- Euclidean 30NN，沿用包含查询自身的约定，即通常为 29 个其他点。
- 二值邻居集合交集 → Jaccard SNN。
- 保留权重 > 1/15；无向图、无自环。
- leidenalg.RBConfigurationVertexPartition，weighted，seed=123，n_iterations=10。
- 不能把 SNN 也按 sample 分开；按切片分开的只有物理图。

低分辨率主网格：

~~~
0.005, 0.01, 0.02, 0.05, 0.1
~~~

先检查连通分量、孤点和度分布。既有未校正 PCA 的归档复核中 R1/R12 各一个连通分量，R2 两个；当前没有“必须先改 k 或剪枝才能得到少于20群”的证据。

目标是找到中等粒度、有组成区别的工作结果，优先审视约 6–12 群；不是强制输出指定群数。若所有低端候选少于6群，可有限补 0.2、0.4；若都超过12群，可在已有低端附近最多补三个对数间隔候选。保存实际网格与触发理由，不展开无限分辨率搜索。

若需要使用其他 Leiden 接口，必须核对目标函数和图；不能把另一种 partition 或不同 resolution 定义当成同一个配置。

### 6.3 选参应考虑什么

对每个候选至少保存：

- 实际群数、各群大小及极小群。
- 固定评估 UID 上的 silhouette，注明使用当前 embedding 还是原始 Hellinger。
- 两环原始组成均值、群内离散程度、群间画像相似性。
- dataset/sample/分析单位支持及单一来源占比。
- 标签传播的低票率、距离异常率和覆盖率。
- 对入围工作候选的稳定性结果。

可以同时保存 silhouette_embedding 和 silhouette_original_Hellinger，但不得混成一个指标，也不以跨表示 silhouette 最高作为唯一选型依据。UMAP 形状不参与参数选择。

工作选择按以下顺序形成书面理由：

1. 排除实际失败、非有限、退化或未收敛的候选。
2. 在可解释的中等粒度候选中，检查新增群是否带来有支持的组成区别。
3. 结合稳定性、覆盖率与跨来源一致性，而非仅看混合程度。
4. 若 K=8/10/12 相对 K=6 没有清晰增益，优先保留 K=6；有增益则记录新增环境及支持。
5. 如只能建议更粗或更细结果，如实说明，保留已完成候选，不凑群数、不假造生物学命名。

source-dominated 是复核标志，不是自动删除规则；极小群也先记录，不为了达到目标群数静默合并。

### 6.4 同 K 比较

每个“模式＋输入”报告 KM 与 Leiden 的 ARI/NMI，同时列出 K。若 Leiden 工作 K 在 2–20，复用或最多补一个相同 K 的 KM，且使用完全相同的表示，提供同 K 诊断。

群数不同的低 ARI不能直接解释成某个算法失败。旧 full-Hellinger KM 与新 NONE_PCA_KM 的差异另列为表示变化，不与整合收益混在一起。

## 7 整合评估与两档参数选择

### 7.1 同时评估来源混合和组成保留

每个 pilot 表示及最终工作表示均输出：

1. 30NN 去自身后的同 dataset、同 sample 邻居比例。
2. SNN 的同来源边比例、同来源边权比例、连通分量和孤点。
3. 局部来源熵或 iLISI 等可复现指标，分别对 dataset、sample 计算。
4. 每个 niche × dataset、niche × sample 的原始 q1/q2 画像、样本量和离散程度。
5. 全体与分来源的标签分配覆盖率，检查困难样本是否被大量拒识。
6. 对切片内空间块抽样/稀疏子集的来源诊断，检查重叠 rings 是否放大同来源近邻比例。

固定相同评估 UID、近邻数和分母。原始 Hellinger 距离可用于衡量校正后的邻居是否混入组成差异很大的 anchors；在 source 内与跨 source 邻居上分别描述。

dataset 与 sample 的类别数不同，它们的熵或 iLISI 数值不能直接横比。随机来源配对比例只作参考基线，不充当 batch 显著性检验。

尽量在独立已有病理区域或预先固定的粗组成层内做条件比较。若病理尚未审定，明确使用的是组成层和描述性检验，不把自动 cluster 当作独立生物学真值。

不要求每个 niche 都包含全部7个 dataset，也不把临床状态、患者差异、采样部位等自动解释为技术 batch。

### 7.2 稳定性

对最终入围的工作配置，而不是每一个参数组合，完成固定种子 124–128 的有限稳定性检查。

- 保持固定表示时，可按相同规则抽取约80%的 reference，重新拟合 KM 或重建 SNN/Leiden。
- 所有配置共用同一套重抽样 UID。
- 报告共同 anchors 上 ARI、匹配后的群重叠、群数及画像变化。
- 明确称为“条件于固定 PCA/Harmony 表示的聚类稳定性”，不能写成全流程或独立患者验证。
- 对 R12 三模式的工作结果补一个按切片内空间块、或有依据的独立分析单位分块的稳定性检查；说明分块方法和可比较集合。
- 若某 Harmony 分支存在具体的收敛/随机性疑虑，仅对对应 R12 分支补有限重拟合；不要无条件对所有参数重复全量 Harmony。

同来源行重抽样的高 ARI 不能替代跨来源支持。全量 Harmony 已使用所有 anchors，任何仅删除聚类 reference 的检查都不是完整的 leave-dataset-out integration。

### 7.3 冻结 DATASET 与 SAMPLE 各自的强度

对两个模式分别比较 WEAK/STRONG，结合：

- 可比较环境的局部混合是否改善；
- 原始组成区别和连续变化是否保留；
- source 特异群是否有可解释证据；
- 覆盖率和稳定性是否付出明显代价。

强度不能根据“恰好分出8群”决定，也不要求越强越好。若较强档没有清晰额外收益或出现难解释的组成合并，选择较保守档。若证据互相冲突，保守档可作为待审工作方案，明确“尚未证实优于 NONE”，继续完成并交付两档比较。

每个模式选出的强度对 R1/R2/R12 和 KM/Leiden 统一使用。保留两个整合模式，最终报告可以推荐其中一个，也可以保留两个有不同优缺点的方案。

### 7.4 小范围 q05 与 mean 来源结构诊断

如果现成、可追溯的 all30 mean 丰度缓存可用，在固定 R12 anchors、rings 和 reference 上做一次 q05/mean 的来源结构对照。仅用来评估输入摘要带来的差异，不替换 q05 主输入，不增加新的主配置，也不运行表达模型或后验重训练。

若需要重新生成大型对象或缺少可追溯缓存，将该诊断标记为未完成并说明，不阻塞18方案，不把它变成新的全项目任务。

## 8 全量标签应用和置信度

所有查询必须处在同一次、同一模式、同一输入的坐标系中。

本轮规定：

1. PCA 在固定 reference 拟合并投影全部共同 anchors。
2. Harmony 对全量投影联合拟合一次。
3. 从同一次返回的 embedding 提取固定 reference 行，拟合聚类。
4. 其他 anchors 也从这同一次返回的 embedding 取行后应用标签。

禁止用未校正 PCA query 查询 Harmony reference；禁止 reference/query 分别跑 Harmony 后拼接。R12 选中的 pilot 若输入、参数和 UID 完全相同，可直接复用，不必再重复拟合。

标签层次：

- reference_native_label：reference 的原生聚类标签。
- candidate_label：全量应用的建议标签。KM 用对应空间最近质心；Leiden 用该空间 reference 的15NN投票。
- common_transfer_label：两种算法均采用同一15NN规则的诊断标签，供比较传播机制。
- accepted_label：通过预先规定的置信度规则的保守应用标签；不通过者显示 UNASSIGNED。
- assignment_status、rejection_reason、label_origin：明确来源与拒识原因。

15NN 沿用等票、reference QC 排除自身、平票按平均距离再按稳定ID处理。顶票比例是近邻一致性，不是后验概率。

工作置信度规则沿用 v2 的已有标志作为起点：

- 对候选标签的邻居支持比例 <0.6，标记低支持。
- 最近 reference 距离超过该表示 reference 留自身外距离的 Q99，标记距离异常。
- KM 要计算“投给最近质心所给标签的票数比例”，不能拿另一个标签的最高票率为其背书。
- reference 原生标签始终保留；对其应用质量使用留自身外诊断。保守 accepted_label 可以按同一规则标为未分配，但不得借此重写 reference 原生聚类。
- 保留原生/候选标签及所有 flags，不通过删除困难点制造漂亮的结果。

每个配置都使用同一全量分母。主比较同时报告：共同 reference 原生标签、全量候选标签、保守标签覆盖率及各方案共同可分配子集。不得仅在各自最容易分配的不同子集上比较并宣布改善。

完整环敏感性优先复用固定标签，在108,345个完整环 anchors上重新汇总组成、来源和覆盖率；明确这是固定标签子集分析，不声称重新拟合验证。

## 9 niche 解释与可视化

### 9.1 环境画像和名称

无论用哪种表示聚类，所有组成画像都回到未经 Harmony 改写的 q1/q2。

每个 niche 输出：

- R1 和 R2 的30类组成；保留两环而不是只看用于聚类的一环。
- 按已冻结映射求和的大类组成，说明它是同一份点估计的汇总。
- 每个 dataset、sample 的原始画像、样本量和群占比。
- 相对全体上皮 anchor 环境的富集/偏低描述；零分母处理要透明。
- 有来源证据的患者数，以及其他分析单位的证据分层计数，分别报告。
- dataset-dominated、sample-dominated、small、redundant、uncertain 等工作标志。

名称先使用“上皮背景中 SPP1 TAM/中性粒相对偏高”“SMC/血管基质相对偏高”等由数据支持的描述。不得因某亚型占比略高，就自动命名为纯 CAF niche、TLS、免疫抑制 niche、肿瘤界面 niche或功能机制。

Fibro/CAF、T/NK、Myeloid 是用户关心的方向，但不能因为这些方向符合预期而人为挑选 K。没有清晰分离就如实报告。

可提供 niche 组成质心的层次树和宏观环境建议。若提出人工合并，保留原始标签、成员表、组成依据和合并映射；本轮不通过无记录的手动合并抹掉来源差异。

### 9.2 图表要求

按 dataset 组织输出，保留 section 层；小 dataset 不省略。

至少包括：

1. R12 五表示的 UMAP 来源对照，按 dataset、sample 及两种聚类分别着色。
2. R1/R2/R12 最终三模式的 UMAP 和来源诊断。
3. 所有主配置的双环组成热图、niche 大小与来源堆叠图。
4. dataset/sample 混合改善与原始组成保留的并列比较。
5. 逐切片空间图：对每个 R1/R2/R12，采用“三模式 × 两算法”的六面板比较；灰色明确显示未分配/不适用。
6. 原生标签与保守分配覆盖率的对照，避免只展示通过置信度筛选的漂亮区域。
7. 工作 niche 之间的原始组成相似性/层次树，以及同 K 算法比较。

同一 embedding 的不同着色共用 UMAP 坐标；不同 embedding 可分别运行相同超参数和 seed 的 UMAP，不能把二维轴直接解释为对应轴。UMAP 只展示，不用于建主图、聚类或选参。

固定每个配置的离散颜色和图例。不同配置的同一数字 cluster ID 不代表同一 niche；如按组成匹配颜色，输出映射和匹配不确定性。不要使用连续色条表示离散群，或重复少数颜色使不同群无法区分。

空间坐标保持正确纵横比和方向。透明底 PNG 和适宜的 PDF；文字与图例清晰。类别过多时按 dataset 分页，不能省略来源说明。连续 QC 图遵循现有配色偏好，零值与缺失值区分。

## 10 单独交付用户的区域复核任务

这项与 niche 计算并行，不是其前置批准。8257个自动证据连通区域不等于8257个独立病理病灶，也不要求用户全部填写。

复用133个优先ROI及现有索引，并按影响 anchors 数、证据冲突和可获得的图像排列。先整理以下区域的证据：

| 优先区域 | v2 涉及共同 anchors |
|---|---:|
| CanDis_2026_Sample8__uncertain_or_conflicting__0001 | 10,161 |
| GSE319536_BC11__uncertain_or_conflicting__0001 | 4,669 |
| CanDis_2026_Sample1__uncertain_or_conflicting__0001 | 4,291 |
| CanDis_2026_Sample3__uncertain_or_conflicting__0001 | 4,269 |
| GSE285715_T1__uncertain_or_conflicting__0001 | 3,913 |

对每个重点ROI尽量提供：

- 原始 H&E 全图、配准 spot 叠加、ROI 边界和局部放大图。
- 同一位置的 mean Epi gate 比例、CopyKAT 原始分类和已有连续 CNV score。
- 必要时从现有 profile 容器提取小型局部图，不运行新的 CopyKAT。
- 已有病理注释、原论文/原始数据来源和配准状态。

先检索服务器现有 H&E、scalefactors、tissue_positions 和病理标注。图像缺失时列出准确的 section、文件类型、已有位置或原始来源；不要生成假 H&E，也不要把普通散点图命名为 H&E。

CopyKAT aneuploid、diploid、unknown 使用固定离散图例；缺失与 diploid 分开。diploid/缺失 CNV 不能单独裁定正常，Epi gate 也不是肿瘤 mask。

大ROI包含多类组织时，支持按 polygon/spot membership 拆成子ROI，保留 parent_roi_id、稳定子ID、成员UID及版本；不要强迫整个大连通区域只有一个标签。

生成 WHAT_USER_NEEDS_TO_DO.md，开头直接告诉用户：

1. 本轮 niche 已完成到哪里，哪些 ROI 仍待审。
2. 优先看哪几个 ROI，每个对应的可点击证据图。
3. 每个 ROI 要回答的具体问题，例如“是否可判定为肿瘤区域、正常上皮、基质/其他组织，还是仍不确定？”
4. 到哪个表填写哪些列：reviewed_label、review_comment、reviewer、review_date；需要拆分时如何标注。
5. 无法判断时如何填写 uncertain，缺 H&E 时需要补的究竟是哪一张图。
6. 完成后把哪个文件交回；后续将影响区域分层/界面解释，本轮不会因此自动运行 MP10。

提供一份空白可填写的 ROI_review_to_fill.tsv。填写示例放在说明文档中并明确标为示例，不写入真实复核表冒充审定。不得填入虚构 reviewer、日期或病理确认。

本轮不必实现或执行最终肿瘤界面、深度和 MP10 推断模块。旧版 all_regions_reviewed 的全局条件不能约束 v3 niche 完成；也不能把旧 MP10 解盲开关改成 true 来规避它。

## 11 交付文件与完成状态

建议目录如下。可根据现有项目规范做等价调整，但最终报告给出真实路径。

~~~
neighborhood_ring_harmony_v3/
  00_preflight/
    RUN_SCOPE.json
    input_contract.json
    input_hashes.tsv
    batch_key_mapping.tsv
    coverage_by_dataset_sample.tsv
    software_and_resources.json
  01_features/
    feature_manifest.json
    reference_anchor_ids.tsv
    pca_manifest.json
  02_R12_pilot/
    NONE/
    DATASET_WEAK/
    DATASET_STRONG/
    SAMPLE_WEAK/
    SAMPLE_STRONG/
    pilot_comparison.tsv
    pilot_selection.md
    selected_harmony_parameters.json
  03_main/
    <18_config_ids>/
      working_config.json
      candidate_parameters.tsv
      reference_native_labels.*
      labels.parquet
      profiles_R1_R2.tsv
      profile_by_dataset_sample.tsv
      support_and_assignment.tsv
      stability.tsv
      graph_or_model_manifest.json
  04_comparison/
    eighteen_scheme_comparison.tsv
    integration_metrics.tsv
    pairwise_ARI_NMI.tsv
    sameK_diagnostics.tsv
    legacy_vs_NONE_PCA.tsv
    niche_dictionary.tsv
    dataset_sample_support.tsv
    complete_ring_sensitivity.tsv
    niche_correspondence.tsv
    figures/
    by_dataset/
  05_user_review/
    WHAT_USER_NEEDS_TO_DO.md
    ROI_priority_review.tsv
    ROI_review_to_fill.tsv
    image_availability.tsv
    evidence/
  reports/
    RUN_SUMMARY.md
    REVIEW_INDEX.md
    METHODS.md
    SCOPE_AND_ACCEPTANCE.md
    review.zip
  code/
    run_niche_only.*
    resume_niche_only.sh
    config.*
    environment_versions.*
  frozen_niche_manifest.json
~~~

所有配置的共同 UID 和 reference UID 必须相同；各种参考/应用/保守标签不得靠当前行号隐式对应。保存至少一种可复现模型或足够明确的质心/图/坐标与参数记录。

eighteen_scheme_comparison.tsv 至少含：

~~~
config_id, integration_mode, batch_key, theta, lambda,
ring_input, representation_id, n_pc, algorithm, k_or_resolution, n_clusters,
n_reference, n_all, fit_status, convergence_status,
stability_type, stability_ARI_median,
low_support_fraction, distance_outlier_fraction, accepted_fraction,
dataset_dominance_summary, sample_dominance_summary,
composition_preservation_summary, recommendation_status
~~~

数值指标与文字评价分开，不能把文字 PASS 当成实际度量。

frozen_niche_manifest 至少记录：

- 18个配置的状态与选参。
- 特征、PCA、reference、batch mapping、代码和标签哈希。
- 被选择与未选择的 Harmony 参数。
- 候选/原生/保守标签定义。
- 本轮推荐的工作 niche 方案及限制。
- 区域复核独立状态。
- stop_after=niche、mp10_authorized_this_run=false。
- 任何失败、偏离或未完成项。

review.zip 放报告、比较表、组成与支持表、重点图和用户复核材料。排除大型原始 counts、全量 c2l/CNV 对象、重复缓存和巨型 Harmony 对象。若过大，按 dataset 拆图包并写清索引；不通过丢失关键比较来压缩。完整标签和模型留服务器。

验收重点：

1. NONE、DATASET、SAMPLE 均实际计算；18个配置都有实际完成或明确失败状态。
2. q05环境、mean Epi gate、共同集合、reference、两环权重与类型顺序保持一致。
3. 两种算法在相同模式和输入下使用同一表示。
4. reference 与全量应用空间来自同一次整合；不做跨空间距离匹配。
5. 明确修订旧K选择规则，保存完整有限候选及工作选择理由。
6. 同时交付来源混合、原始组成保留、稳定性、来源支持和覆盖率。
7. ROI缺少人工结论不阻塞独立niche任务；需要用户做什么写得具体。
8. MP10及来源验证未执行，scope开关仍关闭。

配置失败时只定位对应问题，执行授权范围内的有限修复并继续独立分支。不能用复制另一配置、静默丢样本或全部强制赋群来填满18行。最终仍有失败时，交付已完成结果、真实原因和准确续跑命令，不能写“18方案全部通过”。

## 12 最终给用户的回复

完成后请用中文，按下列顺序回复：

1. 实际完成的样本/anchor范围，以及18方案完成数量。
2. NONE、DATASET、SAMPLE的实测比较：哪个改善了什么，代价是什么，哪些差异仍然存在。
3. 每个模式推荐的KM/Leiden粒度；是否比v2的K2和38–45群更可解释。
4. 主要niche的两环组成、跨dataset/sample支持；哪些是共享候选，哪些仍是来源特异或不确定。
5. 直接链接 RUN_SUMMARY、REVIEW_INDEX、18方案比较表、主要图和回传包。
6. 单独列“你现在需要做什么”：优先ROI、证据图、需填写列、缺失材料、如何交回；不要只说“请人工复核”。
7. 明确本轮止于niche；MP10分析和来源验证留待后续，不主动续跑。

若数据不支持唯一最优模式，可以保留两个工作方案并解释差别。不要只展示一个所谓“最优UMAP”就结束。

## 13 方法依据与项目自定部分

METHODS 中区分文献支持和本项目实施决定。以下链接供实现与方法描述使用，不要求因此引入其他整合工具或重做上游模型。

1. Harmony 原始论文：Korsunsky et al., Nature Methods, 2019, Fast, sensitive and accurate integration of single-cell data with Harmony. DOI: 10.1038/s41592-019-0619-0。
   https://www.nature.com/articles/s41592-019-0619-0

2. 聚合后组成整合的直接空间先例：Multimodal spatial profiling reveals distinct myeloid-defined ENKTL subgroups and tumor-myeloid cooperativity within prognostic inflammatory niche, Leukemia, 2026. DOI: 10.1038/s41375-026-03118-7。其 Methods 使用邻域细胞类型组成 → Hellinger → PCA → Harmony 跨患者整合；它是 CosMx 圆形邻域，不是本项目的 Visium c2l 双环，不据此补写未核实的最终聚类算法或参数。
   https://www.nature.com/articles/s41375-026-03118-7

3. BANKSY 官方 Visium Harmony 示例：空间邻域表示的 PCA 后运行 Harmony，再聚类；支持整合位置和矩阵接口，不要求本轮新增 BANKSY。
   https://prabhakarlab.github.io/Banksy/articles/batch-correction.html

4. Harmony 正式接口与参数：
   https://pati-ni.github.io/harmony/reference/RunHarmony.html
   https://pati-ni.github.io/harmony/reference/RunHarmony.default.html
   https://github.com/immunogenomics/harmony/blob/master/R/ui.R
   https://github.com/immunogenomics/harmony/blob/master/R/harmony_option.R

5. 整合评估：Luecken et al., Benchmarking atlas-level data integration in single-cell genomics, Nature Methods, 2022. DOI: 10.1038/s41592-021-01336-8。借鉴同时评估批次去除与生物信息保留的原则，不把缺乏独立真值的 niche 标签直接当作 benchmark ground truth。
   https://www.nature.com/articles/s41592-021-01336-8

本项目的 Epi gate、独立 R1/R2、双环等权、q05/mean分工、18配置、两档Harmony参数、粒度搜索与置信度工作阈值是预先规定的项目设计，不声称由某篇论文统一规定或验证。

