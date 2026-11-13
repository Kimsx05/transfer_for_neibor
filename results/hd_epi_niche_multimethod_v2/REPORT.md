# Visium HD 上皮相关 niche 四方法比较

2026-11-13T11:43:20.215900+08:00

已完成四方法的 24 个预设结果分支及 1 个条件敏感性，均输出五片全部冻结 anchors 的赋值或明确未定状态。主要工作层仍建议 A NONE K=2 + 连续全组成；CellCharter 作为分子辅助，B/MENDER 的细分保留诊断意义。

## 完成范围与固定分母

| section | anchors | niche_train | evaluable_min | evaluable_max | A_NONE_assigned_min | A_NONE_assigned_max | CC50_assigned | CC50_unresolved | CC200_assigned | CC200_unresolved |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 72171 | 20000 | 72171 | 72171 | 72171 | 72171 | 72155 | 16 | 72057 | 114 |
| D1_Johnson | 97858 | 20000 | 97857 | 97858 | 97857 | 97858 | 97732 | 126 | 97803 | 55 |
| UC_SCC | 263581 | 20000 | 263581 | 263581 | 263575 | 263581 | 263510 | 71 | 262054 | 1527 |
| MPBC1 | 77767 | 20000 | 77767 | 77767 | 77767 | 77767 | 77746 | 21 | 77735 | 32 |
| MPBC2 | 14175 | 0 | 14174 | 14174 | 14174 | 14174 | 14082 | 93 | 13939 | 236 |

evaluable_min/max 是全部25分支的可评价范围；A_NONE_assigned_min/max 是双参考×双尺度四组的已赋值范围。训练外数量=已赋值数−该片训练数；不同方法的具体 source 与未定原因见 coverage_wide.tsv 和逐中心压缩表。niche 描述 anchors/bin 环境，不是真实细胞数。

五片原始对象重读和挂接验证：True。非 anchors 保持 niche 缺失，原 counts/QC/坐标未修改。原525,552个anchor ID不变。

## 结果解释与工作建议

## 本轮结论与使用建议

24 个预先规定的工作/敏感性分支及 1 个由环支持量诊断触发的 MENDER 对照均已完成五片全量赋值，覆盖固定的 525,552 个 anchors。计算完成不等于每个选定分区都适合作为生物学 niche。保留 K-means NONE 的 K=2 作为主要粗结构，同时保留连续原始组成；CellCharter 提供有分子依据但仍受来源与特征选择影响的探索性细分。当前不建议把 adapted Leiden 的数百标签或 MENDER 200 μm 的退化二分作为正式细 niche。没有运行 MP10 关联或通讯分析。

### 1. 四种方法实际回答的问题

| 方法 | 实际实现与正式范围 | 解释地位 |
| --- | --- | --- |
| A | 完整组成 Hellinger + 全非零秩 PCA + MiniBatchKMeans；BLCA/HPCA × 50/200 μm × NONE/WEAK，共 8 分支 | NONE K=2 是可用的粗层；WEAK 是配对敏感性，不替换原始组成 |
| B | 旧 Squidpy raw-P 原生结果只读；新同表示图采用保留频数/边权的完全重复 profile 收缩 + constrained Leiden，共 8 分支 | 必要适配，不冒称 Squidpy 原生；重复输入问题已修复，细碎分群仍只作诊断 |
| MENDER | 25 μm 环宽的 2/8 环计数，官方 normalize_total/log1p/PCA 逻辑，适配图与外部投影，4 分支 | 保留径向证据；当前选定分群不作为通用工作 niche |
| CellCharter | 原始 counts → scVI reference/query → 中心外真实圆内 latent 均值 → 原生 GMM；50/200 μm 各 1 模型；另有固定 K 的特征排除对照，共 4 分支 | 50 μm K=5、200 μm K=9 为探索性分子环境，不能当纯细胞类型 |
| 条件敏感性 | 仅 MENDER BLCA50 增加每环已知类型比例归一化，固定原 resolution | 专门检验环支持量，不是额外生物学主分支 |

所有 niche 模型拟合均使用同一 T_v2=80,000；CellCharter 的 T_expr=200,000 是另一份全类型表达训练名单。训练中心之外：A 使用 predict，B/MENDER 使用明确记录的外部 kNN 投影，CellCharter 使用原生 predict/predict_proba。训练标签不被投影覆盖。全量表保留 candidate、label、source、confidence、distance、evaluable、reason；-1 表示模型未分配，不是注释 Unknown。

### 2. Squidpy 重复组成拆分：修复成立，但不等于细分成立

旧 BLCA50 的 4,608 个完全相同纯 Epi 训练输入确实分到 89 个原生标签；HPCA50 的 11,424 个输入分到 245 个标签。训练组成相似图存在大量断开的分量和零距离并列近邻，不能把这些标签差异解释成组织中的不同环境。

修复使用原训练图 A 及 profile 成员矩阵 G，Q=GᵀAG，组内边保留为相应自环，保留频数及边质量。小图验证了在相同 profile 必须同组的约束下模块度目标等价。全部 25 个最终分支的“完整实际输入相同却分到多个训练标签”均为 0；检查的是各自真实输入，而非仅旧累计 P。证据见 `identical_input_diagnostics.tsv`、`native_graph_diagnostics.tsv` 与 `identical_input_graph_diagnostics.tsv`。

但 B NONE 50 μm 仍有 BLCA 171、HPCA 651 群，WEAK 后分别 371、1,638 群；这说明重复拆分不是唯一的离散图问题。NONE 200 μm 各 6 群也与切片来源/检测量关联。BLCA200 WEAK 虽选到 2 群，却是 525,443 对 105 个已赋值中心，buffer ARI=0；不能把它称为稳定的二类生态环境。k=15/30/50 图敏感性与所有群完整保留，不删小群来改善外观。

### 3. K=2 的粗组成轴有支持，参考和尺度依赖仍明显

NONE A 的四组均仍选择 K=2。以下为全量已赋值 anchors 的平均原始 Epi 份额，并非细胞丰度或注释准确率：

| 参考/半径 | label 0 | label 1 | 解释 |
| --- | ---: | ---: | --- |
| BLCA50 | 0.893 | 0.381 | 较 Epi 占优 vs 较间质相关的复合环境 |
| BLCA200 | 0.793 | 0.366 | 较大范围的同类粗轴 |
| HPCA50 | 0.957 | 0.612 | HPCA 中两组均可 Epi 占优；不能沿用 BLCA 的纯类型名称 |
| HPCA200 | 0.851 | 0.539 | 较大范围的参考依赖粗轴 |

BLCA50 label1 的平均 Fibroblast+Mural 为 0.546，BLCA200 为 0.544；HPCA 不支持 Mural，不能将其当作零比较。A NONE 的条件 buffer ARI 约 0.82–0.99，支持粗分区的可重复性；不代表每个 anchor 的身份真实。

共同 T 上，A 的 50/200 μm ARI 为 BLCA 0.493、HPCA 0.379；同半径跨参考 ARI 为 0.508、0.537。局部与较大环境存在对应，但不是同一分区，也不是独立重复。B/MENDER 与 A 的低 ARI 不能自动称为发现了更丰富的生物结构。完整交叉表、原始组成及组内标准差保存在 tables/。

### 4. MENDER 的径向信息与缺测支持必须分开

50 μm 主计数分支得到 BLCA 32、HPCA 227 群。BLCA50 的纯 Epi labels 14/5/3/22 分别对应平均 input coverage 约 1.00/0.77/0.74/0.67；按环归一化后，这四组中的纯 Epi 主体合并到同一个 label22，合计 36,917 个中心。总体 32→27 群、同 T ARI=0.897，coverage η² 从 0.127 降至 0.012。这个具体拆分主要反映环支持量，而非已证明的生物差异。HPCA50 主分支 coverage η²=0.404，也需特别谨慎。

200 μm 两参考均选 K=2 式的图分区，但分别为 525,402/147 和 525,396/139 个已赋值中心，超过 99.97% 落在一个组。它们的条件稳定性 1.0 与跨参考一致性不能证明有两种稳健生态位。这是本轮有限流程选定候选的退化结果，不推广为 MENDER 方法普遍无效，也不继续为追求更多群无限调参。保留连续分环组成和标签用于诊断，暂不推荐这些离散标签进入正式关联。

### 5. CellCharter 增加了分子环境信息，但细分仍需限制用途

在完成 KL warm-up 修正并全部重拟合后，最终 50 μm 为 K=5，200 μm 为 K=9；此前短 warm-up 结果已明确失效，不能与最终表混用。两组的条件 buffer ARI 分别 0.972、0.681，原生 GMM 未分配中心分别 327、1,964。表中的原始邻居 marker 是真实邻域 bins 的面板检出比例，不是中心自身表达，也没有使用 BANKSY 平滑值。

50 μm 的具体分子证据：label0 和 label4 的邻居上皮面板检出分别为 0.835、0.814；但胶原面板为 0.086、0.469，髓系面板为 0.025、0.180。这支持“上皮信号明显、伴随环境不同”的描述。label2 的尿路上皮面板 0.231、基底/鳞状面板 0.007，与 label4 的 0.036、0.775 不同。不能据此宣称纯细胞群、EMT、恶性或真实混合细胞比例；面板与上游注释/模型特征相关，不是独立真值。

一个直接的新增信息例子：旧 BLCA50 中 4,608 个完全相同纯 Epi 组成的训练中心，在 CC50 中分为 label0=1,808、label2=2,066、label3=726、label4=8。这些中心的分子输入并不相同。label0/2 的基底面板检出为 0.288/0.002，尿路面板为 0.405/0.365；label2 中 2,014 个来自 MPBC1。因而这不是与 Squidpy 相同的“相同完整输入任意拆分”，但新增信息与来源仍纠缠。详见 same_P_marker_case_summaries.tsv，不能仅凭该细分宣称新的跨样本生态位。

200 μm 分子分区与来源的关系更强：section NMI=0.619，中心 nFeature η²=0.387，参考冲突 η²=0.430；50 μm 分别为 0.383、0.143、0.436。UC_SCC 的 74.68% anchors 进入 CC200 label0，D1 的 73.67% 进入 label1，MPBC1 的 74.31% 进入 label6。病理背景可能真实不同，但不能据此排除来源或局部技术质量的影响，更不能把这些差异当作跨患者重复验证。

冻结 Top20 清单中只有 MSMO1、ELOVL5、MXD1 进入主 HVG；移除后按相同 T_expr、同工作 K 重拟合 scVI/GMM。主/排除分支共同 T 的 ARI 为 50 μm 0.543、200 μm 0.575。该结果是特征/模型敏感性，变化也包含重拟合表示的影响；不能全部归因于三基因生物作用，更不能宣称证据已完全独立。未读取 MP10 分数、high/low 或任何下游关联结果。

### 6. WEAK Harmony 的收益、代价与 MPBC2 范围

A 的 section NMI 从 NONE 的 0.100/0.129 (BLCA50/200)、0.068/0.059 (HPCA50/200)，降至 WEAK 的 0.043/0.054、0.039/0.023。与此同时 NONE/WEAK 的原生分区 ARI 仅 0.704/0.691、0.830/0.787。来源关联降低伴随组成分区变化，不能把真实病理构成差异一并视为应去除的 batch。保留 NONE 为主、WEAK 为敏感性；不以跨片混合最大选优。实际 Harmony 达到 10 次迭代上限，日志记录未收敛停止，不声称已充分收敛，也不扩展强校正搜索。

MPBC2：NONE A/B/MENDER 不参加表示/模型拟合，仅外层预测；WEAK 在五片可评价 anchors 上共同学习 Harmony，因此是 transductive；CellCharter 使用官方 load_query_data 对 50,000 个 MPBC2 QC bins 适配后编码全部 422,741 个 QC bins，原四片 latent 数值保持检查最大差为 0。所有分支中 MPBC2 均未参加 niche 拟合、K/resolution 选择。后两种范围不能称为完全独立端到端外推。

### 7. 五片空间复核与推荐

实际打开并检查了四方法空间图、固定 ROI、已有原始单基因图及 D1 可用 H&E 图。不同模型数字 ID 无语义对应，颜色也不能作为跨方法同类证据；ROI 颜色键引用相应全片图的完整 legend-only PDF。

| 切片 | 空间与原始信号证据 | 本轮推荐及限制 |
| --- | --- | --- |
| Cervilla | 粗组成轴在上皮主体与周缘/复合区域间变化；far_UMI_hotspot 的局部粗轴变化在 200 μm 汇入较广的复合环境。CC 增加区域划分，但部分边界邻域检测量较低 | A NONE 双参考粗层 + 连续组成；CC50 作分子辅助。无 supplied 输入的图像区域不赋 niche |
| D1_Johnson | conflict_focus 的既有原始 KRT5/14/17/6A 与 COL1A1/A2 信号同时广泛可见，窄尿路型面板较弱。CC50 label4 和 CC200 label1 的上皮/胶原/髓系检出支持复合环境描述 | 不按 BLCA Fibroblast 或 HPCA Epi 单方定性；A 双参考并列，CC50 辅助复核；仍不能区分所有真实混合、参考适用性与局部测量影响 |
| UC_SCC | 既有 typical ROI 有基底/鳞状相关信号与胶原空间差异；CC50 的基底相关高/低胶原环境可解释一部分新增细分。200 μm 大片区域集中到同一来源相关组 | A NONE 粗层 + CC50 探索，CC200 仅作大环境敏感性；不是与其他片病理匹配的技术重复 |
| MPBC1 | 固定 typical ROI 中 anchors 稀疏分布，A 粗层与 CC 的尿路型相关环境可对应；B/M50 的小标签跳变没有充分证据被解释成多个真实 niche | 使用相对网格、连续组成与 A 粗层；表达辅助可保留，组织范围未知，不做病理边界或组织内比例推断 |
| MPBC2 | 固定 high_information ROI 可见 EPCAM/KRT7/KRT19/UPK2，基底面板弱；局部右上胶原信号保留。原 low_information、sparse_input_edge ROI 无 anchors，未替换为更好看的区域 | 只作可评价局部投影；组织范围未知，query/Harmony 参与范围明确；不作为独立验证队列 |

缺测、边缘、注释 Unknown、非 anchor 与 niche 未分配均保留。全量结果的高赋值率不是身份可信度；尤其 MPBC2 存在较多 feature-exclusion 投影未分配，需使用具体分支的 reason/support。来源、深度、参考适用性与组织构成仍不能完全分离。

### 8. 后续可用层级及普通 Visium 对应

可以进入下一轮方案讨论的是：A NONE K=2 + 原始连续全组成，明确参考/尺度/支持度，并以空间区块及切片作为相关结构处理；不把几十万个重叠 anchors 当独立患者。CellCharter50 可作为探索性分子环境，须同时查看 feature-exclusion 结果和来源/检测量；CC200、B/MENDER 的细分暂不支持稳健的细胞类型特异归因。下一轮即使做 MP10 共定位，也不能据当前结果预先宣称关联独立或因果；本轮只评价适用性，未执行关联。

普通 Visium 正式版本是包含 Epi 的全类型组成、Epi 中心、中心排除；R1/R2 分开特征拼接。HD 的累计 200 μm 圆、MENDER 的 25 μm 同心环、c2l 丰度和 HD 等权 bin 份额不是同一种测量。可在共同四粗类上比较环境方向，不能直接对齐 niche 编号、绝对比例或细胞数。只读对应表已交付，未重跑 63 张普通 Visium。

本轮没有必需输入阻塞，也没有待运行模型分支。仍未解决的是生物学可辨识性、参考依赖、MPBC 组织掩膜缺失、弱 Harmony 未充分收敛，以及表达型细分的特征/来源敏感性。方法不可用或只能探索是本轮的结果，不以增加群数作为完成标准。


## 实际运行矩阵与参数

| method | reference | radius | representation | status | niche_fit | representation_fit | branch_role |
| --- | --- | --- | --- | --- | --- | --- | --- |
| hellinger_kmeans | blca | 50 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| hellinger_kmeans | blca | 50 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| hellinger_kmeans | blca | 200 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| hellinger_kmeans | blca | 200 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| hellinger_kmeans | hpca | 50 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| hellinger_kmeans | hpca | 50 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| hellinger_kmeans | hpca | 200 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| hellinger_kmeans | hpca | 200 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| adapted_leiden | blca | 50 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| adapted_leiden | blca | 50 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| adapted_leiden | blca | 200 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| adapted_leiden | blca | 200 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| adapted_leiden | hpca | 50 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| adapted_leiden | hpca | 50 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| adapted_leiden | hpca | 200 | NONE | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| adapted_leiden | hpca | 200 | WEAK | full_assignments_complete | T_v2 | all_five_sections_transductive | predefined_working |
| MENDER_adapted | blca | 50 | ring_counts | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| MENDER_adapted | blca | 200 | ring_counts | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| MENDER_adapted | hpca | 50 | ring_counts | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| MENDER_adapted | hpca | 200 | ring_counts | full_assignments_complete | T_v2 | T_v2 | predefined_working |
| CellCharter | both_for_interpretation_only | 50 | scVI | full_assignments_complete | T_v2 | reference_T_expr_plus_MPBC2_query_adaptation | predefined_working |
| CellCharter | both_for_interpretation_only | 200 | scVI | full_assignments_complete | T_v2 | reference_T_expr_plus_MPBC2_query_adaptation | predefined_working |
| CellCharter | both_for_interpretation_only | 50 | scVI_without_Top20 | full_assignments_complete | T_v2 | reference_T_expr_plus_MPBC2_query_adaptation | predefined_working |
| CellCharter | both_for_interpretation_only | 200 | scVI_without_Top20 | full_assignments_complete | T_v2 | reference_T_expr_plus_MPBC2_query_adaptation | predefined_working |
| MENDER_adapted | blca | 50 | per_ring_normalized | full_assignments_complete | T_v2 | T_v2 | conditional_quality_sensitivity |

selected_parameter 是有限网格中的计算候选；其生物学使用资格由 model_recommendations.tsv 单独记录，不因入选就认可为工作 niche。

| tag | selected_parameter | n_clusters | n_train | evaluable | unresolved | seed_ARI | ARI_to_full_on_buffer_test |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CellCharter_both_r200_scVI | 9.0 | 9 | 80000 | 525551 | 1964 | NA | 0.6813 |
| CellCharter_both_r200_scVI_without_Top20 | 9.0 | 9 | 80000 | 525551 | 7848 | NA | 0.7406 |
| CellCharter_both_r50_scVI | 5.0 | 5 | 80000 | 525551 | 327 | NA | 0.9723 |
| CellCharter_both_r50_scVI_without_Top20 | 5.0 | 5 | 80000 | 525551 | 2135 | NA | 0.677 |
| MENDER_adapted_blca_r200_ring_counts | 0.005 | 2 | 80000 | 525551 | 3 | 1.0 | 1.0 |
| MENDER_adapted_blca_r50_per_ring_normalized | 0.005 | 27 | 80000 | 525551 | 13 | 0.9639 | 0.9207 |
| MENDER_adapted_blca_r50_ring_counts | 0.005 | 32 | 80000 | 525551 | 15 | 0.9449 | 0.8871 |
| MENDER_adapted_hpca_r200_ring_counts | 0.005 | 2 | 80000 | 525551 | 17 | 0.9622 | 1.0 |
| MENDER_adapted_hpca_r50_ring_counts | 0.1 | 227 | 80000 | 525550 | 351 | 0.8596 | 0.8999 |
| adapted_leiden_blca_r200_NONE | 0.01 | 6 | 80000 | 525551 | 50 | 0.9776 | 0.6729 |
| adapted_leiden_blca_r200_WEAK | 0.005 | 2 | 80000 | 525551 | 4 | 0.9622 | 0.0 |
| adapted_leiden_blca_r50_NONE | 0.005 | 171 | 80000 | 525551 | 87 | 0.9819 | 0.9224 |
| adapted_leiden_blca_r50_WEAK | 0.005 | 371 | 80000 | 525551 | 490 | 0.9806 | 0.9016 |
| adapted_leiden_hpca_r200_NONE | 0.02 | 6 | 80000 | 525551 | 420 | 0.8336 | 0.6837 |
| adapted_leiden_hpca_r200_WEAK | 0.005 | 20 | 80000 | 525551 | 39 | 0.9398 | 0.8465 |
| adapted_leiden_hpca_r50_NONE | 0.2 | 651 | 80000 | 525550 | 389 | 0.8838 | 0.9384 |
| adapted_leiden_hpca_r50_WEAK | 0.005 | 1638 | 80000 | 525550 | 4830 | 0.9657 | 0.8513 |
| hellinger_kmeans_blca_r200_NONE | 2.0 | 2 | 80000 | 525551 | 1 | 0.8957 | 0.8211 |
| hellinger_kmeans_blca_r200_WEAK | 2.0 | 2 | 80000 | 525551 | 1 | 0.9385 | 0.9617 |
| hellinger_kmeans_blca_r50_NONE | 2.0 | 2 | 80000 | 525551 | 3 | 0.9189 | 0.9915 |
| hellinger_kmeans_blca_r50_WEAK | 2.0 | 2 | 80000 | 525551 | 1 | 0.9693 | 0.9852 |
| hellinger_kmeans_hpca_r200_NONE | 2.0 | 2 | 80000 | 525551 | 1 | 0.9894 | 0.9548 |
| hellinger_kmeans_hpca_r200_WEAK | 2.0 | 2 | 80000 | 525551 | 1 | 0.9379 | 0.9715 |
| hellinger_kmeans_hpca_r50_NONE | 2.0 | 2 | 80000 | 525550 | 8 | 0.9701 | 0.9458 |
| hellinger_kmeans_hpca_r50_WEAK | 2.0 | 2 | 80000 | 525550 | 9 | 0.9891 | 0.9866 |

## 输入、实现与验证入口

METHODS.md 描述每种原生步骤及适配边界；config.yaml、versions.json、input_manifest.tsv 记录实际环境、来源和参数。旧METHODS的unpruned文字与代码不一致，本轮继承实际pruned字段。T_v2统一80000个，旧BLCA50仅1个不在共同可评价集合内的中心同片替换；training_overlap.tsv记录与四旧名单的差异。

implementation_checks.json、real_ROI_MENDER_equivalence.json、adapted_graph_pilot.json、final_validation.tsv、loader_validation.tsv 保存针对半径、中心排除、Unknown分母、等价图/环/latent聚合及最终ID挂接的核验。initial scVI 400-epoch warm-up不适配有限epochs的问题已修正，旧结果隔离在diagnostics/pre_warmup_correction，不属于正式结果。详见logs/warmup_correction.json。

重计算在23:00–07:00完成；白天仅做读对象验证、图例修正、解释和打包。原counts/Seurat/参考均未复制。绘图使用实际Arial；PNG300dpi，PDF文字/图例矢量、密集点层栅格。原始单基因/H&E图链接复用上游，不声称是新生成数据。

## 文件导航

| 文件/目录 | 内容 |
| --- | --- |
| tables/frozen_anchors.parquet, T_v2.tsv, T_expr.parquet | 中心、共同niche训练、表达训练名单 |
| all_anchor_assignments/ + tables/assignment_index.tsv | 25个逐中心结果分区；source明确训练/预测/投影/未定 |
| tables/coverage_wide.tsv, section_completion.tsv | 全部组合逐片覆盖与精简汇总 |
| models/ | 工作模型、必要预处理、训练表示及K2辅助层 |
| tables/candidate_metrics.tsv, common_K2_comparisons.tsv | 参数稳定性及共同K2对照 |
| tables/method_comparisons.tsv, partition_cross_tables.tsv | 同T与全量对应，附各自分母 |
| tables/niche_composition.tsv, MENDER_ring_composition.tsv | 全部群原始累计/径向组成，未删小群 |
| tables/quality_marker_support.tsv, technical_association.tsv | 原始邻居marker与质量/来源诊断 |
| tables/identical_input*.tsv, native_graph_diagnostics.tsv | 重复完整输入和图诊断 |
| tables/integration_scope.tsv, expression_fit_diagnostics.tsv | MPBC2参与范围、表达训练状态 |
| tables/fixed_ROI_index.tsv, raw_ROI_figure_sources.tsv | 固定ROI坐标及上游原始证据入口 |
| figures/ | 五片双参考四方法×双尺度、NONE/WEAK、ROI和诊断 |

图中每个模型的数字ID独立。全片图绘制全部中心和全部标签；群数>12时PNG键列出该片每面板最大4群，并明确提示完整legend-only PDF。组成热图最多显示30个最大群，覆盖率见figure_display_scope.tsv；全量数值无截断。ROI使用相应全片图色键，空ROI未换选。

## 读取与复跑

```r
root <- "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2"
ext <- file.path(dirname(root), "hd_bin8_banksy_extension_v1")
source(file.path(ext, "code/load_result.R"))
x <- load_extension(ext, "UC_SCC")
source(file.path(root, "code/load_result.R"))
x <- attach_hd_niches_v2(x, root, "UC_SCC") # x: 既有RNA对象，不覆盖任何旧列
```

```python
from importlib.util import spec_from_file_location, module_from_spec
root = '/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2'
spec = spec_from_file_location("niche_loader", root + "/code/load_result.py")
loader = module_from_spec(spec); spec.loader.exec_module(loader)
metadata = loader.load_niche_metadata(root, "UC_SCC")
adata = loader.attach_niches(adata, root, "UC_SCC")
```

loader按唯一bin ID挂接niche_v2_前缀字段；全片原对象的非anchors保持label缺失。全量压缩表本身的label=-1为未定，candidate仍保留。依赖上游项目原对象与本目录.venv，未另存25份表达对象。

`./rerun.sh` 核对缓存后按夜间窗口续跑；`./rerun.sh --report-only` 仅更新已有结果的报告和review，不启动拟合。review.zip不含完整表达、参考、全图邻接或训练缓存，未上传外部仓库。

新增占用（含隔离CUDA依赖）约 6.791 GiB；review.zip 约 55.9 MiB。
