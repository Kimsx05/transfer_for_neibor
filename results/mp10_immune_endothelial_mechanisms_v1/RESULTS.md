# MP10 免疫与内皮机制分析

## 研究范围与有效覆盖

本轮使用最终 c2l 30 类的单细胞 reference membership、原 MP10 Top20 与冻结的 MP10_cell。原始范围为 40 个合格肿瘤 specimen，含 116,448 个 Epi；未使用 CellChat capped sampling、空间 spot 分数或 LLPS/C5-CAF 定义。全转录组 RNA 的 counts/data 共 90,011 个 features，原对象就地读取。内皮抽检原归一化表达与 `log1p(UMI / 全 RNA library × 10000)` 完全一致；完整分子、分母与 barcode 核验见审计表。

T4463209a/b 是已确认同一患者的两次取材；主分析固定 a，替换分析使用 b。源表的 provisional ID 状态原样保留，用户确认的重复关系另列在 `patient_relationship_audit.tsv`；全部模型均核验未同时纳入 a/b。下面的 n 为标本及记录/派生 patient unit 数，除已确认重复关系外，不声称所有 patient ID 的生物学独立性均已验证。

| 主模型 (每 profile >=50 cells) | 标本 / 记录 patient units | 数据集 | 状态 |
|---|---:|---:|---|
| Epi H50 vs L50 | 39 / 39 | 5 | completed |
| Epi H10 vs L50 | 33 / 33 | 5 | completed |
| Epi H05 vs L50 | 24 / 24 | 5 | completed |
| Treg，Epi MP10 中位数 | 32 / 32 | 5 | completed |
| Tfh/CXCL13-like，中位数 | 23 / 23 | 5 | completed |
| Capillary EC，中位数 | 12 / 12 | 4 | completed |
| Arterial EC，中位数 | 4 / 4 | 1 | not_evaluable：残差自由度 2，低于冻结要求 3 |

>=30 的主要可用范围为 H50/H10/H05：39/36/30，Treg/Tfh/Capillary EC：33/25/14。Arterial EC 在 >=30 时有 5 个 profile、2 个数据集，但调整 dataset 后残差自由度仍为 2，因而该敏感性也不可评估；不能只看有 5 个标本。81 个预设 DE 模型中 68 个完成，13 个明确不可评估。详见 `tables/model_status.tsv`。三种 High 的共同 >=50 标本集为 24；共同 >=30 集为 30，与各方案自身合格集分开分析。

所有已完成设计均满秩。主 receiver 暴露与 dataset 的 R² 分别为 Treg 0.164、Tfh 0.092、Capillary EC 0.130，均保留数据集内变异；这不消除未测量混杂。暴露效应按原始 main-a 标本的固定 SD 表达，中位数 1 SD=0.07852，Q95 1 SD=0.14747；两个暴露分开建模。

## 发送端：NECTIN2 上调较稳定，但不是尾部独占

配对 edgeR 模型已经控制 specimen，完整基因结果和联合主分析 q 均保留。逐一去掉单个 specimen 后，NECTIN2 的平均配对 logCPM 差仍为正 (H50 0.270–0.350，H10 0.427–0.510，H05 0.252–0.329)；这是一项描述性影响检查，不是把留一法结果当作独立重复。重点分子的主结果如下 (log2 FC，括号内为各 comparison 的 BH FDR)：

| 分子 | H50 | H10 | H05 |
|---|---:|---:|---:|
| NECTIN2 | +0.326 (0.000268) | +0.470 (0.000644) | +0.290 (0.0168) |
| VEGFA | +0.191 (0.142) | +0.302 (0.0915) | +0.204 (0.194) |
| PVR | -0.069 (0.624) | -0.209 (0.0998) | -0.282 (0.0920) |
| ANGPT2 | -0.299 (0.0247) | -0.428 (0.0171) | -0.691 (0.000462) |

NECTIN2 在三个主定义中均为正向且通过局部 FDR；H50/H10 的 5 次 LODO 均保持正向且 FDR<0.05。H05 的 LODO 仍均为正向，但去掉 HRA003620 或 InhouseData 后不再达到 0.05。a/b 替换保持方向；>=30 的 H05 效应为 +0.475，FDR=0.00126。

在同一批 24 个共同合格标本中，NECTIN2 的 H50/H10/H05 log2 FC 分别为 +0.252/+0.277/+0.290，FDR 分别为 0.000695/0.0125/0.0168。各自队列中 H10 的更大估计不能全部解释成“10% 是最佳活性阈值”；共同队列显示较平缓的正向趋势。数据支持相对高 Epi 状态与 NECTIN2 表达相关，尚不支持把 H05 定义为独占信号的真实活跃群。

既有 Epi cluster 注释覆盖 101,224/116,448 个细胞 (86.9%)，其余 15,224 个缺少该套旧 cluster 标签；没有为它们新造亚型。亚群内检查只使用有注释且占原 Epi >=5% 的主要 cluster，RSS 中的 Unavailable existing cluster 是缺失标签集合，不是新细胞类型。已有 Epi 表达 cluster 内的方向检查也提供部分支持：H50 中多个主要 cluster 的 NECTIN2 配对变化为正且通过该候选/cluster 家族 FDR；H10/H05 的 cluster 内样本数明显减少，许多结果不可评估或不显著。原 cluster 只是既有表达分群，不是新建恶性标签。

候选并未限于 NECTIN2/VEGFA。完整配体与膜分子结果中，AREG、HBEGF 在 H50 上调 (各自 q≈0.00589/0.00414)，EREG 在 H50/H10 上调 (q≈0.00326/0.0273)，但三者 H05 均不显著；MDK、MIF 在三个主对照中均未达到 FDR<0.05。这些发送端结果全部保留，不能据此直接建立某个 receiver 的作用机制，也不能把不显著等同于无表达或效应为零。

不能忽略质量关联：逐标本 MP10 与 log(nCount) 的 Spearman 中位数约 0.420，与 nFeature 约 0.455。TMM、配对设计与 cluster 内检查减轻部分解释混杂，但不能证明完全分离了测序深度、细胞状态与亚群构成。完整质量差异、单 specimen 影响和各数据集方向均保留；未新造 CNV/恶性标签。

## 免疫轴：受体覆盖支持可接触性，实际下游机制仍有限

人类 NicheNet v2 prior 中，NECTIN2→TIGIT、PVRIG、CD226 三个分支均存在，NECTIN2 也在 ligand–target 矩阵中。PVRL2 是 NECTIN2 的旧名；不能把 literal PVRL2 列缺失说成整个 NECTIN2 轴缺失，也没有补零或重复计数。

主合格标本中受体阳性细胞比例的逐标本中位数：

| Receiver | TIGIT | PVRIG | CD226 |
|---|---:|---:|---:|
| Treg (n=32) | 64.2% | 0.47% | 12.9% |
| Tfh/CXCL13-like (n=23) | 62.4% | 0.45% | 40.4% |

TIGIT/CD226 有 RNA 覆盖；PVRIG 低于预设 5% 和 10% 规则。低检出不是功能不存在，也不能用这些比例推断受体占有率或竞争胜负。Tfh 的 TIGIT 高不等于耗竭。

Treg 主 DE 只有 2 个负关联基因达到局部 FDR：LINC01871、ABI3，均 q≈0.0472；Tfh 有 2 个正关联基因 HEXD-IT1、AC069277.1，均 q≈0.0164。四个 receiver 主检验家族的联合 q 中，没有基因达到 0.05。完整正负结果仍保留，不能把同一 TME 复制为 High/Low 两份来增加“响应”。

Hallmark 层面存在协调变化。Treg 的氧化磷酸化负向 (NES=-2.738) 与 G2M 正向 (NES=+2.398) 在 5 次实际 LODO 及 Q95 中均保持方向，LODO 的相应局部 q 均<0.05。Treg 的 TNF/NFKB、Tfh 的 TNF/NFKB 和氧化磷酸化更依赖数据集：部分 LODO 方向翻转，Tfh 的 Q95 主敏感性也不显著。这些是控制设计后的基因排名富集，不能等同于 NECTIN2 引起的患者级因果反应。

CollecTRI、DoRothEA 和 PROGENy 的主结果、每个 profile 的内部 ULM 结果及样本层模型 P/q 分开保存。CollecTRI/PROGENy 的主 receiver 样本层检验没有 FDR<0.05 的结果，不能把内部富集 P 当作患者间关联 P。

## 内皮轴：有状态变化，未确立 VEGF 驱动的血管生成

Capillary EC 中 FLT1/KDR/NRP1/NRP2 的中位阳性比例分别约 87.2%/58.0%/74.0%/30.0%；Arterial EC 的 FLT1/KDR 也可检出，但主表达关联模型因覆盖不足不可评估。FLT1_KDR 复合体与 NRP1/NRP2 cofactor 定义单独交付，未把复合体名称当作一个实测基因。

Capillary EC 无单基因主 FDR<0.05，但 Hallmark MYC targets、TNF/NFKB、氧化磷酸化、mTORC1、UPR 等为正向。MYC V1 和 TNF/NFKB 在可运行的实际 LODO 与 Q95 中保持正向显著；剔除 PRJNA662018 后模型因残差自由度不足而 NE。该数据集贡献了 12 个主合格标本中的 7 个，其余为 CNP0000460 的 3 个、GSE326225 和 InhouseData 各 1 个；因此仍不能排除主要数据集驱动，不能把其余三次 LODO 当作四个数据集的独立验证。氧化磷酸化部分 LODO 失去显著性。

这不是“血管生成增强”的充分证据。VEGFA 自身配对 DE 未达 FDR<0.05，PROGENy 主 receiver 通路没有显著样本关联；Hallmark ANGIOGENESIS 主模型 NES=-1.359、q=0.202，Q95 敏感性反而为负向 NES=-1.761、q=0.0165。ANGPT2 在 Epi 的三个对照中下调，也应与 VEGFA 等其他候选并列解释。

丰度单独建模后，主要 lineage 分母的四个 receiver 均未通过主 BH q<0.05。Capillary EC 的点估计偏负，但不能将该趋势、RNA 程序和 TF 活性合并成一个“活化”标签。

## TF 与 pySCENIC

Epi 的 CollecTRI 显示 SREBF1/SREBF2 活性与 High 状态强相关。去掉 Top20 后效应减小但仍显著：H50/H10/H05 中 SREBF2 的 M20-masked 效应分别约 0.728/1.441/1.894，q≈1.47e-6/1.29e-9/5.60e-9。SREBF1 也保持正向。这说明相关调控程序超出了直接评分基因的贡献，但仍不是完全独立于状态定义的生物学验证。

不能直接据此画 SREBF2→NECTIN2 的确定因果箭头。具体 TF→候选配体 prior 和共同下游 targets 分表保存。CollecTRI 中 NECTIN2 的已记录调控者包括 CREB1、FOS、JUN、SP1、SP3；M20-masked 的 H50 中 SP1/SP3 活动为正且局部 q≈0.047/0.0127，但 H10/H05 不显著；CREB1/FOS/JUN 的对应主活动检验也未通过 FDR。因此暂不能把 SREBF 程序与 NECTIN2 上调连接为一条已闭合调控链。

建网基因在冻结 motif 排名库中的映射率为 Epi 14,551/20,216 (72.0%)、T 谱系 11,235/13,806 (81.4%)、内皮 12,415/15,590 (79.6%)。重点配体/受体及 SREBF1/2 均在参考中；其他未映射基因与不满足模块覆盖要求的情况单列，未得到 regulon 不作为反证。

内皮共享网络经 motif 筛选得到 196 个 regulon，7,376 个全量细胞已完成 AUCell。Capillary EC 主模型无 regulon 达到样本层 FDR<0.05 (最低 q≈0.112)。预先按 CollecTRI 选出的内皮 5 个 TF 没有对应 motif-pruned regulon，第二种子虽然完成了预设目标推断，其关键 TF 连接稳定性仍应标 not_evaluable，不能声称已验证稳定。

Epi motif 网络得到 166 个 regulon，T 谱系得到 236 个。Epi 的 SREBF1 regulon 含 21 个 targets，其中 ACSL1、SCD、MSMO1 为 Top20 成员，masked 后保留 18 个。SREBF2 虽在 motif 数据库中，但未获得 Epi motif-pruned regulon，应记为该网络不可评估，而不是与 CollecTRI 相反。T 谱系的 44,448 个全量细胞已完成 AUCell。Treg/Tfh 各有 228 个可检验 regulon，主模型均无样本层 FDR<0.05 (最低 q 分别约 0.392/0.977)，与 CollecTRI、DoRothEA 主 receiver 检验均不显著的结果一致。Treg 在剔除 InhouseData 的单次 LODO 中出现 ZNF226 信号，其余门槛、Q95、患者替换和 LODO 均没有显著 regulon；该孤立结果不能代替主结论。T 谱系预选 TF 没有 motif-pruned regulon，第二种子关键连接一致性也为 not_evaluable。Epi 的 116,448 个全量细胞已完成 Full/M20-masked AUCell，各主模型有 160 个可检验 regulon。Full 的 H50/H10/H05 分别有 93/98/98 个局部 FDR<0.05；masked 分别为 95/102/97 个，正负方向均交付。它们属于同一表达数据上的程序关联，不是独立的机制验证。

SREBF1 的 pySCENIC 结果显示直接评分重叠贡献很大：

| 对照 | Full AUCell 差 / q | M20-masked AUCell 差 / q |
|---|---|---|
| H50 | +0.01611 / 2.35e-5 | +0.000516 / 0.391 |
| H10 | +0.03831 / 2.08e-6 | +0.003515 / 0.0124 |
| H05 | +0.04647 / 2.50e-5 | +0.004821 / 0.0296 |

DoRothEA A/B 与 CollecTRI 的 SREBF1/2 在三个 masked 对照中均保持正向且 q<0.05；pySCENIC 的 masked SREBF1 只在 H10/H05 的各自主队列达到该门槛。masked SREBF1 的 LODO 方向均为正，但 H10 剔除 HRA003620/InhouseData、H05 剔除 HRA003620/PRJNA662018 后 q>0.05；共同 24 标本中的 H10 也为 q≈0.0628。因此支持部分程序一致性，不能称所有网络和所有敏感性均稳健。

第二种子中，SREBF1 的 20 条可比较非自 GRN 连接均恢复，两个种子的正权重 Spearman 相关约 0.797；KLF13 恢复 15/16，正权重相关约 0.607，ZDHHC24 连接未恢复。regulon 中的 TF 自身成员用于评分，不能当作推断出的自调控边，已从恢复率分母和网络箭头排除。其余预选 Epi TF 缺少可比较 motif 连接，仍标不可评估。连接恢复只表示该推断的技术一致性，不代表直接结合或患者层因果效应。

## NicheNet 两个方向

**Epi→TME**：目标始终来自各 receiver 自身的连续 MP10 关联。Treg 只有 2 个可映射负关联 targets，属于低信息量；Tfh 的 2 个局部显著非编码基因没有进入该 ligand–target prior；Capillary EC 为零 target；Arterial EC 模型 NE。没有放宽 FDR、借用身份 marker 或 Epi DEG 来补齐 receiver targets。因此目前没有可靠的正向 NicheNet 活动证据证明 NECTIN2 或 VEGF 的机制。预设的门槛、Q95、患者替换和 LODO 检查也没有解决该信息量限制：各个可运行 receiver 设置的可映射联合 targets 均仅为 0–2 个。

稀疏 Treg targets 下，native NicheNet 对常数 prior 预测的并列处理可产生较高 AUPR。原生值和排序原样保存，并标记常数预测、低信息量；它们不是候选获胜证据。另提供明确命名的非恒定 prior 诊断排序，不替换原生结果。

**TME→Epi**：target/background 同时移除 Top20 后，H50/H10/H05 的正向 target 数分别为 1,531/916/527，能够进行完整探索。主 TME-focused 正向排名如下：

| 对照 | M20-masked 前三 |
|---|---|
| H50 | TNFSF10、NUCB2、SERPINE1 |
| H10 | TNFSF10、NUCB2、MPZL3 |
| H05 | TNFSF10、NUCB2、MPZL3 |

TNFSF10 在 Full 与 masked 三个主对照中均居首，但所有敏感性并不完全相同：masked 的 H50 排名范围 1–3，H10 为 1，H05 为 1–24。NUCB2 的 H05 排名范围为 2–60。不同设置不是独立重复；完整 sender-agnostic、focused、阴性目标、联合目标及真实 TME 来源均交付。

这些排名是可能解释 Epi 状态的上游线索，不证明这些 TME 细胞在“High 微环境”上调了配体。非经典候选还需要核查蛋白定位、分泌/接触形式与实际受体路径，不能只按 prior 排名挑功能实验。

## 建议与未闭合证据

NECTIN2 是目前较适合优先继续核验的发送端分子：配对表达与多项敏感性支持较完整。后续应先确认 Epi 蛋白表达、与 Treg/Tfh 的实际接触及受体分支，再开展能够区分 TIGIT/PVRIG/CD226 的功能验证。当前数据不足以宣称 TIGIT 主导或免疫抑制已经成立。

内皮方向应保留 VEGFA/VEGFB/ANGPT2 等比较，但暂不把 VEGFA–FLT1 定为已验证血管生成轴；优先结合 HD 邻近关系、KDR 与实际内皮响应。反向 TNFSF10/NUCB2 等可作为等待空间和路径核查的探索候选。HD 支持全部标 pending；本任务没有读取尚未完成的空间共定位结论。

没有冻结新的 MP10 生物活性阈值，没有从多种设置中挑最小 P，也没有将 top5% 自动命名为真实活跃细胞。

## 交付与运行限制

完整表、图、直接绘图数据、方法、输入/版本/哈希与仅从表重画的入口见 README、FIGURE_INDEX 和 module_status.tsv。核心 651 项一致性检查和 37 项 SCENIC 检查全部通过。50 张图共交付 198 个文件，已完成 PDF 渲染、矢量/字体/页边界、PNG 300 dpi 检查及六张总览图和重点图的人工视觉检查；有图例的图同时保留完整、主体和图例 PDF。可运行模块均已完成，未将低覆盖或低信息量写成阴性机制结论。

实际运行中发现服务器全局 Rprofile 会将 shell 设置的单线程覆盖为 24。早期部分 R 阶段因此未满足要求的总线程上限；问题已记录，后续使用 `--vanilla`、运行时 BLAS/OMP/RcppParallel 控制及单核亲和性修复。没有调整或终止 HD 作业，也没有因此修改统计设计或筛选结果。该执行偏差见 `logs/resource_constraint_incident.json`。
