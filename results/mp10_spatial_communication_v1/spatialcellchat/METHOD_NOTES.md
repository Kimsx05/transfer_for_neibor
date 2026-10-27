# SpatialCellChat 原生位置推断与软归属

固定官方 commit `be88f300464cf970b36a6165ea7d2e41a47a7511`，独立包版本 0.1.0。项目 Rlib 隔离安装，不改旧 CellChat 1.6.1。MERINGUE 依赖来源官方 c0d28fec4659e891dd81db620f58318536ce2d58；为避开大型示例数据，仅安装其正式 R/src 代码。RcppExports 由官方源码的 Rcpp attributes 生成，未改算法。

每片完整背景独立推断；实测 counts 按全基因 library size 归一化至 10000 后 log1p（与共享输入导出层同口径）；不调用会覆盖表达的 ALRA preProcessing。CellChatDB.human 显式保留 Secreted Signaling、ECM-Receptor、Cell-Cell Contact 三类 protein LR。原生 MERINGUE 空间变化筛选 adjusted P<0.05，至少一个 LR 成分空间变化，复合物全部亚基可测。无手选通路或 LR。

位置模型使用原生 computeCommunProb，mass action Hill Kh=0.5、n=1，AGONIST/ANTAGONIST 默认保留，距离倒数衰减。整片只有一个 bookkeeping idents 标签，Treg 等所有类型的组成均在后续 Q 聚合中保留；不依赖 dominant type 或 hard-label min.cells。坐标按平台 array geometry 标定 um。预设范围：分泌和 ECM 为 2×片内 median nearest-neighbour distance，Cell-Cell Contact 为 1.05×dNN；tol=1e-6 um（仅避免恰好范围边界的格点被浮点误排除；dNN 取8位小数），scale.distance=1.01/min nearest-neighbour distance。该 commit 实际源码将 ECM 与 Secreted 一同作 diffusible；首片3×dNN敏感性因此同时放宽分泌与ECM，接触范围不变。

位置支持采用原生 filterProbability(nboot=100,seed.use=666,thresh=.05)：源码从100个发送位置的现有邻接概率采样估计95%分位阈值，并非完整置换患者差异检验。然后原生 filterCommunication(min.cells=NULL,min.links=10,min.cells.sr=10) 作位置链接、发送与接收位置覆盖筛选。稀疏显式零先清理并同步 net$tmp$prob.cell 和 net$prob.cell。未进行 group-level permutation，不制造 group P。

冻结位置矩阵后保存原生概率（含同spot），聚合时复制并去除同spot对角项。每个LR用公开 helper `computeAvgCommunProb_LR_Visium_Avg`，F=Q' P_cross Q，K=I(Q>0)' I(P_cross>0) I(Q>0)，B=F/K，无支持为0。主 total_mean_strength=sum_LR B，不做 group P 二次筛选；类型无组成覆盖为NE。major/fine分别从同一P聚合，绝不相加fine均值冒充major。保留逐LR F/K/B、总flow、同spot排除贡献、类型支持位置和有效丰度。源码helper已与小矩阵显式公式核对，正式每LR再次做公式一致性断言。

`total_mean_strength` 不是每个真实细胞平均。Visium 结果是细胞类型比例加权的跨spot通讯归属，依赖去卷积和混合表达；同spot排除仅为聚合选择，非重新求解禁止同位置通讯；邻spot contact类别不证明真实膜接触。不同算法单位/归一化不互比。总强度不代表免疫抑制强度。

表中 source/target_effective_abundance 为 `sum(Q)`（归一化组成权重的spot等效质量），不是后验细胞数；source/target_support_spots 为 `Q>0` 位置数，也不等于细胞数量或已确认患者数。Treg 即使从未dominant仍可有非零Q并参与全背景聚合。

数值检查 `numerical_checks.tsv` 验证每层级全类型 flow 总和等于对应冻结位置网络跨spot总flow（因为Q行和为1），同spot排除量亦守恒。major和fine主均值可以不同，这是不同K分母的预期结果。

特征筛选采用该版本 MERINGUE 的默认 Delaunay 邻接（`getSpatialNeighbors(filterDist=NA)`），单侧更强空间自相关的 Moran 检验后对本片 signaling genes 作 BH 调整。这只是进入通讯模型的空间变化特征筛选，不是 Epi H/L 差异检验。

可核查源码：[SpatialCellChat 固定版本](https://github.com/jinworks/SpatialCellChat/tree/be88f300464cf970b36a6165ea7d2e41a47a7511)、[位置模型与 Visium helper](https://github.com/jinworks/SpatialCellChat/blob/be88f300464cf970b36a6165ea7d2e41a47a7511/R/modeling.R)、[MERINGUE 固定版本](https://github.com/JEFworks-Lab/MERINGUE/tree/c0d28fec4659e891dd81db620f58318536ce2d58)。
