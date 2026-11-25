# 四片普通 Visium 的 H05 标签映射方法

## 冻结范围与复核

切片固定为 GSE171351_B1、GSE238145_B17、GSE319536_BC01、GSE285715_T1。项目挂载根目录未见 AGENTS.md 或顶层 README。依据本地 full_v1、p95_front、immune/endothelial、section_specificity 和 scope_frozen_resume 的 README、方法、manifest、脚本及日志定位输入。空间读数取 full_v1 的 EpiMP10_Top20_AUCell，含混合 spot 的 RNA 信号。全片作定位，主要比较沿用已存 query。其上游先定义 posterior mean Epithelial/all30 ≥ 0.25 的历史资格，再使用原评分可用性与既有 QC。未使用 q05、R12、common anchor 或新支持度重选总体。完整 barcode 后缀原样保留。每片检查原始计数列和总 UMI、独立上游坐标/评分表、旧 query 与 30 列 mean 比例。所有输入以绝对路径记在 input_manifest.tsv。

单细胞标签从 40 个原合格 specimen 的完整上皮状态池读取。原方法在每个完整 specimen 计算 type 7 Q95，严格大于为 H05。其余全部原合格上皮为 other Epi。原始状态池与 Epi_all 以 specimen 加 barcode 核对。T4295396 一枚处在序列化阈值边界的细胞在后续 pseudobulk H05 表中被纳入，原始完整池为 FALSE。本轮使用原始标记；差异可查 membership_version_difference.tsv。原合格上皮外的细胞仅在已有非上皮大类注释时以其原大类进入候选。

## 参考和映射

从上述 40 个 specimen 的已注释细胞出发。预先留出 CNP0000460::P01 整个来源组及 InhouseData::T4463209 整个来源组。后者 a/b 是已确认同患者重复样本；前者依旧来源清单保守合组，其患者身份未经独立核实。剩余 37 个 specimen 内按 specimen × 标签、每层最多 35 个不放回，NumPy default_rng(20261125)，固定 9,957 个参考细胞。抽样权重改变各类先验，不能将支持度当作细胞数、比例或校准概率。真实切片、去 Top20 和伪 spot 全部使用同一参考 ID。参考取 minimally_harmonized 原始 RNA 整数 counts，根据 obs dataset、orig.ident 和细胞 ID 核对；空间取原始基因乘 spot 稀疏计数，重复 symbol 求和。未读 Harmony、插补表达或 K32 训练标签。

使用已安装的 Seurat。参考和 query 各自按 RNA library 总量 LogNormalize 至 10,000。先在四片共有且至少 5 个参考细胞检出的基因中，用参考 VST 选 1,200 个候选，追加可测 Top20；各 query 再与非零基因取交集。含下划线或竖线的特征保留计数用于库深，但不进入相似性，避免 Seurat 的自动改名。ScaleData、25 PC、FindTransferAnchors(pcaproject, k.anchor=5, k.score=30, k.filter=NA, Annoy 50 trees, L2 默认)，TransferData(pcaproject, dims=1:25, k.weight=50)。每次固定随机种子。逐特征核对实际 PCA loadings 和 anchor.features。主分析显式容许全部可测定义成员，便于考察共享 RNA 所产生的一致性。消融从原始参考和 query 同时去掉 20 行再归一化，保留其余同一候选特征，重新 scaling、PCA、anchors 与预测。未用消融重算 MP10 或细胞真值。逐 spot 检查十类支持度非负且和为一。

## 逐片比较与定位

H05_support 为高状态类的原始支持度。Epi_support 为两类 Epi 的支持度之和。H05_conditional_support 为二者比值，分母为零记 NA，分母 < 0.10 记低分母。主统计保留低分母。各切片、各读数单独在冻结 query 内计算 Pearson、Spearman、缺失、零分母和低分母覆盖。在每一对观测的同一完整集合中，分别将 MP10 和映射读数对截距、c2l Epi fraction 和 log1p(nCount_RNA) 拟合 OLS，再对两组残差计算 Pearson 和 Spearman。此处残差 Spearman 是先线性残差化再排序。表中另列与单个协变量的相关及仅调该项的残差 Pearson，不作 spot 独立性 P 值或空间 CI。去 Top20 的成对效应比较使用两个分支共同可用的 spots。

展示规则在计算效应前冻结。每片 query 内 MP10 ≥ Q75 且原始支持度 ≤ Q25，或反向 MP10 ≤ Q25 且支持度 ≥ Q75，保留并列。用原 array 六邻接和同一 tissue_component 合并相连 spots，输出位置、区域大小和原始 UMI。分位数均为线性分位数。原始支持度 ≥ Q95 另存描述性排序标志，不能据此给 H05 阳性率。全片非上皮诊断预定为 c2l Epi fraction < 0.10 且优势大类为其他谱系，汇总支持度、低分母及 marker 原始读数。另列原 0.25 门槛以下且其他大类占优的位置，作为 B1 缺少 < 0.10 区时的补充描述，四片统一保留。区位沿用正式 inside/outside 标签。低 Epi 比例与正常来源均不能证明胆固醇生物学为阴性。

## 来源留出的有限伪 spot 检查

复用旧 TEST 的按深度稀释 counts 和完整组成细胞 ID。仅取上述留出组中有原冻结状态定义的 P01T、T4463209a/b。完整连接属于同 specimen 的原成员并要求成员数等于 N_total，同一 mixture 的两深度逐个核对细胞多重集合。每 specimen 固定随机最多 80 个合格 medium mixture，配同一 mixture 的 low 稀释，保留旧场景构成。真实 H05 个数与占全部组成细胞的份额只由原始 H05 布尔值计算。逐来源报告真零分布、支持度对含量的秩相关、配对降深度变化。伪 spot 库限于旧共享基因。最终配对深度分析在两组 query 共同检出的相同 859 个映射特征上运行，两次固定相同种子及参考 PCA 输入，只分别构建 query 投影、anchors 和预测。初次按深度各自检测筛选的过程在 map.log 中，最终固定特征的两次运行见 pseudo_paired.log。重复细胞及旧场景选择影响可推广性，这属于旧 TEST 的重新标记检查。

## 绘图与资源

jins-routine-plotting 没有 Visium 组织叠加或相关散点模板。沿用项目 p95_front 与 section_specificity 的 Matplotlib 样式、真实 Arial、白底、viridis、独立矢量 PDF 和 300 dpi PNG，有图例另存 panel-only 和 legend-only PDF。空间坐标以原组织像素对齐，组织图为观测 spot 范围加边距的同一视野；原评分仅着色冻结 query。散点显示全部有效点，红线是固定带宽（片内 MP10 全距的 8%）高斯核加权平均。没有对效应范围裁切。顺序执行阶段，BLAS/OpenMP 以单线程为主，运行入口限制最多四个 CPU，不升级共享环境。仅 PCA/选定特征和小型诊断为稠密数据。运行资源见 logs/*.resources。

Seurat 方法入口见 [空间整合示例](https://satijalab.org/seurat/articles/spatial_vignette.html) 与 [TransferData 文档](https://satijalab.org/seurat/reference/transferdata)。


H05 原始支持度图每片用 0 到该片实测最大值，便于辨认 B1 的小动态范围；跨片绝对强度须读色标。条件支持度和 Epi 支持度统一使用 0–1。另附全片低分母位置单图。
