# 方法与可评价性

## 冻结参考、样本和空间评分

沿用 v1 阅读的 Fig. 1、3A/B/J、ST-seq Methods 和勘误。原研究对每片进行 BANKSY 聚类，再用 MIA 标记富集组合定义区室。肿瘤生态区为 TNN+TNS，TIL 生态区为 TIN+TIS。Fig. 3B 对比免疫治疗应答组；其差值、轮廓和带宽没有完整披露。本分析使用项目的 μm 和 200 μm 规则。[原文](https://doi.org/10.1186/s12943-026-02682-x) · [勘误](https://doi.org/10.1186/s12943-026-02777-5)

主标签直接读取 v1。其 BANKSY 设置为 λ=0.2、k_geom=6、逐片 2,000 HVG、最多 30 PC、Leiden 分辨率 0.7。每 cluster 的上调 marker 与参考 E/I/S 区室 marker 在共同可测基因中做单侧超几何检验，主阈值 p<0.05。八类是组合结果，cluster 编号只在本片解释。N 指该次检验未富集。Unknown 记录未分类。E 仅代表上皮参考富集，细胞比例、恶性身份和 MP10 各自独立。去 Top20 的 v1 标签是固定敏感性。已对齐的邻域复用 SAMPLE_R12_LEIDEN 的 accepted_label 0–8，仅在已赋值 anchor 内比较。[v1 方法](../mp10_visium_niche_boundary_v1/METHODS.md) · [30 类参考映射](audit/v1_reference_mapping.tsv) · [R12 字典](audit/R12_dictionary.tsv)

逐片、按完整 barcode 连接 H5AD 原始计数、v1 空间表、组织块、图像和 c2l posterior mean。重复 symbol 在完整基因轴上求和后排序。复用冻结 EpiMP10_Top20_AUCell，降序排名包括零，ties 取平均秩，每片 K=ceil(0.05×本片基因数)。对可测 Top20 成员累加 max(0,K−rank) 并除以相同成员数在最佳秩次上的总分。未在 E 域重筛基因。表达 QC 为已有 center_qc、UMI>0、nGene≥150。D_MP10 为 all_E ∩ 表达 QC ∩ 分数可计算。零分及 Top20 低检出均保留。抽取各片原 query 的至多 12 点，无原 query 时从需覆盖的有分点抽取。用原脚本相同函数复核。交叉表按同一 barcode 直接计数。旧显示条件为 query 且已评分；旧肿瘤主分析另有肿瘤来源、inside 限制。新 P95 使用同片全部 D_MP10 线性第95百分位、≥保留所有并列，正常的标签仅为本片 E 域描述。[冻结基因](audit/v1_top20.tsv) · [评分交叉表](audit/score_cross_cohort.tsv)

## 双边界和局部单位

全部 E 包含 ENN/ENS/EIN/EIS，第二套仅 ENN/ENS。使用组织块内观测 hex 一阶邻接，一对异侧 spot 形成一条有限共享边。Unknown、空洞及截断边不生成轮廓。保留岛和断裂。按组织块和当前同侧连通分量求点到线段最小欧氏距离，内部负、外部正。无本分量观测接触点时为 NA。strict 要求最近观测边距离严格小于最近缺测/截断边距离（1e−6 容差）。border_spot_distance 则为同侧分量内到最近接触 spot 中心的距离，中心可正负零，带别按显式 mask 判定。只重复 200 μm 的连续原始/调整终点。

主带内侧 −200≤d≤0 与 d<−200，外侧 0<d≤200、d>200，100/300 μm 另算。全部 E 外侧 MP10 域为空；ENN/ENS 外侧拆分 EIN/EIS 和 non-E，后者只作 TME 描述。每个观测位置面积为 √3×0.005 mm²。网格以阵列物理坐标原点按 floor(x/w)、floor(y/w) 划分，各组织块独立，按实测 QC 点数求面积。局部表保留 B1/B2 的 n、分数和、P95 个数、免疫丰度分子和总丰度。同片分量配对要求每个分量至少 3+3 个点，按 min(nB1,nB2) 加权，合计每带仍需 20 个点。分量形态表包括无深部的岛。分量外围补表将 G1 点仅分配到同块、同外侧分量最近共享边所属 E 岛，等距时沿用几何索引返回的一条边。窗口外围按位置只计一次。共同域 C = D_MP10 ∩ ENN/ENS ∩ 两套距离有效，每张片分别分带, 双方均过 20+20 才计算配对差。[边段与外侧类型](boundary/) · [片内效应和分母](tables/section_effects.tsv.gz)

## 免疫指标与盲态发现

各谱系丰度取相应名下 c2l posterior mean 之和；组成比例取谱系和/全部 30 类和，密度取谱系和/实测面积。全片、E 内、100/200/300 μm 内外过渡带与同块 non-E 外带分列。记录 E 内丰度占全片丰度份额及除以 E 面积份额的比。P_T = log2[(E 内 T 密度+1)/(外带 T 密度+1)]，伪计数单位 estimated cells/mm²。缺外围为 NA；外带仅有 1–4 个实测位置时，密度比另标为低支持。CD8 合并 exhausted 与 GZMK/effector-memory 两个 posterior，T 排除 NK。c2l 已有验证标记提示亚型可辨认度有限。另列两种 CD8 与 NK 后验相关绝对值均小于 0.90 的描述子集，这只是一项共线性筛查，不能证明细胞身份正确。[实际分组成员](audit/c2l_axes.tsv) · [参考验证意见](audit/c2l_reference_usage.tsv) · [逐片 CD8 审计](audit/CD8_identifiability.tsv)

W 为实际同块邻居的行标准化二元权重。R1+R2 去重并排除中心，同 spot 用单位权重单列。孤点不进入跨点指标。令 x、y 为相应谱系 spot fraction，共现为 mean(x·Wy)；双向 CLQ 为 [ΣxWy/Σx]/mean(y)。接触比 proxy 为 ΣxWE/ΣxWx。z 为在有效非孤点内的标准分数。全局 Moran 为 mean(z·Wz)，局部 T Moran 为 z_T·Wz_T，双变量 Moran 为 mean(z_x·Wz_y)，Lee L 为 mean(Wz_x·Wz_y)。常数标准化或零分母返回 NA。普通 Visium 的这些量衡量混合 spot 的邻接结构，不能定位单细胞接触。筛后局部分布报告 IQR、中位数和 MAD；F 与 H 使用空间支持的窗口标签，H 为两种已支持模式面积比例的 Shannon bits，与未决面积份额一起读取。[完整词典](audit/metric_dictionary_with_states.tsv)

发现只读取分开的免疫表、几何连接和来源。肿瘤片发现，正常仅投影。按预设优先级选供给、进入、共现、聚集、异质性等 8 个候选；缺失>20% 或与先保留者 |Spearman r|>0.85 时剔除。全片 T 密度、R1 E→T CLQ 和 500 μm 窗 E 内 T 比例 IQR 先 log1p，再以肿瘤样本均值和标准差标准化，无删异常值或零填充。前 80% 方差、最多 5 PC，固定 PC 空间比较 K=1–5 对角 GMM，协方差正则项 0.05，20 次初始化。另做 10 个不同初始化与 80 次按 analysis_unit 重采样；重采样只用有单位映射的标本，保留单位内所有切片，使用原 PCA 坐标。考察 BIC、silhouette、ARI、k-means、来源 AMI、QC η²。8 个特征另拟合一维 GMM，Scott 带宽 ×0.8/1/1.25，峰 prominence ≥最大密度 6%，60 次单位重采样。原 cluster 与解释标签分别保存。[模型选择](classification/model_selection.tsv) · [峰数与稳定性](classification/multimodality.tsv)

## 冻结规则与推断

500 μm 窗基本门槛 5 个 E 点、外围支持 4 点。进入样窗需 E 内 T 密度≥150、T fraction≥0.075、E 的 R1 邻居 T fraction≥0.08，且有外围时 P_T≥log2(2/3)。外围不足时须另有更强 E 内证据，密度≥250 且邻域≥0.10。低浸润样窗来自两种条件：E 内密度<150 且 T fraction<0.075，或有充足外围、外密度≥150、P_T<log2(2/3) 且内密度<250。出现两种低进入证据时优先记外围停留。1000 μm 所需 E/外围点数为 20/16。

窗口需同块有另一个同状态窗，通过实际 E–E 邻接跨窗相连。具此支持的面积须覆盖全片 E 面积的≥45%。在符合 E 点数要求的窗口 E 面积中，两种状态各≥20%且各≥3 窗时判 mixed；否则单种状态占≥50%且≥3 窗时判相应主导格局。1000 μm 最低 2 窗。证据不足记 uncertain，无可用窗记 not_evaluable。15/25% 仅变更 mixed 份额。阈值来自盲态跨切片局部四分位与物理可解释要求，具体依据写入 JSON。冻结标签哈希后才接入 MP10。[阈值与原因](provisional_classification.json) · [原始标签和局部面积](classification/section_classification.tsv)

每片 D = mean(B1)−mean(B2)；同片 OLS 的 B1 指示系数调整 E_fraction 和 log1p UMI，设计秩不足记 NA。P95 终点为两带概率差×100，按有效标签分母检查 20+20；去 Top20 的 P95 仅使用主分型域已有 P95 的交集。曲线 100 μm 箱每片≥5 点，在同侧连续箱内用 1:2:1 邻箱均值作显示，中心化减同片 D_MP10 均值。跨片按单位内等权、来源内单位等权、来源间等权，展示实际切片、单位和来源数。组间差仅在含双方的来源内先求，再跨来源汇总；一患者可贡献多个状态，来源内的对比影响量按同一单位合并计算。跨来源 95% t 区间以来源数−1 为自由度，至少 3 个来源才提供常规 P。稀疏来源保留点估计并列单臂少于 2 的标志。每片结果属描述，未用 pooled-spot P 推断。6 个主要成对连续比较共用 BH，缺组比较仍占预定族规模。每种片内终点及格局的总体均值另外按两个边界作 BH。其他终点、敏感性和组成分别成族，八种比例在同一相关组成族，Unknown 另作覆盖。连续关联按每个终点将全部候选指标作 BH；来源内至少 4 个有效单位才估 Spearman，95% 来源内区间为近似 Fisher 变换，跨来源区间是来源均值的 t 区间。关联是探索性统计，许多候选共享分子分母。

样本 high/low 在每片 D_MP10 中取 MP10 中位数，同一单位重复片等权，再依来源内单位中位数区分 high、low、tie，邻近正常不入组。排除 B1 的敏感性只取目标距离可确定且在 B1 以外的合格 E 点，要求≥20 点。各格局比例中同一单位若包含不同状态则按其切片份额平分一个单位的权重。H&E 提供形态参照。TLS、临床治疗、片内 E/I/S、c2l 和 MP10 的独立性含义各异；本分析未取得独立病理浸润标签。[逐来源文献与已知字段](audit/source_studies.tsv) · [逐标本核对](audit/section_clinical.tsv)


逐片临床核对中，GSE319536 按 BC 编号从论文补表 S1 填写取材后的方案、分期及病理 TLS。NPJ 全部为高等级 pT1，原文 Table 1 的具体膀胱部位未强行对应本地 sample1–7。其余来源按可验证的样本或研究层级记录。外部临床信息未进入 PCA、阈值或结果筛选。[MIBC 补表 S1](https://doi.org/10.1158/2159-8290.CD-26-0099) · [NPJ 原文](https://doi.org/10.1038/s41698-026-01417-x) · [Gouin](https://doi.org/10.1038/s41467-021-25103-7) · [Batf3](https://doi.org/10.1016/j.celrep.2024.114141) · [ZFP36](https://doi.org/10.1073/pnas.2505125123) · [SPP1/CXCL9](https://doi.org/10.1158/2159-8290.CD-25-1774)
