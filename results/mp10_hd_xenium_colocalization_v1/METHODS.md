# 方法与复现条件

## 范围与身份体系

三片 HD、GSE326226 五个 Xenium 样本独立分析。GEO 系列记载五个肿瘤，缺少与单细胞参考或其他 HD 切片的患者级对照。样本名未用于推定配对。项目根未找到 AGENTS.md 或 README，采用实际挂载的 highres_VisiumHD/analysis 中报告、清单、成员、代码和日志以及 subtype/v2 的审阅版分类。SpatialBC 仓库提供 R 脚本，脚本引用外部 RDS；本地五份 H5 为表达输入，未取得可核对原始 cell_id 的作者最终标签。

32→30 沿用冻结映射，FOLR2 TAM 合入 C1QC TAM，Naive B 与 Memory B 合并。Pericyte/SMC 归 Mural。分类表同时保留合并理由和可检测 marker。祖先大类与未解析对象在不同层级矩阵中各占互斥类别；TME 重点谱系包括 T/NK、Myeloid、Fibroblast 和 Endothelial。对于其他谱系展示可靠大类。供体组成和样本抽样保存在参考成员表。身份程序无需 MP10 或空间邻近。候选、竞争标签、支持度及原始 marker 数分别保存。全片计数只对本平台分析单位有效。

## HD 表达与 RCTD

逐片从完整 supplied 8 µm 原始 UMI 重建 16 µm 矩阵，使用与上游区域一致的零原点、行列整除 2、原组织组件及格内连通成员。以条码和固定区域 ID 对齐；逐基因及全矩阵 UMI 守恒，逐区域 nFeature、成员数、质心及格坐标复核。完整成员即旧 T16_members 的全组织网格，该名称来自旧分析；建格不使用任何 T 资格筛选。原图右为 x、下为 y，图上保留 µm。8 µm BANKSY 与双参考 SingleR 用于组织解释。旧重建细胞没有转贴到新网格。D1 与 UC 有强角蛋白与基质混合背景，Cervilla 的上皮标记构成不同。观察窗口是提供组织坐标覆盖范围。

RCTD 采用本地 spacexr 2.2.1 的 full 模式，允许一个 bin 中多谱系表达贡献。重用冻存的 fine_reference.rds（30 类、原始计数、已登记训练成员；原对象剔除 TIGIT），逐片 fitBulk、choose_sigma_c，然后分 256 个 bin 的有界批次运行原包 fitPixels。所有批次共享该片校准。UMI_min=100，UMI_min_sigma=300，counts_MIN=10，sigma 校准随机种子 123、N_fit=200，其他参数继承 create.RCTD 默认。未通过算法过滤的格仍在全片表。原生权重分批保留，展示前按行归一化；大类权重为同大类子类和。竞争大类权重至少 0.25 且与第一名差小于 0.25 时标 mixed。大类候选权重至少 0.50、优势至少 0.12 并检出对应原始 marker 才作大类支持。细类限父类匹配、最高权重 ≥0.40、优势 ≥0.12、检出至少一个谱系内 marker 和 anchor。Treg 另需至少两个 Treg panel 成员或 FOXP3 检出。这些是展示用的操作规则，权重不表示实际细胞计数。原始权重、第二候选、覆盖分母和 major-only 全部可追溯。

Guo 等 Fig. 5k 用于 16 µm、参考注释及邻域分析的总体框架。本研究的 full、标签门槛、50 µm 与 1,000 次置换均为本项目配置。[Guo 等](https://doi.org/10.1038/s41467-024-54685-1) · [Seurat HD](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette)

## 原冻结 Epi64

原历史 Epi64 原样复用 measurement_eligible、MP10、H05、全部区域及 8 µm 成员。18,085 基因背景、平均并列秩 K905 及原始分母按上游对象；逐片复核 H05 等于合格区域原分数严格超过 linear Q95。圆按实际成员质心到邻居质心距离构建，200 µm 主、50 µm 补充。边同时符合原组件、观察窗口连通和直线穿越检查。将成员表通过原 8 µm 条码显式内连接；共享任何原条码即整对剔除。覆盖按同组件圆中 supplied 8 µm 格数×64/(πr²)，上限 1。每个区域占一个中心，与同片其他合格区域比较。剔除边数、剩余邻居和覆盖在 epi64_centers 中，覆盖 <0.5 或有效邻居为零记不可估计。旧 T32/T16 资格不套用于全亚群邻居。

## Xenium 注释与状态

仅取 feature_type=Gene Expression。逐片核验唯一 feature ID、符号和条码，逐细胞 Gene Expression 总数与 GEO transcript_counts 相等。质心/面积以原 CSV 为准，原多边形有独立 ID 核验。面积最小/最大 0.5% 标记复核，退化边界标未知；面积离群上皮不进入状态评分。每张切片另存混合谱系、面积及局部多边形复核。表达簇用归一化表达、900 HVG、15 PCs、邻居图与 resolution 0.45 建立，用于审阅同簇 marker 和身份覆盖。

新版转移参考使用项目最小基因协调的 scRNA raw counts H5AD，按冻结 30 类 ID 精确对齐。每细类最多 450 个细胞，上皮 1,200，样本轮流抽取。H05/nonH05 训练在全部原资格、冻结 MP10 和每样本原 H05 名单计数下保留真值；分组从完整可评分对象上恢复，逐样本核对 n_H05、原 Q95 和可定位的原正例，不在抽样或 panel 遮罩后重定阈值。nonH05 含全部其他合格上皮。另作状态训练每类上限 1,800、验证每样本最多 240；保存每个成员用途及来源。GSE326225 在主参考内。Seurat 在 4,966 个参考和共同 panel 的实际共有符号上重新 LogNormalize、选最多 1,500 个 HVG、ScaleData 和 PCA（最多 20 维），分别进行大类、谱系内细类及已支持上皮的状态转移。FindTransferAnchors 的 pcaproject 后用 TransferData 保留连续支持度，基因未测值不填零。分割细胞首先求大类。大类至少一枚原 marker、支持 ≥0.50、领先 ≥0.12。竞争得分 ≥0.25 且差 <0.25、两谱系各有原 marker 记 mixed。有大类支持但细类原 marker/anchor 或竞争差不足保留 major-only。细类得分 ≥0.50 且优势 ≥0.15，加与 HD 同源的 marker/anchor 规则。上述阈值不作为校准概率。跨类原 marker 覆盖和细胞深度影响保守性。[10x 标签转移](https://www.10xgenomics.com/analysis-guides/xenium-cell-type-annotation) · [GSE326226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326226) · [SpatialBC](https://github.com/VividKai/SpatialBC)

2026-11-25 从官方 MSigDB 下载并冻结 v2026.2.Hs 的 JSON，保留 SHA256。原生成员与真实 panel 的交集为 Reactome 10/27、Hallmark 35/74。主程序去除 SREBF1/2 后保留 8 个下游成员，保留 10 基因及 SREBF1/2 单基因作为对照。MP10_covered 为冻结 Top20 的精确交集 9 基因，ACAT2 保留原符号。三组读数均仅在合格上皮报告。UCell 使用 5,001 个五样本共同生物基因背景、maxRank=min(1500,N)、average ties、missing_genes=skip。零表达参与并列秩，超过 maxRank 的秩截断，按该版本 UCell 的 U 统计量归一化。含 TF 的补充读数和 8 基因派生读数有共享成员。H05_transfer 是从 scRNA 原 H05 真值转移所得连续支持度。逐样本按主读数严格 >type-7 Q95 标记 Chol-high，不拆 ties。至少 5 个超阈值、>2 个不同分数、最高分非零、超阈值中至少半数检出两枚以上 8 基因成员才生成分组。其他合格上皮写 other Epi。支持度与 RNA 程序均不测定胆固醇浓度、通量或供给。

[Reactome MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_CHOLESTEROL_BIOSYNTHESIS.html) · [Hallmark MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_CHOLESTEROL_HOMEOSTASIS.html) · [Reactome R-HSA-191273](https://reactome.org/content/detail/R-HSA-191273)

## 遮罩及降深验证

验证按整数据来源留出，训练排除同来源所有样本，并核对训练与测试样本及已记录患者 ID 无交集。grouped fold 的测试按完整来源内各样本抽取，保留近似自然阳性比例。遮罩后在参考实际有的 4,966 个基因背景重算两种实测分数，另从其余来源重新学习 H05 转移。每个验证细胞再以固定随机种子抽取一个真实 Xenium 样本及其中一枚合格上皮的检测深度，超出目标时用无放回多元超几何降采样，浅于目标的输入保留。分别输出 Pearson、Spearman、AUROC、average precision、阳性率和全 PR 曲线。不同来源阳性分母与 AP 需一起看。降深只能近似稀疏性。真实 Xenium 检查单基因主导、多基因共同检出、总 transcript、nFeature、面积与分数关系。跨方法分歧细胞和具有至少 5 个合格上皮的 80 µm 固定格保留实测计数及状态秩差。

## 空间统计

逐切片构建稀疏物理半径图，无自环。HD 窗口沿完整 supplied 8 µm 组织格及原组件。原阵列与 µm 验证仿射变换残差 <0.1 格。Xenium 将原分割多边形绘制到 8 µm 格，2 格半径方形闭运算容许普通细胞间隙，由该窗口的八邻连通组件界定组织片段。它是近似观察窗。两端同组件，沿线每半格检查是否落在观察窗，拒绝横越明显组织缺口的边。每张图固定后再筛选分析标签。Xenium 50/200 µm，同片窗口面积覆盖为圆形卷积；覆盖 <0.5 不进入中心比较。HD 仅作一次同物理尺度敏感性，用 Cervilla 中心坐标中位数最近点为中心取 1,024 µm 方窗，8 µm 原 BLCA 大类及同一大类投票聚合到 16 µm，并列 Mixed。

Squidpy 1.8.4 gr.nhood_enrichment，n_perms=1000、rng=123、normalization=none、handle_nan=keep，单线程内核。在相同调用中截取 native permutation moments 计算随机期望和总体 SD，与返回的 Z-score 核对。观察边按有向邻接槽计；对称矩阵的两个方向不能相加后称为无向计数。无代表节点的类型不可估，观察/期望分母为零时 NA。邻域图展示的对角线留白，完整数值可查。所有邻域组成以全部独立邻居为分母；同时给出已知/有细类支持数、每平方毫米观察面积密度和平均邻居总数。各中心等权，报告 high 减 other 的片内差异。全局类型置换仅用于描述性富集，未用全局置换测试给定上皮状态的特异性。无患者级 P 或独立性假设。不同平台仅比共同关系方向与评估覆盖。

## 运行与绘图

本轮可在白天执行。每个入口绑到至多 16 个逻辑 CPU，BLAS、OpenMP、Numba 及 data.table 单线程。各样本分别运行，原始表达稀疏，仅按需保留计算缓存。图例颜色随 shared/colors.json 固定。小提琴为 jins-routine-plotting 的原生 R 模板与导出器。空间图、原多边形、散点平滑和 Z-score 热图无对应模板，用相同白底、Arial 字体及独立 PNG/PDF 自定义绘制。

存档成员复核从 results/mp10_immune_endothelial_mechanisms_v1/tables/pseudobulk_membership.tsv.gz 直接恢复 116,448 个 Epi_all 和 5,841 个 H05 条码。InhouseData::T4295396 的已归档 H05 有 34 个，早期分位表记 33 个，唯一差异在 Q95 的序列化精度内。以显式成员为准；本轮身份参考、状态训练、验证抽样中的全部条码逐一核对，差异条码未进入状态训练或验证。表见 shared/H05_membership_check.tsv。

Xenium 各取 64 个合格上皮，将全部 5,001 基因的原始表达按平均并列秩独立重算四套 UCell 读数，逐值与实际 R 输出核对。误差见 reports/UCell_native_check.tsv。R 运行日志中有 stack imbalance 警告，保留原日志；独立数值核验与逐 ID 图/矩阵回读通过。

补充审阅在未定大类且高参考支持的候选中统计另一组既有谱系状态 marker，供复核使用。原条件下的 Unknown 均保留。用相同 H05 上皮参考、panel 特征工程和锚点转移参考来源标签，按细胞的 H05 转移和 8 基因程序各自的样本内 20/80% 秩分组检查来源权重与原始上皮标志。该诊断不能确定真实 Xenium 细胞对应的参考患者。

补充的 Epi64 邻居模型读数对拟合通过且与中心表达来源独立的邻居取归一化 RCTD 权重均值。区域覆盖和有效拟合数另存；这组连续值表示表达混合成分。
原 Xenium 边界多边形的拓扑、面积与质心距离另行审核；外边界简化可能导致自交。保留现有分割，在敏感性分析中排除拓扑异常的两侧端点并重新计算保留上皮的 Q95。使用原稀疏图，详见 reports/boundary_review.tsv 和 boundary_sensitivity.tsv。

D1 和 UC 在空闲配额释放后并行处理尚未计算的批次。辅助进程复用该片已选定的 sigma 和原包 Q 表；启动时在同样 counts 和同样 reference 上重建 bulk 平台校准，先对一整个已完成批次的 256 行全部权重核对（最大绝对误差须 <1e-8），之后才写入剩余批次。每个输出先写独立临时文件再原子重命名。批次位置或是否调用辅助进程不会重新学习分类规则。
