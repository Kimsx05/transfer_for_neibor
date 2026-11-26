# 冻结定义和方法

本轮的分析单位是切片。固定四片全片用于空间映射与定位，片内主要统计沿用旧 pilot 主表 query 的完整 section_id 和 barcode。历史资格来自 posterior mean Epithelial / all30 ≥0.25、既有 QC 与原评分可用性；保留所有实测 MP10 零值。检查输入 manifest、空间计数列之和、独立参考身份、原 niche 及两模型的逐键连接。项目实际根目录无顶层 AGENTS.md 或 README，沿用本地 pilot 的 README、METHODS、REPORT、manifest 和代码。新结果写到独立 results 目录。

## 方案 B

先从原 abundance_mean.csv.gz 读入 30 列绝对丰度，以其和核对旧 Epi fraction。主丰度为 log1p(posterior mean absolute abundance)。大类在原始丰度上求和再 log1p。敏感性分别使用原始绝对丰度、占 all30 均值总量的份额、占非上皮 29 类均值总量的份额。中心 MP10 是既有 EpiMP10_Top20_AUCell，来源为混合 spot 的实测 RNA。已有 cell2location 输出只读。q05 只用于还原原外围环境，单独标记。

## 方案 A

参考结构借鉴 Chu 等人的两次 SCT 标签转移。论文的空间共定位方法写 Pearson，作者 v1.0 代码第 144–149 行分别求 anchors 和支持度，第 205 行指定 Spearman。本轮使用本项目的参考、冻结 query 与注释，Pearson 为主、Spearman 为敏感性。[论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12284189/) · [作者归档代码](https://github.com/Cthfudan/hcc/blob/v1.0/Fig2%26S4-S5/spatial_seurat.R)

上皮参考逐 ID 复用旧 pilot，核对原完整状态池的 H05 布尔值。原状态定义在每个合格 specimen 的全体 Epi 上按 type-7 Q95 严格大于，other Epi 含全部其余原合格细胞。TME 参考从同一 37 个未留出 specimen 的原候选细胞读取既有最终细注释，按 specimen × 细亚群最多 35 个不放回，种子 20261125。两个固定参考共 2,449 个上皮细胞与 19,269 个 TME 细胞。保留 CNP0000460::P01 和 InhouseData::T4463209 完整来源留出，已知同患者不同取材的 a/b 同留出。抽样后的类别构成见 reference_sampling.tsv；预测支持度受这些构成影响。

读取 minimally harmonized 单细胞整数 counts 和各片全基因 Visium counts，重复空间 symbol 先求和。每个参考和每片 query 分别运行 SCT v2，最多用 1,000 个观测拟合模型，参考常规选择 1,500 个可变基因，conserve.memory=TRUE、do.correct.umi=FALSE，采用已安装环境的原生拟合。仅以候选中实测共有特征计算 30 个 PC；舍弃可能在 Seurat 自动改名时发生混淆的含下划线或竖线的候选特征，原库总量不变。FindTransferAnchors 使用 SCT、pcaproject、参考 PCA、参考模型重算 query 残差，k.anchor=5、k.score=30、k.filter=NA、Annoy 50 棵树。TransferData 使用投影的 30 个 PC 和 k.weight=50。固定映射种子 20261126，各实际 PCA/anchor 特征写入 features。MP10 成员参与常规竞争，不追加。H05_separate_transfer_score 是二类上皮参考的连续支持度。两类上皮分数之和为一，二者的负相关由构造决定；TME 参考的 29 类亦竞争单位总支持度。共同大类的映射分数为对应细亚群分数之和，保留全参考分母。低 TME 诊断以全片或 query 各自非上皮 29 类原始均值丰度之和的 Q25 为界，保留并列，所有诊断均不改变主总体。

## 实测功能

优先评分 Reactome 胆固醇生物合成（R-HSA-191273）和 Hallmark 的胆固醇稳态、脂肪酸代谢。实际成员来自项目已冻结 MSigDB 2026.1.Hs 的 gene_sets.tsv.gz，版本与每个成员均存档。补充独立描述项目内 SpaCET 1.4.0 随包 GOBP.gmt 的 CHOLESTEROL_EFFLUX。每个集合计算原版和删去固定 MP10 Top20 后的版本，各过程单列。若少于 5 个成员可测则预设为信息不足，同时列实际检测数。[Reactome](https://reactome.org/content/detail/R-HSA-191273) · [Hallmark 稳态](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_CHOLESTEROL_HOMEOSTASIS.html) · [Hallmark 脂肪酸代谢](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_FATTY_ACID_METABOLISM.html)

沿用本项目 MP10_spatial_score_pilot_core19_v1 的确定性 normalized AUCell 秩积分口径。背景为该片完整 counts 基因列表，K = ceiling(0.05 × 基因数)，每个 spot 以原始 RNA 整数倒序排名，并列取平均。分子为各可测成员的 max(0,K−rank) 之和，分母为同数目成员处于最前排序时的积分。稀疏实现仅需对非零项排序，零表达成员的全背景平均名次在 K 之外，贡献为零。每片抽取固定首、中、末 spot 对照全基因原公式。此评分无随机步骤，不使用 Seurat 预测表达、HVG 子集或参考转移的功能读数。原 MP10 沿用已有数值。去交集仅重算功能模块，参考映射不因功能集合重建。覆盖表另列模块和两种实际映射特征的交集。

## 统计和空间尺度

每种方法的完整方阵保存每对读数的有效 n、Pearson、Spearman 和两种残差相关。方阵以不重复的上三角保存，图省略对角线。直接比较在每片每类同时具有 H05、MP10、TME mapping score、c2l 均值丰度、Epi fraction 和深度的原 query 上完成。分别对相关两端用截距、Epi fraction 和 log1p(nCount_RNA) 做 OLS，再计算残差 Pearson；残差 Spearman 指线性残差化后重新排序求相关。常数和缺失相关记 NA。保留原始和调整值，无独立 spot P 值或显著性标记。跨片仅汇报方向与效应范围。

按细胞类型对比四种组合。固定另一端在两个可能的口径，各自求替换某端导致的相关变化，两次变化绝对值取平均，跨类型再取中位数，用于描述更敏感的一端。两个端点的平均有符号变化之和恰等于 A−B。该运算描述方法差异，无法分离测量误差或因果效应。功能评分另直接与两套 TME 读数计算相关。每片 query 的 H05 或 MP10 分数 ≥Q75 为 high、≤Q25 为 low，保留并列，若上下阈值重合，仅报告连续结果和未分组分布。输出组大小、中位数、均值与差值；两个相反分位尾部另连接组织位置、原始单基因 UMI、深度与功能得分。

原 niche 使用 SAMPLE_R12_LEIDEN 的 accepted_label（0–8），展示编号加一，颜色沿用既有色板。以原 spot_uid 和 section/barcode 连接，缺失或负值单列，未分配不计入各 niche 分布。原 R1 为一跳六邻接；R2 为恰二跳、剔除中心和 R1，原纳入至少分别 4 和 8 个有效邻居。R12 对原 q05/all30 邻居比例均值分别开方再拼接除以 √2，继承冻结标签，不重聚类。为区分尺度，在相同已分配中心上比较中心实测读数与中心 q05/all30、原 R1 或 R2 比例的相关。该匹配尺度对照沿用原 q05 空间环境；方案 B 主表使用 posterior mean 绝对丰度。原外围画像只解释定义。

## 绘图与计算

每片正式输出独立矢量 PDF 与 300 dpi PNG。高低区域使用 jins-routine-plotting 的小提琴内嵌箱线模板，其他图种沿用旧 pilot 的 Arial、viridis 和白底风格。相关热图统一 −1 至 1，灰色表示缺失或省略的对角线。散点显示全部有效点，趋势线固定高斯带宽为 x 全距的 8%。预设 Treg、SPP1 TAM、C1QC TAM 与 Matrix_CAF 作空间示例，各自保留两种 TME 读数，其余亚群在矩阵中完整展示。niche 汇总图使用各读数在片内已分配 query 的平均百分位，原始读数和 n 见统计表。原组织图仅缩放背景分辨率，保留原像素坐标。所有中间计数为短期稀疏文件；数值工作使用单线程 BLAS、CPU 亲和性合计至多 4 个核心。GNU time 记录各阶段时间和峰值 RSS。


各支持度空间图色标从 0 到本片该类实际最大值（全零时显示 0–1），跨片绝对强度须读色标。低 TME 首位亚群按最大分数统计，精确并列时用表中固定次序。

运行记录中，首次标签轴检查发现 R 的分组变量同名问题，修正后以独立的 barcode 顺序断言重新运行。正式计算的八组预测和特征表均写完；该进程在尾部读取脚本时另报了解析错误。保存的脚本已通过语法解析，八组输出均重新读入按原轴复核，见 logs/map_verify.log。原生 SCT 对部分基因的离散度估计有迭代上限提示，见 logs/map_warnings.txt。分数轴、非负性、类别和与动态范围的实测核验见 tables/score_QC.tsv；常数列保留为不可估计。
