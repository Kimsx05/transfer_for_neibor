# 亚群距离增量 v2

## 冻结输入与目标
使用 full_v1 的 63 片、完整 barcode、query、坐标、组织块、评分与身份。60 片肿瘤研究来源和3片邻近正常分别汇总。主分数为原 Top20 AUCell，成员及与官方集合的交集在 gene_overlap.tsv。复用 MSigDB 2026.1.Hs 的全部55项 UCell，未改变评分或成员。

目标为 Matrix_CAF、iCAF、rCAF、vCAF、C1QC TAM、SPP1 TAM、Treg，以及合并 CAF/TAM。以所有30项 posterior mean 为分母，在 full_v1 的有效组织集合内取 P95/P90，>=阈值且>0，保留并列。目标全片统一定义。平均丰度中位数<1为描述性标记，不删点。各单因子的边际 q05 单列，合并目标没有把边际 q05 相加。P80历史表原样复制。

最近目标中心距离使用 nominal Visium 100 μm 点间距的已有 x_um/y_um，同组织块内精确近邻。同spot为0，无同块目标为NA；精确外缘距离不作前提。新增距离对全组织点计算，统计主集合仍严格取冻结query。谱系内比例单独保存，谱系分母<1另作标记；它们未用于目标选择。

## 局部界面
固定保守规则在程序关联计算前设定。上皮侧要求 Epithelial/all30≥0.5，原 epithelial 标记分数大于 fibroblast 标记分数，且至少满足作者 Tumor Cell/Urothelium、旧 Tumor/Normal urothelium 身份、原至少2个上皮标记被检出之一。间质侧要求 CAF/all30≥0.2，fibroblast 标记分数较高，且满足作者 Tumor Stroma/Normal Stroma、旧 Non-tumor stromal region、至少2个成纤维标记检出之一。相反的明确作者身份作否决，QC不合格不作任何侧别。规则没有使用MP10、EMT或其他待检程序。

采用原array六邻接和v3.4组织块，只有观测到的异侧邻接才生成界面。共享hex边以两点中点为中心，垂直于连线，长100/√3 μm。距离取到有限线段的精确最短距离，确定上皮侧为正、间质侧为负，未知侧保留unsigned距离和side_uncertain。每块独立，无界面不回填；接触点为一圈，下一层真实邻接为两圈。组织切缘、缺点或空洞不生成异侧接触边。直线距离可穿过同块内部非观测区域，未宣称是可通行路径距离或严格外缘距离。完整v3.4几何状态继续保留。原始作者区分使病灶内部 Tumor Cell–Tumor Stroma 接触可被识别。H&E图中青色为局部候选接触，橙色为观测病灶整体外缘。此规则为组成及标记支持的候选几何，未进行独立病理复核。

## 模型与空间相关性
每个section、每个轴独立计算原始Pearson、Spearman以及调整相关。协变量为上皮组成、log1p(nCount_RNA)、lesion_location哑变量。对MP10、程序与距离采用相同完整配对集合；两模型评分层仅使用原epi_score_valid且两层均有效的交集。至少20对才能估计，本轮比full_v1的主相关显示门槛保守。

交互模型一次固定为标准化Program ~ 截距 + composition + depth + lesion_location + 标准化MP10 + 标准化distance + MP10×distance。距离主项采用线性；因此交互只检验线性斜率调制。逐点数值标准化仅用于模型，原评分保留。每个距离轴单独建模，恒定/秩亏导致不可估计时保留NA。实际距离分箱曲线显示均值化的片内中位程序和该箱实际query距离中位数，供观察非线性；曲线不构成额外非线性检验。

正式片内不确定性使用按组织块及固定500 μm方格划分的空间簇。至少8个块，399次Rademacher符号乘以块聚合回归影响量，固定seed 20261103，得wild-score区间与零假设经验P。这保留块内共同变化，依赖块间近似独立。少于8块只有效应，无空间P/CI。未报告普通独立spot P值。R12和dataset-local内交互仅在至少100query的层内估计，并保存实际配对和空间块数。没有在稀疏层中重定义目标。局部标签按dataset限定。

每片先得效应。按已有单位或核验患者合并重复切片，相关取Fisher-z后等权，交互标准化系数取等权均值；再按dataset等权。区间通过399次嵌套重采样dataset与其内单位计算。跨dataset检验对dataset效应作零均值t检验，至少3个dataset；小样本和异质性需结合方向数解读。BH分别按来源、人口层级、MP10距离/程序距离/MP10程序/交互/界面族及主P95亚群与敏感性参照分层。完整结果均保留。within-niche图仅为描述性均值，正式片内空间推断以明细的块数和P为准。

## Reactome与单细胞
每片用完整RNA LogNormalize表达，残差化MP10对上皮组成、深度和位置；用基因表达与残差MP10的标准化协方差排序。片内秩百分位先按已有单位平均，再dataset等权。61片有足够query建立排序，肿瘤来源进入Reactome细化。MSigDB 2026.1.Hs完整Reactome，匹配15–500基因且覆盖≥50%的所有通路均测试，fgseaSimple 10,000基因集置换、seed固定。完整阴性和排除范围保留。它是基因排序层面的富集，P/FDR没有患者重复验证含义。部分dataset有较多并列秩，富集需结合此限制；本轮没有针对结果换排名或挑通路。

单细胞读取本项目正式reference_v2_raw_counts_minimally_harmonized.h5ad。沿用reference_celltype_v2的七个亚群及Epithelial，orig.ident为58个已有样本标识。sample_id字段大量缺失，不用缺失值推断映射。每样本每类型固定随机最多100细胞，seed 20261103，合计19,425细胞。完整基因轴上的55套官方UCell，参数同空间；结果按样本–类型中位数汇总。这是分层抽样来源支持；没有构造空间–单细胞患者配对，也没有新算单细胞MP10来替代冻结空间评分。CancerSEA-X是可选项，本轮未纳入。

## 输出
tables/effects.tsv.gz 为逐片完整效应，summary.tsv和unit_effects.tsv.gz提供身份分层汇总；targets、interface、niche、strata、niche_interactions、gene_overlap及Reactome/单细胞表供复核。各切片 *_program_band 为组织与Epi-query分开的55程序分布，*_curves 为组成及Epi程序，所有niche中MP10均来自有效Epi queries，至少5点。全部原始值及新几何在spots/。图册按dataset收录每片目标和界面；overview.pdf包含全部55程序热图和分箱曲线。绘图遵循jins-routine-plotting的白底Arial/矢量规范；空间和效应矩阵无匹配模板，使用专用Matplotlib实现。程序用8进程、秩计算6进程、单细胞1进程、绘图6进程，BLAS/OMP为1，避免全量对象复制到各进程。
