# v3方法与冻结规则

## 输入与范围
原NMF output_full_cc_seed5_v2的123699个上皮分析细胞，包含Normal及原已纳入的Cycling，不是TME atlas。沿用57个dataset|orig.ident biological samples、5 datasets、22436 Normal。原Seurat、NMF、MP top50、usage、v1/v2标签/降维只读；没有重跑NMF、UCell、PCA、Harmony、UMAP或BAL。患者映射仍derived/provisional，不声称独立患者重复。精确文件/版本、大小/mtime、cell ID、软件版本、图/分数校验见00_audit。

## 分区与high
M_RAW与M_Z按cell ID精确核验，14列固定。Z与FULL全局mean/sample-SD显式列运算一致。两套原SNN复用，原矩阵SHA256核对；k20、annoy/50trees、Euclidean、l2.normFALSE、prune1/15。Seurat algorithm4、leidenbase、RBConfigurationVertexPartition、modularity.fxn1、n.iter10、group.singletonsTRUE、seed42；固定12个resolution×2表示，每次前设seed。原RAW r.4和Z r.2/.4/.6都按成员集合复现，允许群号置换。只有resolution改变。

在新扫描前冻结high逐细胞逻辑及精确二进制阈值：NMF使用原usage_raw>=保存的16.820936842046201，共12370；UCell使用FULL原MP10 quantile(type7,p.90)=0.25970159376059676，>=并列纳入，共12372，交集7002。这是v3新定义的UCell-high，非v1/v2旧规则。旧C8 TP2014/3937核验通过。q80/q95只对预先确定的全局F1冠军描述，不进入主排行榜。

## 覆盖、候选与对照
每个非空单群计TP/FP/FN/TN，P=TP/n_cluster，R=TP/n_high，F1=2TP/(n_cluster+n_high)，enrichment=P/(n_high/n_FULL)。每分区分别选两种high F1冠军，相同者合并；排序F1、min(P,R)、P降序，最后canonical数字群号升序。全部群和两种high的分区内/全局P-R Pareto前沿保留。全局推荐按同样指标，最终representation/resolution确定性排序，在DE之前冻结；不按marker替换候选，不合并群冲榜。

局部竞争群从对应原SNN计算sum(SNN[C,K])/(n_C*n_K)最大者，排除自身、canonical群号破同分。所有候选均使用同一图规则，没有逐候选切换质心。另保存NMF/UCell最大群外TP承载群，并按冠军所属high（双冠军优先NMF）最多加一个不同于最近邻的carrier对照。对照资格不足不更换。

sample/dataset/source覆盖使用完整来源池，无30细胞过滤。high>0且候选=0的R为0；high=0的R为NA；候选=0的P为NA。pooled和sample等权分开，NA/零召回分母全保存。来源占候选比例与候选占来源比例分别输出。

## RNA marker与资源
原RNA90011×123699稀疏counts重新核验非负整数；已有data与普通LogNormalize(count/library×10000)核对，展示与计数模型分开。原Seurat只加载一次，按FULL sample总量和每partition所需sample×cluster顺序稀疏聚合；rest由同sample总量减candidate，核验counts非负及守恒。没有dense全基因×cell矩阵。线程/worker≤2；重型图、聚合、模型顺序执行，BLAS1/fgsea1；估计RNA峰值约30–40GiB、输出缓存约数GiB，实际资源见00_audit及报告。

模型与v2一致：每侧≥30细胞的完整sample对，至少3对，design~sample+group，Right reference、groupLeft正向；满秩/剩余自由度/非零库核验。filterByExpr(group=group,min.count10,min.total.count15,large.n10,min.prop.7)，TMM，estimateDisp robustTRUE，glmQLFit robustTRUE legacyFALSE，glmQLFTest。每比较所有通过过滤基因单独BH。NOT_TESTABLE保留原因、资格和描述；不是无差异。完整基因/设计/效应/检测率/方向/输入成员关系均保存。只按两侧完整成员hash及完全相同模型配置复用；同细胞集意味着同sample资格。Jaccard高不复用DE。

短名单：FDR<.05、log2FC≥.5、目标样本等权检测率≥.1、表达正方向≥.7；按dataset方向、sample方向、检测率差、logFC、ID排序，至多20。展示最多3，优先有已知gene_symbol、同时通过rest与固定最近竞争群门槛的非输入基因，兼顾已冻结功能类别和不同基因名前缀；不足不补。局部NOT_TESTABLE时若展示rest基因必须标为仅rest支持。额外carrier的完整效应保留，不能因此改冠军。基因前缀只是展示去冗余辅助，非正式基因家族注释。最终关键候选与旧C8还在真正共同合格sample上汇报描述性表达方向。

## 统一功能分析
复用冻结目录11A_msigdb_gene_sets_cache.qs，msigdbr26.1.0缓存签名7d3852357e591330；提取9377个人类GO:BP/Reactome集合。缓存未独立记录数据库发行版本，故按精确文件SHA256锁定，不宣称额外版本。排名sign(log2FC)*sqrt(F)，保存每个受检gene_id；映射到已知SYMBOL/唯一ENSEMBL符号，非有限/未映射记录排除，重复符号统一取logCPM最高、gene_id破同分，不选最大效应。fgseaMultilevel sampleSize101、nPermSimple1000、minSize10、maxSize500、eps1e-4、nproc1、seed42、scoreType std。P值低于1e-4不继续追求极小值，NA保留；每比较所有两库通路一起BH，非全候选搜索统一FDR。全表NES/FDR/leading edge及输入基因贡献保存。未做ORA。

展示固定区分fatty_acid、cholesterol_sterol、other_lipid、response_stress、proliferation、epithelial六类通路名关键词；正/负NES分别展示、每类最多3，类别可重叠，不能相加计独立通路数。完整9377集合的实际受检富集表一直保留。所有marker、富集和非输入基因都是同批发现RNA证据，不是外部验证；MP10参与建群/high选择带来定义性描述及选择偏倚。

## 图与参考
原C8重复模型测试基因集与v2完全相同（15902基因），log2FC最大差1.42e-9、FDR最大差9.06e-9，表达与检测率一致；微差来自行排序下的浮点计算。

同一分支共享原UMAP连续底图，RNA用已冻结umap_allgene；dataset取全局坐标子集，不重拟合。全部5dataset保留。非负连续plasma，精确0黑、NA灰，FULL范围固定；usage仅显示log1p，统计始终raw。RNA DotPlot是sample等权mean LogNormalize/detection，非usage；无网格。每图透明PDF和300dpi PNG，容器背景也核验；图数据/来源索引保存。成员流向的纵坐标为群大小堆积布局，不是表达或UMAP距离。

- [fgsea官方教程](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [edgeR官方入口](https://bioconductor.org/packages/edgeR/)，执行API以本地保存help为准。
- [Seurat FindClusters](https://satijalab.org/seurat/reference/findclusters)
