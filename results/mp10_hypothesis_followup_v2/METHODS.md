# METHODS — actual execution

## Design and units

冻结时间、版本和技术修订见 `config/analysis_plan.yaml`、`M0_FROZEN.json`、`technical_revision_*.json`。这是已看过发现性结果后的定向内部检验，不是前瞻性独立验证。H05主，H10/H50敏感性；全部沿用原 exact c2l reference membership、未舍入分组和 L50，不重新划组，不修改评分。主profile≥50，敏感性≥30；≥4记录患者单位、满秩、残差df≥3、有效暴露变化。主T4463209a、替换b敏感性；其余映射保留项目 provisional 状态，不能将记录单位当作已确认独立患者。

原始计数配对DE继续引用 edgeR TMM robust QL `~ specimen + state` 的全基因结果；不为本轮候选重做两基因FDR。新增样本程序模型对 High−Low 差值用 OLS HC3，双侧t检验和同一HC3残差df的95% CI。跨细胞模型为 `receiver ~ dataset + standardized H05 sender`；H05−L50 sender差值及其他设置为登记敏感性。连续模型 h≥1−1e−10 时 HC3 NE；cycling为quasi-binomial，使用各模型离散度和残差df t区间。NE保留NA，真实零不替代缺失。

## Scores, gene sets, masking

复用 v1 `objects/scoring_matrices.qs` 的真实全谱系 TMM-logCPM 背景与reference profiles；没有截取候选并集重归一化。GSVA 2.4.9 ssGSEA，alpha=.25，normalize=FALSE。每集合在main-a Epi_all≥50或接收端原参考集固定mean/SD；所有High、LODO及患者替换均沿用该尺度，不中心化每个delta。

正式13过程+1EC运输补充均来自冻结MSigDB2026.1.Hs/Reactome资源；完整成员和GO非编码项不静默移除。程序评分≥5可用成员且≥50%预期集合；GSEA minSize15/maxSize500，signed全基因DE排序，seed123，nPermSimple10000，eps1e−10，scoreType std。EC整套Hallmark因数值诊断统一提高至100000，保留10000结果和警告；49可检验+1小集合NE全部登记，BH计划50。极小P为计算下限，不解释精确尾概率。

Full为正式过程描述主结果；M20仅从下游集合及排序/评分宇宙移除Top20，不改变RNA原counts、library、TMM或分组。接收端Full主，M20仅明确敏感性。旧自编面板来源单列，不冒称正式集合。基因标识ACAT2未替换成SOAT2，NECTIN2/PVRL2不重复计数。

## Direct sterol comparison

同一H05/L50配对中，D=(zCHOL_H05−zCHOL_L50)−(zFA_H05−zFA_L50)。主CHOL_SYNTH对FA_ACYL_COA；Bloch补充。每个mask两项一个家族。相同评分算法和预设参考；两个集合各自参考SD而非通量单位。200次固定seed123匹配：可用基因数相同，在全Epi参考背景标准化的平均logCPM和检出比例上最近邻无放回匹配，随机匹配顺序，所有结果保留。它检验集合构成敏感性，不是200独立重复。

## SCENIC stage audit and correction

原pySCENIC0.12.1 seed123 adjacency及精确训练cache复用；8577模块重建而不再跑GRNBoost2。rho=.03、非dropout-mask、仅activating、min_genes20；importance .75/.90、top50targets、top5/10/50TF-per-target；ranking5000、AUC .05、NES≥3，原hg38/refseq-r80/mc9nr 10kb ranking及匹配v9 annotation，motifsimilarityFDR .001。逐模块ID含TF/context/targets哈希，原始分母保留。

实际软件 `missing_fraction >= .20` 排除，即正好80%也排除。所有motif先在同输入下补算诊断再逐门核验。唯一无歧义HGNC/NCBI历史别名与ranking链接可修正；保留未映射、歧义和原名，不删基因降低分母。安装HGNC20170418与NCBI2025-09-24参考来源和hash均登记；PRAG1→SGK223获得同一NCBI157285支持。没有为SREBF2挑新ranking或放宽阈值。

修正后regulon逆映射回表达原名，保留原TF-self成员语义。AUCell在原20216训练基因顺序、完整116448 Epi normalizedRNA上按2000cell chunks计算，seed123/aucThreshold .05/1worker。Full与M20；原SREBF1前3000cells逐数值校验maxabs0。活性差异是样本配对证据，不是ChIP、成熟蛋白或独立数据复现。

## mTOR sources and QC

Epi读取v1修复后的M0031完整Hallmark家族；EC M0073为EC自身rawRNA计数对Epi Q50的dataset调整关联，保留历史primary标签；空间B3为混合spotRNA UCell，Epi位置仅决定query，不能将RNA归属上皮。空间条件区间和dataset-level q独立列出。

EC过滤交叉诊断固定两组样本与两种filterByExpr设计，逐基因输出原检出/CPM/保留/leading-edge/Top20。对原3905 EC用原RNA counts核对内皮标记、EPCAM/keratin共表达、原QC字段及单细胞计数贡献。该标记不是doublet判定；没有执行ambient RNA校正。单个代谢基因表达不被定义为污染。

## FDR, reproducibility and scope

完整登记1328行（1278主登记+50EC Hallmark诊断），所有family计划数见 `FDR_family_summary.tsv`。每High/mask/setting程序家族13；直接D家族2；跨细胞primary2、EC探索5、Treg次要6；恢复SCENIC每setting/mask/state2。敏感性另成family，NE保持NA而BH使用计划n；不选最小P。旧RNA/TF/GSEA的原q及来源另存。

`qc/validation.tsv`独立Python重算HC3、quasi、BH和约束；`tables/counts_membership_validation.tsv`校验原整数计数及精确成员；`qc/render_audit.tsv`核对18图62文件，PDF矢量Arial、约12pt、300dpiPNG，带图例分离PDF。未运行新的HD、CNV、CellChat、NicheNet或全GRN，未改变原始或历史目录。

实际项目：`/home/data/t070721/codex_workspace/Bladder Metabolism`。版本见 `logs/sessionInfo_registry.txt`、`sessionInfo_scores.txt`、`environment_versions.tsv`；代码、参数、成员与参考哈希见checkpoint及各缓存metadata。原始大对象只读；工作峰值RSS来自恢复AUCell约44.61GiB，内部线程1，总并发计算不超过4。
