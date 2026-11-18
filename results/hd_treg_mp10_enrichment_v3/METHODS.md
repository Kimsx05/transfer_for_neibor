# 方法与解释边界

本轮是v1/v2之后的探索性区域空间分析，统计单位为T32和Epi64聚合区域。五片逐片处理，患者关系保留unknown，MPBC1/2原组织掩膜仍未知。输入与代码路径见input_manifest.tsv和definitions_audit.tsv。只读取原分数、区域与成员，不新增全表达副本。

## 冻结定义

T32来自完整supplied表达网格的32μm非重叠聚合，不先挑T类bins。T_eligible=nFeature≥325，CD3D/CD3E/CD3G/TRAC/CD2至少2个检出且UMI≥3。连续Treg程序使用全部T_eligible，包含计算零值；FOXP3-supported仅额外要求raw_FOXP3>0。严格标签沿用v1：T_eligible、11基因检出≥3、core≥1、score>0.19955389201641047。阈值及labels不被本轮校准替换。

11基因：FOXP3, IL2RA, CTLA4, IKZF2, RTKN2, CCR8, TNFRSF18, TNFRSF4, ICOS, LAYN, MAGEH1。完整背景UCell average ties/maxRank1300；分数直接复用。FOXP3原始非零但程序为零可能由秩窗口造成，示例按原v2秩审计保存。

MP10 Top20：HMGCS1, MSMO1, INSIG1, FDFT1, SQLE, FDPS, IDI1, ACSS2, ELOVL5, LDLR, MVD, SCD, HMGCR, ACAT2, SLC20A1, ACSL3, ACSL1, ADIPOR2, MXD1, NPC1。复用Epi64资格、18085完整基因背景、平均并列秩K905与历史归一化。H05/H10为合格Epi64每片MP10严格大于linear分位数q95/q90，不拆并列。

## 校准分母修复

复用v2的expanded_mixtures.parquet及原partition。每个background×competitor×depth的训练阴性条件，分母为所有T_eligible阴性，不要求先过multigene/core。先应用完整候选规则再计FP。

允许FP数=floor(0.05×分母)。将预门槛合格负例score排序，找到使严格`score>cutoff`的FP不超过允许数的最小非负cutoff，所有非空训练条件取最大值。检验前一个可达较低cutoff会违反至少一条件。零分母率NA；训练外的旧sample_holdout只评价，不参与选择。Wilson区间是模拟样本条件下的描述区间，不代表患者独立区间。结果0.27019811718420683仅作为校准候选，不更新任何HD标签。21600缓存逐条件输出及分子分母均保存。

## 物理邻域与来源

复用原μm坐标(向右x、向下y)，同section和原输入组件。质心距离≤50或200μm，两个累计圆。成员映射重新join核验：若T32和Epi64共享任一原bin，整对去掉；不拆分表达，不改变Epi64分数。输入组件不等于独立病理组织块。

U_r=T_eligible且原supplied-circle coverage≥0.5且≥1独立合格Epi64。没有可测Epi不纳入非high组。优先拒绝顺序coverage不足→半径内无合格Epi→所有候选共享来源；另存每条独立失败数，不能把重叠失败相加。

## FE、PR、OR与限制置换

唯一中心列联表a=阳性且near_high，b=阳性且有Epi但非near_high，c/d为其余可评价T对应两组。N=a+b+c+d。

FE=[a/(a+c)]/[(a+b)/N]；PR=[a/(a+c)]/[b/(b+d)]；OR=ad/(bc)。不加伪计数；合法0、正数除0的Inf、0/0的NA分别保留。

固定坐标、Epi分数/high、资格与来源mask；只在U_r中的T标签行置换。分层=原component×组件内nCount是否高于中位数×T-panelUMI是否高于中位数×Epi邻居总数是否高于中位数×coverage≥0.8。切点不用标签、MP10分数/high或效应。两种标签同时按相同行置换，保持交集结构；9999次，seed=123+切片序号×10+radius，分块批量128。小层不依据结果合并。

至少20个处于混合标签层的可交换对象、5个可移动阳性、2个混合层才给正式P；没有足够支持仍保留计数和描述性null。双侧P=(1+sum(abs(null-null_mean)≥abs(observed-null_mean)-1e-12))/(B+1)，离散ties纳入。Z=(O-null_mean)/null_sd，sd为样本标准差；零SD的Z/P均NA。边数O使用同一个有效配对图，其中心按边数贡献，区别于FE。

主FE家族计划10项(FOXP3×H05×5×2)，对可计算P统一BH，所有计划行保留。H10/strict其他FE共30项辅助家族；边统计40项一个辅助家族。其他分析不另增正式P。非交换部分会限制零模型解释；置换不保留完整空间自相关。

## 空间区块、热点与连续读数

400/800μm固定原点方块，component分开，按块有放回抽999次(seed123)，中心不重复展开表达。FE/PR/OR记录每个指标NA/Inf抽样比例、无阳性/无high比例；至少12块、3个阳性块、≥95%有限draw才报告有限draw的2.5–97.5百分位区间。剩余无效draw不隐去，区间是这些条件下的有限draw区间；不满足则NA。普通bootstrap尾部未作为P。

热点预先定义为同component阳性中心间≤64μm连接的连通组；逐组移除其阳性中心，仅作影响诊断。另复用v1固定candidate ROI原256μm框合并范围，剔除其中全部T中心，未按新的512μm绘图视野扩大；原邻域读数保持不变。这一额外诊断延续历史固定ROI，不以新效应选择范围。

连续主读数是独立Epi64 MP10等权均值；H05/H10邻居比例为辅助。保留Spearman及秩标准化线性效应；候选协变量固定为T-lineage、T侧UMI、Epi侧平均深度、Epi邻居数、coverage和原背景。按固定顺序对现有设计残差方差>0.1才纳入，以统一处理常数/共线性，非按P删变量。Treg自变量不因strict或FOXP3标签筛选。秩在原U_r固定后进行区块重采样；至少50中心/12块/95%有效draw/足够自变量残差变异才报告区间。全部实际设计矩阵保存；没有新连续正式P，v1原P/q单独归档。

## 距离曲线与K/L

适配曲线的总体固定为U200，每10μm从0至200累计，同来源mask和9999限制标签置换。指标为独立high Epi邻接边数/阳性中心，及观察/随机均值比；点式2.5–97.5 null包络不是同时包络，不在峰处追加主终点。

标准形式K/L仅在前三片可核验窗口做描述性估计。窗口为原supplied in-tissue 8μm格的并集，非bbox、非完整组织面积；验证阵列到μm等距近似。为形成不共享表达的两点型集合，全局剔除与任意合格T32共享来源的Epi64，记录对象损失。8μm像素distance transform给border距离(中心距减半格)，以距观察窗边界≥r的A作为内点。每组件K_AB(r)=观察面积×内点A至B的配对数/(内点A数×全部B数)，L=sqrt(K/pi)。至少10A、2B、5内点A，否则NA。

这是透明自实现的homogeneous border描述估计，不是调用未经修改的spatstat Kcross检验；观察强度非均匀未被该公式解决，L-r>0不代表特异吸引。随机标记适配曲线是另一个条件比较，不是对标准K的inhom校正。MPBC没有原组织窗，不执行标准K/L。参考[Kcross文档](https://search.r-project.org/CRAN/refmans/spatstat.explore/html/Kcross.html)。

## 双向mixing

观察总体预先取coverage合格T(包括无Epi边者)和其组件中全部合格Epi(包括无边high)。这与FE的U_r不同，不能互换分母。A为标签T、B为highEpi以及反方向各算。I_AB每对一次且整对共享来源排除；I_AA同类无序点对，无自环。MS=I_AB/I_AA；NMS=[I_AB/(n_A n_B)]/[I_AA/(n_A(n_A-1)/2)]。n_A/n_B包含无边点。

无A/B或无AA分母→NA；有分母且AB=0才是0。999次同分层、同pairmask的标签随机基线仅描述，不给额外正式P。记录退化null比例。此处为区域混合，非细胞接触；参考[SPIAT](https://trigosteam.github.io/SPIAT/reference/mixing_score_summary.html)，项目对来源与观察机会的限制是额外适配。

## 实现与图件

SciPy KD-tree只用于必要点对，不构造N×N稠密距离矩阵；既有稀疏配对复用。总并发≤4，数值BLAS/OpenMP设1，Python affinity至4个CPU；有两项单线程绘图短时并行。此次是既有表的有限计算，复用用户先前明确直接计算授权，无重上游作业；rerun默认仍遵守夜间计算窗口，已有授权需显式HD_ALLOW_DAYTIME=1。

图按jins-routine-plotting公共规范，分布直接使用其Violin_distribution模板；空间、森林、距离图无适配模板，使用matplotlib独立图。白底Arial，PNG300dpi，PDF文字/图例矢量、密集点栅格。有图例输出完整/主体/图例。直接重绘入口不调用评分/关联。

项目FE不是[Squidpy co-occurrence](https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_co_occurrence.html)原生统计；限制标签置换也不冒称[Squidpy默认neighborhood enrichment](https://squidpy.readthedocs.io/en/latest/api/squidpy.gr.nhood_enrichment.html)。这些官方资料用于明确统计单位与零模型差异。
