# 方法与执行记录

## 输入、分母与研究单位

本轮纳入 Cervilla、D1_Johnson、UC_SCC、MPBC1、MPBC2。定向核查 data、已有 asset manifest 和加载入口，没有发现可另行纳入的独立 HD 切片；不同分辨率文件不是新增切片。普通 Visium63 未运行。输入表为 `tables/input_manifest.tsv`，每片原始路径与文件身份保存在 `prepared.json`。5 张切片不等于5个已确认独立患者，患者映射仍未知；MPBC为两个额外MIBC微乳头/尿路上皮混合标本，不能按论文标题归为NMIBC。

复用完整 supplied 8 µm raw integer counts，共18,085基因。Cervilla/D1/UC有原组织范围；MPBC1/2来源为上轮按固定原点从实际2 µm counts聚合的8 µm，原组织掩膜保持未知。没有将 filtered 文件名或表达阈值当作组织掩膜。聚合前未删除低UMI bins。确认组织外才排除；未知组织归属保留。所有结果是区域测量，非真实Treg细胞数或密度。

坐标、物理尺度和Cervilla/D1图像变换复用已验证的节点表与 v2 transform。UC使用已有hires图及scalefactors；MPBC仅相对阵列坐标。x向右、y向下；未提供图像/表达位置保持缺测。支持面积是实际成员bins×64 µm²，MPBC不能称为组织面积。输入连通组件不是病理分区。

## 区域构建与来源分离

T主测量在原阵列0点的32 µm非重叠格内按8连通、同输入组件汇总；16 µm为唯一尺度敏感性。保存全部区域，不只旧T/NK或有标签区域。每个独立分支核对逐基因和全局 `included+excluded=input`。未知、实测零、测量不合格独立记录。

上皮支持使用未经平滑的10基因面板(EPCAM,KRT7,KRT19,UPK2,KRT5,KRT14,KRT6A,KRT6B,KRT17,COL17A1)：至少2个检出且任一冻结参考预测Epi，或至少3个检出。无Treg/MP10参与此选择。仅将支持bins在固定64 µm格内连通汇总，不跨缺口。名称为“上皮富集区域MP10”，不是已分离纯上皮RNA。所有成员映射以唯一bin ID保留。

位置为成员bin等权质心。50/200 µm为精确累计欧氏半径；不跨输入组件，一个Epi区域只贡献一次。主终点排除任何共享原始8 µm bins的T/Epi区域对；包含共享来源单独作为共同检测诊断。不能把评分复制到成员后当作独立观测。邻域需至少1个独立可评分Epi区域，并且圆内实际合格输入网格支持面积/圆面积≥0.5(截顶1)；这是可观察输入支持，不是已证实的组织覆盖率。50 µm缺少目标时不扩大半径补足。

## Treg signature与参考校准

使用最终30类型精细参考 `subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad` 的raw counts；复用旧冻结参考成员20,446，18,289训练/2,157留出，91个留出Treg。留出按orig.ident样本，不称患者独立验证。类型背景保留Tfh/CXCL13-like、普通/活化CD4、CD8、NK、上皮及其他非T谱系。

只在训练部分选基因：候选中Treg检出率≥0.08，平均每万UMI表达>1.2倍最高竞争T/NK类，且>1.5倍合并非T背景。最终11基因：FOXP3,IL2RA,CTLA4,IKZF2,RTKN2,CCR8,TNFRSF18,TNFRSF4,ICOS,LAYN,MAGEH1。LRRC32未通过；排除MP10及预定通讯读数TIGIT,NECTIN2,CD96,NECTIN1,PVR,CD226。未按HD关联选基因。

文献复核为De Simone et al., Immunity2016 Table S4的完整309基因，来自实际补充Excel；保存原始表、来源及逐基因修改。剔除预定读数/MP10并要求参考和HD共同可测后，评分版本283基因。原文为CRC/NSCLC肿瘤Treg，不是BLCA专属分类器；含活化和增殖相关基因，单独解释。不是自行拼接marker冒充原文signature。

参考验证含样本留出、预期500/2000/5000 UMI thinning、leave-one-gene-out、无Treg的70%上皮+30%竞争T RNA阴性混合以及10%/30% Treg RNA混合。每组100次。混合比例指RNA贡献而非细胞比例。训练期特征冻结后才用于留出；不将验证集重新选择特征。

## UCell与测量资格

实际R4.5.3、UCell2.14.0；软件清单及sessionInfo附后。按五片T支持32 µm区域的最低中位nFeature向下取整100，限制500–3000，得到共同maxRank=1300。保留平均并列秩，包括零值；使用完整基因背景。没有逐片调到最有利参数。

原生公开ScoreSignatures_UCell在80个完整表达示例核验。全量使用UCell2.14 `u_stat_signature_list`原生评分core、明确maxRank=1300，对在完整18,085背景计算的精确秩分块执行。原因是公共precalc.ranks会从缓存最大值推断cap，精选基因/半整数并列秩缓存不适于该推断；不是换用自定义评分器。数值比较修复了浮点32位约3e-8差异，double实现与原生core/历史MP10达到<1e-12。正式区域字段使用原生UCell结果。未保存全基因标准化或稠密排名矩阵。

T可测量：nFeature≥325、≥2个T-lineage基因、T-lineage UMI≥3。这不是单细胞QC删除规则，而是本轮区域评分资格；所有原始区域保留。Treg相关连续分数为主，合格T区中的计算零值也进入主分析。描述性Treg-like层级另需主signature≥3基因、核心FOXP3/IL2RA/CTLA4/IKZF2中≥1、分数超过训练竞争T的Q95=0.199553892。这个层级不等于可信单细胞身份，不是每片强取最高百分比。主分析不要求达到这个Treg-like层级。

## MP10与冻结顺序

Top20严格复用项目冻结文件：HMGCS1,MSMO1,INSIG1,FDFT1,SQLE,FDPS,IDI1,ACSS2,ELOVL5,LDLR,MVD,SCD,HMGCR,ACAT2,SLC20A1,ACSL3,ACSL1,ADIPOR2,MXD1,NPC1。历史公式：全18,085基因平均并列秩，sum(max(905-rank,0))/sum(905-1,...,905-20)。无默认AUCell替代、无平滑或HVG裁剪。每个评分区域≥905检出基因、≥5个Top20成员检出；原始Top20未更换。仅为区域测量资格，不把旧“5细胞”改成“5 bins”。合格计算零保留，不合格不是Low。

单次半深度测试和分数-检测量关系用于测量限制记录。全片合格上皮区域type7分位数生成L50≤Q50、H50>Q50、H10>Q90、H05>Q95，不拆并列。连续值仍为主。

配置于2026-11-17 00:28 CST冻结，早于正式空间关联；参考与几何先行。ROI使用T证据、空间分散及T强度匹配对照，无MP10关系参与。24个最终ROI清单在关联计算已启动但结果尚未查看时保存，明确记录在models/roi_freeze.json；不是在关联图上追选重叠热点。少量候选(Cervilla/MPBC各2)全部展示；D1/UC按预定分数范围及空间分布选4并配2对照。

## 逐片统计

因变量为独立Epi评分区域MP10等权均值，自变量为连续主Treg UCell。效应是标准化平均秩OLS斜率，不是表达倍数或细胞比例。未调整和固定协变量调整均保留。固定顺序：总T score、T UMI、Epi平均UMI、Epi支持面积、输入圆覆盖、原始Epi-Fib表达背景。常量或VIF≥10的后续变量删除，具体实际协变量写每行，不按方向或显著性选择。Tfh/Activated_T附加条件作为单独敏感性。

空间400 µm区块根据最大邻域及边际相关诊断事先固定，800 µm为唯一更大区块敏感性。至少50个T区域、12个区块及足够不同分数才做推断。CI为499次空间区块pairs bootstrap；P为限制零假设模型的HC2残差、区块Rademacher wild bootstrap(499次，seed123)。不是普通bootstrap系数尾部比例，也不全片逐bin打乱。区块模型仍是片内探索性推断，边界重叠和未测混杂不能完全消除。

主家族为5片×2半径，10行全部保留、7行可检验，BH对非NA P执行。文献、16 µm、包含共享来源、总T、Tfh、活化T、反向分析、更大区块各自单列家族；全结果保留。反向分析以Epi区域为单位，不当独立重复。无法确认患者独立性且测量支持不同，不输出合并总体P或meta效应。

## 旧细胞、SpottedPy与通讯前检查

Cervilla复用上轮1,209个细胞限定Epi区域MP10，不重跑分割或细胞评分。通过原2 µm bin-to-cell映射及原阵列floor/4恢复到8 µm来源，保守排除任何共享8 µm格。除原覆盖分析外，使用完全相同247个T中心和相同协变量比较两种MP10来源；50 µm共同仅7个，描述性。其他片没有可直接复用的同类完整细胞限定MP10区域结果，不补造。

SpottedPy0.1.7源代码锁定SHA，原样提取calculateDistancesHelper执行，不安装/替换其旧依赖；这是原生距离helper适配调用，不声称运行Gi*。预定Treg-like位置与Epi H50/H10/L50位置作物理µm最近距离，每次一个源位置、同输入组件、排除共享原始bins，未找到或距离超出可观察边缘记NA/截尾，不填最大值。仅辅助可视化，不代替主统计。

冻结scRNA候选为TIGIT_NECTIN2、CD96_NECTIN1，来源mp10_cellchat_high_definition_sensitivity_v1 H05表，CellChatDB.human1.6.1。保留数据库编码的ligand/receptor字段，但TIGIT/CD96受体型分子与上皮结合伙伴的生物学方向不能由该字段顺序证明。只核对raw表达、来源分离、50/200 µm邻近和MP10组描述，不算配体受体得分。没有运行CellChat/LIANA/NicheNet/TF等。

## 绘图、资源、验证

采用已安装jins-routine-plotting规范。其常规模板不足以覆盖本次空间图册/森林图布局，使用matplotlib自定义布局，Arial约12 pt、白底300dpi PNG、PDF文字和图例矢量、密集点及H&E栅格。有图例输出full/panel_only/legend_only PDF。黑比例尺64 µm，红圆50 µm；ROI256 µm。橙色空心方块标中心T区域实际成员bins，青色空心方块标可评分Epi区域成员bins；不画跨缺口的凸包或插值边界。空白是无输入，灰底为supplied，不可评分在覆盖图单列；分数零不等于缺测。全片底图可抽稀显示，所有分析和评分使用全量合格区域。

默认总CPU亲和性4，BLAS/OpenMP/Numba1、R --vanilla及运行时线程控制，稀疏计数/空间索引。用户明确“没有必要等，直接算就可以了”后本轮允许日间继续，记录runtime_authorization.json；没有升级共享环境或重复上游BANKSY/SingleR。最大记录单进程RSS约5.1 GiB，原生UCell全量14分09秒。资源表列实际计时步骤，不把多进程峰值相加冒充实测总峰值。

实际重读所有区域/成员/邻域、核对唯一ID、坐标质心、counts来源、逐基因守恒、评分来源、缺测、共享bin排除、50/200几何与同组件、统计区域不复制；Python/R loader实际读回五片。图像检查见ROI_VISUAL_REVIEW.md。上游原counts/QC不修改。

## 方法来源

- [UCell参数](https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html)：以本机2.14实现核对，不套用默认maxRank。
- [De Simone2016](https://pubmed.ncbi.nlm.nih.gov/27851914/)；[完整XML与附表入口](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5119953/fullTextXML)。本地sources/DeSimone_tables.xlsx Table S4为309基因来源。
- [Oh2020 BLCA T细胞](https://pmc.ncbi.nlm.nih.gov/articles/PMC7321885/)：背景参考，不将自组panel称其完整signature。
- [SpottedPy作者源码](https://github.com/secrierlab/SpottedPy)：本地源码快照与实际helper执行记录。
