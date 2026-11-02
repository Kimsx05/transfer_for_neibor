# MP10 空间距离先导 (v1)

## 取片与输入

每个 dataset 直接采用 v8 已登记的代表切片。v8 在有至少 30 个 E0 中心的切片中选择中心数最接近 dataset 中位数的一片，同距按 section 排序。7 片均有 H&E、RNA、30 因子后验均值、v3.4 位置和两个候选侧。记录为肿瘤研究标本；局部候选位置及匿名样本不能作为独立确诊患者身份。选片理由、身份及实际覆盖分别见 [名单](pilot_sections.tsv)、[输入文件](tables/input_manifest.tsv) 和 [覆盖表](tables/coverage.tsv)。

以 section_id 和完整 barcode 一对一连接逐点表。保留原 position_id、spot_uid、图像像素 x/y、已有名义微米 x_um/y_um 和 tissue_component_id。H&E 面板直接使用 v3.4 匹配的低分辨率图及 geometry.json.gz 的像素边段。沿用原 v3.4 lesion_location、region_role、精确 signed_distance_um (内部负、外部正) 和 distance_band_200um，原有 NA 保留。

## 作者代码与项目适配

Oral 和 cSCC 脚本按定义的成纤维目标点取 min(Euclidean)，cSCC 对上皮支持的点用 FeatureScatter 作连续分数–距离图。原目标来自样本内的簇或区室筛选，cSCC 示例另用 Epithelial > 0.2。Lung 脚本将区室、组成与 AddModuleScore 程序图并置。三个原脚本原样保存在 authors/。本项目保留中心距离、连续图和并列展示，改用冻结的膀胱项目分数、目标、坐标单位及组织碎片限制。没有执行原作者整套工作流。[论文](https://www.nature.com/articles/s42255-026-01514-y) · [原脚本](https://zenodo.org/records/18790615)

**主分数**沿用区域工作流的 mixed-spot `EpiMP10_Top20_AUCell`。主查询取该区域接口标为 `analysis_eligible_legacy` 的 Epi anchors 且原始 MP10 可用，直接读名单。R12 的 common_anchor、v8 更严格的 E0 (`score_eligible`) 与旧 RCTD anchor 另列，没有用这些环覆盖条件删去合格的原始分数查询。common_anchor 和 E0 子集只在 [描述表](tables/anchor_subset.tsv) 中单独比较主效应。另存 spot_MP10_model_universe、epi_MP10_model_universe、epi_score_valid、原评分定义、覆盖及状态；Epi-component 图只显示其原有效点。原 20 基因覆盖按片见 [原分数覆盖](tables/raw_MP10_coverage.tsv)。

**目标**采用 v2_30type 的 abundance_mean.csv.gz。每个 spot 以所有 30 项 posterior mean 的和为分母。CAF 汇总 Matrix_CAF、iCAF、rCAF、vCAF；TAM 汇总 C1QC TAM 与 SPP1 TAM；Treg 使用单一 Treg 因子。FOLR2 没有独立空间因子。Treg/T-NK 按项目 major_map 的全部 T_NK 因子作分母。分区组成是各有效 spot 比例的等权平均。估计丰度保留原模型单位。target 中位估计丰度 <1 的探索性标记仅便于阅读，没有去除点或改变阈值。

每片用正式组织表中 center_qc 有效、未标 quality_excluded 且 30 个估计有限非负、总和为正的 spots，计算组成线性分位数。所有 >=q80 且 >0 的点入选，包括并列。全零、缺失或无变异不定义 target。主分析为 P80 全组织且含零；另分别做 P80 全组织仅正距离、P90 全组织含零、P80 仅 margin target 含零。margin 限制直接用 inner_margin / outer_margin，不在边界带内重新估阈值；查询位置不限边界带。

每片每个组织连通片用 cKDTree 精确欧氏搜索，在各有目标的连通片抽取至多两个查询与逐点 min(Euclidean) 核对。重叠点为 0，没有同片同块目标时为 NA，并保存 reason 和 query_is_target。病灶边界距离缺失不影响最近目标计算。按片报告 Pearson r 与 Spearman ρ，正值表示离目标越远 MP10 越高。散点线为片内线性描述，箱线图以 0、(0,200]、(200,500]、(500,1000]、>1000 μm 分组；距离仅为分箱消除浮点端点误差而四舍五入到 6 位小数。未计算 spot 独立性假设下的 P 值或跨片合并相关。

## 分区、niche 和程序

R12 为 SAMPLE_R12_LEIDEN 的 0–8；旧图的 R12_niche1–9 等于原标签 +1。分母为该区域有有效 R12 标签的合格 Epi anchors；另报 labelled/eligible 覆盖。全组织 dataset-local 标签用其有效组织 spot 分母，candidate meta-niche 独立保留，K12 没有并入。未覆盖记 NA；单元 MP10 中位数至少需 5 点。分层距离 r 至少需 20 个可用 query。两者为描述性显示下限。环组成表直接保留冻结的 [R1/R2 原 profile](tables/source_R12_profiles.tsv)，是原 R12 参考汇总，源于冻结的 q05 环组成；仅用于解释标签背景，与本批中心 spot 的 posterior mean 分开显示。

msigdbr 从本机 human MSigDB 2026.1.Hs 读取 Hallmark 50 项、C2 CP:REACTOME 的 3 项及 C7 IMMUNESIGDB 的 2 项。保留官方名称与完整成员，按 symbol 去重、与逐片全 RNA universe 交集。RNA 来自原始计数的逐片 Matrix Market 导出，沿用完整 gene axis。msigdbr 包版 26.1.1；UCell 2.14.0 使用 maxRank 1500、平均并列秩及 50 spot 分块。前两张固定切片另在同一计数的 Seurat LogNormalize 上运行 AddModuleScore (默认 24 bins、100 个对照、seed 1)，仅限 Reactome 胆固醇生物合成、Hallmark EMT 和缺氧。三组分数分别存储，以全片形态秩相关及各距离梯度方向比较。C7 DN 表示 activated 相对高，UP 表示 resting 相对高。所有程序都是 spot 转录程序，原值没有除以细胞丰度。[msigdbr 文档](https://igordot.github.io/msigdbr/articles/msigdbr-intro.html)

功能图中的距离曲线显示实际距离组内 UCell 中位数 (至少 5 点)，效应表与热图显示连续距离 Pearson r。胆固醇一致性仅在 Epi queries 比较 MP10，同一组织 spot 内比较两种方法及 Treg 组成时使用全片。主空间图使用项目 v3.4、R12 和谱系配色；没有找到本机可调用的单独定稿绘图 skill，适配脚本实现了本轮白底、矢量 PDF、300 dpi PNG 及拆分图例要求。连续图通常按该层/模块跨 7 片的第 99 百分位封顶，仅影响显示，指定比例或效应图另有标尺，AddModuleScore 用两片合并的第 1–99 百分位。PDF 保留原生 H&E 位图。

## 入口与复用

逐点整合表在 spots/<section>.parquet，包含原分区、所有 30 因子的 mean 与 fraction、原评分、niche、三种 target 范围及距离、程序分数。各 target 字段名采用 `目标_P80_all`、`目标_P90_all`、`目标_P80_margin` 前缀。数值与缺失原因逐点保存。主小表为 [target](tables/targets.tsv)、[距离](tables/associations.tsv)、[分区组成](tables/region_composition.tsv)、[niche](tables/region_niche.tsv)、[程序](tables/program_distance.tsv)。图索引给出完整图、图主体、仅图例和 PNG 的实际位置；review.zip 仅带完整图、报告和核心表。

在本项目根目录依次运行 scripts/select.py、scripts/distance.py、scripts/context.py、scripts/score.R、scripts/figures.py、scripts/programs.py、scripts/overview.py (路径均以 mp10_paper_distance_pilot_v1/ 为前缀)。Python 用 python3.12 -B，R 用 Rscript；设置 OMP_NUM_THREADS、OPENBLAS_NUM_THREADS、MKL_NUM_THREADS 为 1。逐片读取和评分，原输入只读。
