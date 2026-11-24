# 方法与复现

在当前项目根运行，阶段脚本位于 code/。本任务输出位于 results/mp10_endothelial_hypothesis_validation_v1，公共 accession 数据位于 cache/mp10_public_resources_v1。优先路径 /home/jinshengxi/bladder_metabolism 不存在，工作区及上级路径未发现适用 AGENTS.md，两个命名的资源 JSON 未在工作区找到。直接读取本地 v1/v2 的方法、脚本、输入清单、模型和结果表。冻结规则见 code/design.txt。没有复制或打开全量 Seurat/H5AD。工作线程和 BLAS / OpenMP 设置为 1，Rscript --vanilla，推断脚本绑定单个 CPU；通讯先验读取另在共同 flock 下完成。

## 输入及文献定义

原始 Epi MP10 Top20、Q50_z/Q95_z 及样本选择均读自 v1。保留 ACAT2。T4463209a/b 主取 a，替换取 b。每个受者亚型 profile 按受者类型/标本唯一，静脉等原标注保持。总内皮把同标本四种互斥亚型 raw counts 求和。组成的分母为该标本全部内皮。仅在已有 ≥50 或 ≥30 门槛、≥4 个记录单位、满秩及残差 df≥3 时作推断。旧毛细 50 Hallmark 和 v2 M20 直接复用，其定量值见 legacy_full_m20.tsv。新增总内皮及可估静脉使用 filterByExpr(实际设计)、TMM、robust edgeR QL 和有符号 sqrt(max(F,0))，Hallmark 2026.1.Hs / fgseaMultilevel，15–500 个成员，按模型 50 通路 BH。

- Fan 等的 Data S1 第二表为红肌相对白肌升高且与既有 ATF4 靶基因重叠的集合；与同研究真实 Atf4 KO 分开。小鼠骨骼肌内皮/运动、人内皮实验提供生态背景。[Cell Metabolism](https://doi.org/10.1016/j.cmet.2021.07.015)
- Torrence 等 Figure 1 source data 1，第 2 与第 3 表取稳定 Ensembl ID 交集，得到原作者 77 个 mTORC1/胰岛素 ATF4 生长响应基因，第 5 表为 774 个 ATF4 依赖 tunicamycin 响应；共有 61。全集并不全是合成酶。另以 GO:0006418、GO:0003333、GO:0008652 含子项限定功能子集，属于次要功能读出。原始小鼠集合保留于 literature_original_sets.tsv，投射属于从 MEF 外推。[eLife](https://doi.org/10.7554/eLife.63326)
- FOXO1 操作能影响内皮 MYC 与代谢；本项目未测失活。[Nature](https://doi.org/10.1038/nature16498)
- RELA 敲低×TNF、可及胆固醇及 SREBP2 加工见原人 HUVEC 研究。SREBF2/KLF6 的定向细胞因子及 IFN/MHC-I 竞争解释独立核查；MHC-I 次要集合为 GO:0002474 含子项。IL6/CXCL1/CXCL8 为选定读出。[eLife](https://doi.org/10.7554/eLife.79529)；[JLR](https://doi.org/10.1016/j.jlr.2023.100411)
- PDHA1 模型复用同一 18 文库设计。选取原文 CXCL1/2/6/8 位点作补充读出。未将 RNAs 当作 PDK4 周转、PDH 通量或组蛋白乙酰化。[PNAS](https://doi.org/10.1073/pnas.2218150120)
- QPC 条件性缺失小鼠 P15 肺 CD31⁺ 内皮是功能受限参照，RNA 反应与原研究功能实验分开解释。[Nature Metabolism](https://doi.org/10.1038/s42255-018-0011-x)

## 外部数据

每个 GEO 归档分别下载处理矩阵及 family SOFT，临时文件校验后原子发布。传输截断或 wrapper HTML 见 logs/resource_attempts.tsv。恢复采用有界重试、字节范围和 gzip/归档结构验证，成功结果记录大小与 SHA256。外部样本、物种、标题、处理字段、原始列名、最终组别均在 external_sample_map.tsv。外部是文库数；HUVEC 供体块没有可靠列级映射，未按重复号配对。
 
| 归档 | 输入、映射及模型 |
|---|---|
| GSE174785 | 20 个小鼠 kallisto 文本。按 GSM 的补充文件名对应每列，核实 ENSMUST，使用 GENCODE mouse M25 唯一 transcript→gene 映射；129367/129439 个转录本可映射，单库 TPM 归属比例 99.43%–99.88%。基因 TPM 求和，log2(TPM+.5) limma-trend。分别为氧化肌/糖酵解肌 WT−KO、差异中的差异以及各基因型背景差。|
| GSE201466 | HUVEC 18 库，SOFT description 前缀唯一匹配 CPM 列，CTRL/RELA×0/4/10 h。 |
| GSE207787 | HUVEC 18 库，description 唯一匹配，CTRL/PDHA1/SREBF2×0/16 h。PDHA1 六库独立成组。 |
| GSE207919 | HUVEC 12 库，20,908 行，按 GSM 标题方括号内标签匹配打乱的 12 列，CTRL/KLF6×0/16 h。 |
| GSE121770 | WT/KO 各 4 库，以 GSM 标题匹配列；实际 ENSMUSG、整数 HTSeq counts、TMM/robust edgeR QL。 |

三个 CPM 归档使用 log2(CPM+.5) limma-trend、robust eBayes。非 counts 输入在至少最小组文库数中表达≥1 后保留。使用 HGNC exact approved 优先，其次唯一 previous/alias；多行同符号 CPM 求均值，避免重复基因座膨胀。鼠人取 babelgene 22.9 全鼠表 ≥3 来源支持并且两边 Entrez 各唯一，留存映射表；未按字母大小写推断。模型分别输出 KD−CTRL、CTRL 时间响应、(KD刺激−KD基线)−(CTRL刺激−CTRL基线)，外部全基因表保留原方向及每个比较的基因 BH。

每个外部比较在映射后按有符号统计量排，正负各取最多 100 个成员。并未把所有成员声称为 DEG。KD 及其交互在投射时反号以表示依赖；ATF4 为 WT−KO，QPC 为 KO−WT，TNF 为 CTRL 晚−早。所有比较统一在五实验与项目各自可检验的交集背景上用 fgsea（最少 10 个成员）。同一受者/遮罩的 24 对比×2 端单一家族 BH，缺失按计划数保留。全基因 Spearman 和符号符合率为描述性。固定排名、方向、UP/DN 及外部显著性均有表。另对 QPC 自身测试原 Hallmark 家族。文献三程序集在项目完整排名上富集，单独 BH；ISR 尺寸上限放宽到 1000。

## 标本与患者层

新增每个内皮类型/门槛使用主 Q50 设计可测的全基因背景，逐标本对基因表达作平均并列秩，再计算 (rank−.5)/N；同一文库内以 counts 排名与其 CPM 排名一致。单端取均值，外部两端分别达覆盖一半方可取 UP−DN。去 MP10、去两 ATF4 程序共有、去 E2F/G2M 成员分别进行。原固定 mTOR NE 未更改。患者效应单位为排名分数 / 冻结 Epi 暴露 SD，非平台间原始分数比较。

OLS 模型按标本拟合，调整 dataset，经典同方差 t 区间；两个 singleton 导致 HC3 不适定，小样本区间依赖同方差假定。偏 Pearson、先秩化再调整 dataset 的偏 Spearman、留一患者斜率另存。敏感性包括 Q95、a/b 替换、≥30、单独调整深度(每细胞 counts)、周期(E2F/G2M 平均)、总内皮毛细比例。另以同标本 Epi_all 已有 M20 的 CHOL_SYNTH、FA_ELONG、mTORC1 作为连续暴露。主家族 7 程序、比较通路 4 个、外部 24 个双端评分分别 BH。不同内皮类型/门槛/遮罩/敏感性分开；三个 Epi 配对暴露在核心程序家族内合并 21 项。后续功能定义及竞争读出单列次要家族。无高低二分。

## 通讯、作图与复现

在 M0073 的 12 个标本中，按标本计算配体/受体阳性细胞比例，中位数≥10% 为可检出，候选使用原 NicheNet 人 v2 LR 与 ligand–target 先验。按四个实际 leading-edge 分别计算原先验 corrected AUPR，常数先验单列；Epi 合格排序与受体合格全域排序都保留。分母为 10,865 个可测且进入先验的基因。来源/受体表、既有双向 CellChat、反方向 M20 原排名和空间混合 spot 汇总均在 tables。真实 KO 与 BLCA 相关属于不同实验，患者层得分与 M0073 排名、通讯使用同一份 BLCA RNA。

单线程按 fetch.py、recover.py、literature.py、supplements.py 获取资源；GENCODE M25 表来自 UCSC mm10 wgEncodeGencodeAttrsVM25.txt.gz。随后依序运行 external.R（内含 maps.R）、project.R、scores.R、projection.R、communication.R（共享 flock）、predictions.R、figures.R。common.R 管理线程与依赖，fit.R 管理患者模型。已有任务请使用新目录副本运行代码以避免覆盖交付。绘图使用 jins-routine-plotting 的 Arial、组成模板和 jin_export，森林图、点图、方向图按数据类型单独制图。日志、文库审计、成员覆盖、模型状态、数据和图见同目录。

交付核验使用 verify.R、receipt.py。Full/M20 各有 50 个预设 Hallmark 条目、49 个可估条目。PANCREAS_BETA_CELLS 在 Full 只映射 11 个，低于原 15 个门槛；两版均以 NE 入表，其他通路保留各自原分析给出的 q。表中的 NE 与方向相反分开保存。图 05 的依赖方向相当于 CTRL−KD；原始 KD−CTRL 及其交互保存在外部基因统计表。患者排名评分中，百分位按 0–1 表示。

JLR Supplementary Data S2 的三张工作表已核对，ChIP 表保留 KLF6 启动子行供查证。QPC 的第二份补充 PDF 传输中断，附录读取状态和 SHA256 见 public_resource_manifest.tsv；五组实际表达输入、样本注释和所用文献基因全集均已读入。补充材料最终读取记录、签名来源、文件校验分别由 reference_checks.R 和 receipt.py 输出。
