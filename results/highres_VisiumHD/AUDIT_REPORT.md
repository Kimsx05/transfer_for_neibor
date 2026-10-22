# Visium HD 恢复、下载与输入审计

审计日期：2026-10-21（Asia/Shanghai）。工作目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD`。本轮只做输入恢复、计数/空间/影像审计和可测性描述；未运行分割、细胞注释、MP10 高低分组、差异分析或整合。

## 关键结论

| 数据 / section | 已审计分辨率 | 当前路线 |
|---|---|---|
| Cervilla / visiumhd_Bladder | 2 μm | bin 可用；原始 H&E 未定位，核级 pilot 仅有条件 |
| GSE300692 / D1_Johnson | 2、8、16 μm（同一 section） | bin 可用；需补原始高分辨 H&E 后再细胞重建 |
| GSE299573 / UC_SCC | 8 μm | 当前仅 bin；旧 derived 16 μm 不视为新的实测分辨率 |
| GSE318916 / MPBC1、MPBC2 | 各 2 μm | 表达与相对阵列网格可用；缺组织位置/影像/配准 |

所有 5 个 section 的 Top20 均需分别看 feature 与非零表达。逐基因精确计数和组织内分母见 `mp10_feature_detection.tsv`。MPBC1/2 的 in_tissue 未知，绝不以 filtered 文件名或相对网格替代组织掩膜。

## 恢复、迁移与预算

原始解压输入已实体移动到 `data/`；同文件系统 rename 保留 inode，不复制数据。旧入口保留兼容软链接；GSE299573 的四条失效输入链接已修复。迁移和修复日志为 `logs/relocations.jsonl`、`logs/link_repairs.jsonl`。历史分析结果未移动、删除或重算。

本地盘点先按文件名与旧 manifest 定位，并扩大到用户项目、生物数据及临时目录；完整搜索范围见脚本/日志。未找到项目根或本轮路径适用的 AGENTS.md。阅读了既有高分辨率审计的输入、结构、冻结与来源记录；其中历史“禁止联网”是旧任务范围，本轮按用户明确授权下载。

启动时 quota 约 458/500 GiB（硬上限 505 GiB），约 42 GiB 软配额余量；文件系统空闲约 12 TiB，/tmp 约 526 GiB；可用内存约 1.1 TiB。共享服务器 swap 已占满，因此使用最多 2 个数值线程、逐 bin 块稀疏读取，不生成密集全矩阵或完整 AnnData 副本。Cervilla 完整审计实测峰值 RSS 约 6.37 GiB（`logs/audit_complete_inputs.log` 的 Maximum resident set size），逐样本释放内存。本轮没有因磁盘资源阻断的数据集。

| 数据 | 本轮 data 内保留大小（MiB） | 获取方式 |
|---|---:|---|
| Cervilla | 337.36 | ZIP Range 选择性下载 |
| GSE300692 | 1589.18 | 移动既有输入 + 一次本地 TAR.GZ 遍历恢复小影像 |
| GSE299573 | 223.41 | 移动既有 spatial + 补下载 8 μm H5 |
| GSE318916 | 370.83 | 移动两份既有 H5；不重复下载 |

本轮记录到的大文件 HTTP 有效载荷总计 **491.32 MiB**（含 Range 探测和重复预读；不含 TLS/HTTP 头）。小型元数据另不足数 MiB；精确成功访问字节见网络日志及 metadata 文件。新下载最终保留 **483.87 MiB**。工作目录 data 合计 **2.462 GiB**，大部分为已有文件原地迁移。

D1_Johnson 的 27,503,110,814-byte TAR.GZ（25.614 GiB）保留旧位置引用；本轮仅从本地连续读取一次并消费到 gzip EOF，CRC 通过。未重新下载、未提取 BAM/大型可视化包。该归档不计入新工作目录的新增数据体积。临时下载峰值预算控制在约 1 GB 新增量；最终无压缩包副本或 `.part` 遗留。

Cervilla 官方 ZIP 为 3,362,702,746 bytes，官方 MD5 `e3f5e6dacae3b4c53f8007dabc64d9a6` 与用户线索一致。每个范围请求要求 HTTP 206 且 Content-Range 精确匹配；只提取必要成员，逐成员 CRC 和本地 SHA256 通过。**未下载完整 ZIP，因此没有验证整包 MD5。** 不获取 CytAssist_Bladder 数据集。

## 表达、空间与影像

| section / 分辨率 | GEX features × bins | counts 中位数 | genes 中位数 | barcode 匹配 | hires μm/像素 |
|---|---:|---:|---:|---:|---:|
| D1_Johnson / 2 μm | 18,085 × 6,362,240 | 8 | 8 | 100.0% | 5.591 |
| D1_Johnson / 8 μm | 18,085 × 400,147 | 152 | 140 | 100.0% | 5.591 |
| D1_Johnson / 16 μm | 18,085 × 100,828 | 625 | 526 | 100.0% | 5.591 |
| visiumhd_Bladder / 2 μm | 18,085 × 4,722,911 | 6 | 6 | 100.0% | 1.643 |
| UC_SCC / 8 μm | 18,085 × 695,646 | 157 | 126 | 100.0% | 1.231 |
| MPBC1 / 2 μm | 18,085 × 11,222,500 | 7 | 7 | NA | NA |
| MPBC2 / 2 μm | 18,085 × 9,781,598 | 1 | 1 | NA | NA |

按 100,000 bins/块读取 CSC；检查全体存储值、CSC 索引和 pointer，Gene Expression 单独统计；未发现负值、非有限值或非整数值。整数性仅为数值检查，实验 counts 归属由 GEO 说明、原始 Space Ranger 输出路径、H5 header 和报告共同支持。未使用 Thor 推断/扩散/插补表达。重复 symbol/ID、重复 barcode/坐标、零计数 bins、完整分位数与仿射系数均保存在矩阵审计表。

空间表按 barcode 重排连接，不依赖行顺序。array_row/array_col 为阵列索引，pxl_row/col_in_fullres 为 source-image 像素，叠图使用 fullres × tissue_hires_scalef。独立按 seed=123 取样拟合 grid→pixel 仿射，数值残差与镜像/旋转方向记录在表中；仿射行列式本身不作为错误判据。可视化只能排查明显错配，不能证明亚核级配准准确。

所有已定位 H&E 风格图均为处理后的 hires/lowres 或 alignment preview，未确认原始全分辨 H&E。D1_Johnson 的 web summary 明确引用 `JohnsonA1_D1_Scan1.qptiff` 和 `H1-PF3RWFN-D1-fiducials-image-registration.json`，但它们没有包含在当前归档中；CytAssist TIFF 是另一类影像，不能替代该原始 H&E。图像实际格式/宽高见 asset manifest 与 `metadata/image_inventory.json`。

Cervilla web summary 记录 Space Ranger 3.0.0、GRCh38-2020-A、Probe Set v2.0；原始输入路径指向 `Bladder_HD_he.tif` 与 `Bladder.json`，但这些文件未包含在 ZIP。Cervilla 的 6000×4711 hires，依据 scalefactors 为约 1.64 μm/像素，属于下采样图。仅在明确核形态仍可辨、配准通过人工复核的局部区域考虑有条件 pilot；优先找回原始注册 H&E。另有 1,274 个 in_tissue bins（0.02697%）位于提供的 hires 图下边界之外：图高 4711，组织内 y 最大 4879.343 pixels，约 168.343 pixels 超界。它们保留在 counts/MP10 统计中，未静默过滤；依赖影像的 pilot 应避开这些缺失视野并先核查原始影像。其余 99.9730% 的组织内 bins 在图像边界内。不得插值放大并声称恢复核细节。D1_Johnson 的约 5.59 μm/像素预览不作为合格分割输入。

ROI 用 ROI1/2/3 命名，固定 seed=123，覆盖不同坐标位置与边缘；不能确认上皮/淋巴/间质身份，因此未赋细胞标签。所有图仅作描述；MPBC1/2 只画 barcode 相对网格和 counts/genes 分布，不能画出可信 H&E 叠图。

## MP10 冻结与检测

冻结来源：`../MP10_spatial_multisection_top20_vs_core19_v1/01_signatures/01_MP10_Top20_vs_Core19_gene_table.csv` 的 DirectTop20，按 MP10_rank 排序；SHA256、期望名单和匹配结果见 `metadata/mp10_lock.json`。另用旧正式 `highres_ST/highres_spatial_audit_v3_GSE300692/01_MP10_lock/frozen_EpiMP10.tsv` 的 top20 交叉核对，其来源为正式 `ST/mapping/07_EpiMP_gene_basis.qs`；两处与用户给出的 20 基因完全一致，未替换基因、未修改定义。

| section / 分辨率 | feature 存在 | 所有输入 bins 实际非零基因 | 组织内实际非零基因 |
|---|---:|---:|---:|
| D1_Johnson / 2 μm | 20/20 | 20/20 | 20/20 |
| D1_Johnson / 8 μm | 20/20 | 20/20 | 20/20 |
| D1_Johnson / 16 μm | 20/20 | 20/20 | 20/20 |
| visiumhd_Bladder / 2 μm | 20/20 | 20/20 | 20/20 |
| UC_SCC / 8 μm | 20/20 | 20/20 | 20/20 |
| MPBC1 / 2 μm | 20/20 | 20/20 | 不可判定 |
| MPBC2 / 2 μm | 20/20 | 20/20 | 不可判定 |

组织内比例的分母是 supplied matrix 中按 barcode 匹配且 in_tissue=1 的 bins；空间表中所有组织内 bins 另列，避免把 filtered-matrix 分母当作整张 capture array。MPBC1/2 仅报告 all supplied bins 的检测和总计数，组织内列留 NA；候补数据不能因此宣布完成组织内可测性验证。

## 来源冲突、细胞结果和下一阶段

D1_Johnson：GEO 声称 Space Ranger 3.1.0；实际 H5 与 web summary 为 3.0.1。本报告保留两种来源，按实际输入版本审计。probe_set.csv 为 Visium Human Transcriptome Probe Set v2.0、GRCh38-2020-A。2/8/16 μm 的 filtered 集合略有不同，总计数差异不能直接解读为数据丢失。

MPBC1/2：GSE318916 标题是 HGT1 研究，但原始论文明确将两例 MPBC 标为额外的 MIBC 混合 micropapillary/urothelial 手术病例，不能沿用旧记录的 NMIBC 归类。论文称 Space Ranger 3.1.3；实际 H5 分别为 MPBC1=3.0.0、MPBC2=4.0.1，版本不一致须向作者核实。样本到患者的独立匿名 ID 未公开；跨数据集患者重叠/相邻切片关系 unknown。

论文报道 MPBC 使用 Bin2cell（mpp=0.5、prob_thresh=0.1、nms_thresh=0.5），随后做归一化/log/scaling。当前 GEO 仍只提供两份 2 μm H5；本轮未在数据开放说明定位作者 cell h5ad、核/细胞边界或 bin-to-cell map。若以后获取 h5ad，必须区分 X 的变换值与未经改变的聚合 counts，并检查 cell ID 一致性；扩张距离/软件版本未能从公开结果确认，记 unknown。

当前五个 section 均没有可直接复用且来源/映射可核验的细胞矩阵。Thor 论文使用推断表达，不等同实验 counts；当前 D1_Johnson 归档中没有该推断矩阵/细胞边界。checked 作者仓库仅有分析代码，不以“软件能生成”当作“已有结果可复用”。

官方 `spaceranger segment` 可独立对合格 H&E 做核分割，不需要先下载 FASTQ；但其图像核结果仍需与 bins 配准并聚合 counts。`spaceranger count` 重跑是另一路线，需要测序及完整实验/影像输入。本轮五个 section 均不具备本地完整 count 重跑输入。D1_Johnson 的人类原始测序未公开；GSE299573 已核实 SRX29144058 / SRR33944517，ENA 列出两端 FASTQ 合计 22,692,009,528 bytes，仅登记不下载。

建议下一阶段优先找回 Cervilla 和 D1_Johnson 的原始注册 H&E；取得合格输入后再决定 Bin2cell 或官方核分割+实际 counts 汇总。GSE299573 保持 8 μm bin 层，不能反推原始 2 μm。MPBC1/2 保留候补状态，优先请求作者 cell counts、边界和映射或原始配准输入。

## 访问受阻与交付检查

服务器默认大小写代理配置冲突导致最初 HTTPS EOF；随后显式使用现有本地代理成功获取 Zenodo/GEO。CORA Dataverse 镜像返回 HTTP 200 的反机器人 HTML，而非 JSON，记录为不可验证；未试图规避其人机验证。此限制不影响 Zenodo 原始 record 的核心成员获取。核心成员完成后，额外的 6,134,282-byte probe_set.csv 请求连续遇到 TLS EOF（包括官方文件直链复核），因此仅登记远端成员；版本已由 web summary 核实。没有将该小文件故障作为阻断其他数据的理由。Github 个别初次 TLS/分支错误已重试，成功和失败证据保留。未以某一个失效来源终止其他审计。

交付物：四个 TSV（asset_manifest、matrix_spatial_audit、mp10_feature_detection、segmentation_readiness）、本报告、figures 下轻量 PNG、code 下可复跑脚本，以及 metadata/logs 中来源快照、校验值和验证结果。没有大型报告包。`logs/validation.json` 记录表间交叉检查；`logs/software_versions.json` 记录复用 Python 3.12 环境。

复跑：在本目录执行 `bash code/run_audit.sh`。下载脚本只取本轮最小输入；ZIP 支持按已压缩成员断点恢复、CRC 后清理自身临时片段，GEO H5 按范围块续传；已有文件和审计检查点均按状态复用。若输入发生变化，审计签名失效后重查。图示抽样 seed=123。

## 核验来源

- [Cervilla Zenodo record](https://zenodo.org/records/17999961)；[作者代码仓库](https://github.com/scervilla/SpatialBenchmarking)。
- [D1_Johnson GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9066380)；[Thor 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12325965/)。
- [UC_SCC GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9040999)；[ENA accession](https://www.ebi.ac.uk/ena/browser/view/SRR33944517)。
- [GSE318916](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318916)；[MPBC 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13232474/)。
- [10x 官方 segment 输入](https://www.10xgenomics.com/support/software/space-ranger/4.0/analysis/running-pipelines/space-ranger-segment)；[Bin2cell 官方示例](https://bin2cell.readthedocs.io/en/latest/notebooks/demo.html)。
