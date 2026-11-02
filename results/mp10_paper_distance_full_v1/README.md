# Epithelial MP10 全量空间距离分析

[全量总览 PDF](full_overview.pdf) 收录覆盖、逐片和数据集效应、固定调整、匹配评分层、分区与功能结果及七个预先登记的空间示例。[结果](RESULTS.md) · [方法](METHODS.md) · [逐片名单和可评估覆盖](full_sections.tsv) · [图索引](figure_index.tsv) · [审阅包](review.zip)

## 图册与逐点数据

[Dataset 图册索引](tables/atlas_index.tsv) 指向完整矢量 PDF。图索引提供 atlas 页码、每张完整 PDF、白底 300 dpi PNG、主体和图例文件。每片 9 页，依次为组织背景、三个目标、连续 MP10 与实距、分区组成与 R12、两页重点程序图、Treg 全组织背景、距离分箱程序走势、全部程序效应。总览按来源单列邻近正常。

spots/<section>.parquet 为逐片整合表，含完整 ID、原图坐标、μm 坐标、v3.4 几何、身份、MP10 原评分及有效标记、30 项组成与 mean、各 target/距离/原因、R12 和 dataset-local 标签、标准转录程序。scores/ 包含可复用的逐片评分缓存。服务器保留全图和缓存；review.zip 提供报告、总览、核心表及七张代表片的连续距离图。

## 数据表

| 内容 | 入口 |
|---|---|
| 登记、模块覆盖及不可用原因 | [源文件](tables/input_manifest.tsv)、[分来源覆盖](tables/cohort_coverage.tsv)、[不可估计效应明细](tables/unavailable_effects.tsv) |
| 目标定义、定位、距离和查询子集 | [阈值和绝对丰度](tables/targets.tsv)、[分区与 niche 定位](tables/target_localization.tsv)、[四个距离条件](tables/associations.tsv)、[距离覆盖](tables/distance_coverage.tsv)、[实距分箱](tables/distance_bins.tsv)、[分区内效应](tables/context_distance.tsv)、[anchor 子集](tables/anchor_subset.tsv) |
| 新增的固定对照 | [完整配对集调整前后](tables/adjusted.tsv)、[同一配对集的模型评分层](tables/matched_layers.tsv) |
| 分区与 niche | [区域组成及 MP10](tables/region_composition.tsv)、[两种 niche 的分母、覆盖、区域分布与 Epi-query MP10](tables/region_niche.tsv)、[原 R1/R2 profile](tables/source_R12_profiles.tsv) |
| 功能 | [锁定基因集](tables/gene_sets.tsv.gz)、[逐片匹配](tables/gene_set_coverage.tsv)、[55 项距离效应](tables/program_distance.tsv)、[分区功能、按统计对象](tables/region_programs.tsv)、[胆固醇一致性](tables/cholesterol_MP10.tsv)、[Treg 全组织参照](tables/Treg_program_composition.tsv)、[原两片 AMS 对照](tables/method_comparison.tsv) |
| 身份和跨数据集汇总 | [切片长表](tables/section_effects.tsv.gz)、[单位级效应和成员切片](tables/unit_effects.tsv.gz)、[按单位/患者等权](tables/unit_equal.tsv)、[dataset 分布和有效数](tables/dataset_effects.tsv.gz)、[数据集等权](tables/dataset_equal.tsv) |

汇总长表按 source_context、family、condition、target、feature 和 statistic 区分统计量。dataset 等权表的 population 区分全部切片、可用单位、核验患者。在 dataset 效应表中，all_sections 的 dataset_effect 为切片中位数；analysis_units 和 verified_patients 为单位等权 Fisher-z 均值。未估计效应在逐片表保留 NA；跨单位表仅汇总实际可用值。

## 复用

在项目根目录运行 `bash mp10_paper_distance_full_v1/scripts/run.sh`。默认从冻结 manifest 动态发现切片。可通过 MP10_PROJECT、MP10_MANIFEST、MP10_OUTPUT 指定项目、输入清单和输出目录。只在缓存对应同一套冻结输入和参数时复用；新的输入版本应使用新的输出目录。成员文件直接来自原 pilot，无联网更新。完整运行需要同一项目数据、R 用户库和定稿绘图 skill。配置见 [参数](parameters.tsv)。
