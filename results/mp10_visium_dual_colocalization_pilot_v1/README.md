# 普通 Visium 四片：分开映射、实测功能与原 niche

[报告（HTML）](REPORT.html) · [报告（Markdown）](REPORT.md) · [独立图件](figures.html) · [详细方法](METHODS.md)

实际结果目录为 `/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_visium_dual_colocalization_pilot_v1`。从项目根目录使用 `bash results/mp10_visium_dual_colocalization_pilot_v1/run.sh` 可顺序重建。默认在本结果目录内重新计算并写入。若输入及参数完全相同，可设 DUAL_RESUME=1 复用已完成的逐片映射。输入未变时亦可直接用已保存的 spots 主表和统计表重画。修改输入或参数应另建版本目录。当前报告结论是本次冻结输入的解读。

| 入口 | 内容 |
|---|---|
| input_manifest.tsv、frozen_query.tsv.gz、tables/join_audit.tsv | 只读数据来源、实际 query barcode 与原 niche 覆盖 |
| reference_manifest.tsv.gz、reference_sampling.tsv、annotation_mapping.tsv | 每个参考细胞 ID、来源与精确类别、抽样前后数量、29 个细亚群及大类 |
| gene_sets.tsv、tables/gene_coverage.tsv、features/ | 完整基因成员、可测与检出、Top20 重叠和实际 PCA/anchor 特征交集 |
| spots/切片名.tsv.gz | 每个全片 spot 的数值主表及坐标。新上皮分数 H05_separate_transfer_score，旧 H05_support 与 H05_conditional_support 作为历史对照 |
| spots/切片名.Epi.tsv.gz、spots/切片名.TME.tsv.gz | 两参考原始类别支持度。TME_score__ 前缀为 Seurat 分类支持度，所有 29 类之和为 1 |
| tables/A_pairwise.tsv.gz、tables/B_pairwise.tsv.gz | 每片完整对称矩阵的非重复上三角，含有效 n、原始 r/ρ 和残差 r/ρ；表的左右两列是轴名 |
| tables/four_way.tsv、tables/endpoint_swaps.tsv.gz、tables/endpoint_swap_summary.tsv | 四组合的同一总体效应、两端替换差异和类型谱一致性。每种组合含细亚群与共同大类 |
| tables/abundance_sensitivity.tsv.gz | 绝对均值、占 all30 比例和占非上皮 29 类比例；主 log1p 均值在 four_way 表 |
| tables/functional_correlations.tsv、tables/functional_regions.tsv、tables/functional_contrasts.tsv | 两种上皮读数与实测功能关系、片内 Q25/Q75、真实分组人数、组差 |
| tables/functional_TME.tsv.gz、tables/discordant_expression.tsv | 直接功能—细胞关系及两个相反尾部的原始 UMI、组织位置、深度 |
| tables/niche_scores.tsv、tables/peripheral_vs_center.tsv.gz | 各原 niche 的读数、相同中心上的同 spot 和外围 R1/R2 对照 |
| tables/old_vs_new.tsv、tables/low_TME.tsv、tables/score_QC.tsv、tables/dynamic_ranges.tsv、tables/diagnostic_covariates.tsv | 旧对照、低丰度诊断、范围与类别和、组成与深度关系 |
| figure_index.tsv、figures/、figures.html | 正式独立 PNG/PDF，图例分离版本随相应主图保存 |
| logs/、resources.tsv、status.tsv、config.json、versions*.txt、code/ | 运行日志、峰值 RSS、交付核验、参数、环境和复现代码 |

逐 spot 表中 c2l_mean__ 是原始 posterior mean absolute abundance，log1p_c2l_mean__ 是主分析变换。major_ 开头的类型为原始均值先求和，或子类映射分数直接求和。Synthesis、Homeostasis、Fatty_acid、Efflux 各自带 full 或 minus20 后缀；minus20 仅指从功能基因集去掉重叠成员。niche 按原色板显示为 1–9，accepted_label 保留上游编号。未分配不纳入各 niche 的分母。报告页面使用相对链接，复制查看时保留目录结构。


补充文件中，tables/functional_detection.tsv 给出每片 query 的实测成员稀疏度与零分比例，spots/*.original_rings.tsv.gz 保存同一中心的 q05/all30 及原 R1/R2 环比例，tables/original_niche_profiles.tsv 保存原标签的全队列画像。完整 spot 表也保留 mean_all30_fraction__ 和 mean_TME29_fraction__ 两种不同分母。
