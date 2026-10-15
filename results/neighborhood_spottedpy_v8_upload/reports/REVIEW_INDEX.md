# v8 离线审阅入口

状态：`V8_COMPLETE_AWAITING_REVIEW`。先读完整报告，再核查全量证据和覆盖。

- [实际结果与全部目标](RUN_SUMMARY.md)
- [方法](METHODS.md)
- [限制](WHAT_IS_STILL_UNCERTAIN.md)
- [续跑](RESUME.md)
- [自动验证](validation_results.tsv)
- [九niche证据](../07_synthesis/niche_geometry_evidence_all9.tsv)
- [29类细胞证据](../07_synthesis/cell_geometry_evidence_all29.tsv)
- [所有dataset证据](../07_synthesis/dataset_evidence_all.tsv)
- [全部不可推断结果及原因](../07_synthesis/unavailable_results.tsv)
- [点级分区索引](../02_hotspots/point_hotspot_statistics_manifest.tsv)
- [全部热点覆盖](../02_hotspots/hotspot_coverage_all.tsv)
- [两层匹配结果](../05_sensitivity/two_score_layers_matched_units.tsv)
- [匹配敏感性](../05_sensitivity/matched_sensitivity_effects.tsv)

点级原始parquet留在服务器，不放入轻量包；索引所列相对path是服务器项目根目录相对路径。

## 总体图

图上星号为关系检验q_scope。局部热点P/q只用于定义区域。
- [CELL_MP10_EPI_LOCAL_BH05_heatmap](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_heatmap.pdf)
- [CELL_MP10_EPI_LOCAL_BH05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_k_sensitivity.pdf)
- [CELL_MP10_EPI_LOCAL_BH05_unit_boxes](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_LOCAL_BH05_unit_boxes.pdf)
- [CELL_MP10_EPI_SOURCE_P05_heatmap](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_heatmap.pdf)
- [CELL_MP10_EPI_SOURCE_P05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_k_sensitivity.pdf)
- [CELL_MP10_EPI_SOURCE_P05_unit_boxes](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_EPI_SOURCE_P05_unit_boxes.pdf)
- [CELL_MP10_SPOT_LOCAL_BH05_heatmap](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_heatmap.pdf)
- [CELL_MP10_SPOT_LOCAL_BH05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_k_sensitivity.pdf)
- [CELL_MP10_SPOT_LOCAL_BH05_unit_boxes](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_LOCAL_BH05_unit_boxes.pdf)
- [CELL_MP10_SPOT_SOURCE_P05_heatmap](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_heatmap.pdf)
- [CELL_MP10_SPOT_SOURCE_P05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_k_sensitivity.pdf)
- [CELL_MP10_SPOT_SOURCE_P05_unit_boxes](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/CELL_MP10_SPOT_SOURCE_P05_unit_boxes.pdf)
- [NICHE_MP10_EPI_LOCAL_BH05_heatmap](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_heatmap.pdf)
- [NICHE_MP10_EPI_LOCAL_BH05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_k_sensitivity.pdf)
- [NICHE_MP10_EPI_LOCAL_BH05_unit_boxes](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_LOCAL_BH05_unit_boxes.pdf)
- [NICHE_MP10_EPI_SOURCE_P05_heatmap](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_heatmap.pdf)
- [NICHE_MP10_EPI_SOURCE_P05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_k_sensitivity.pdf)
- [NICHE_MP10_EPI_SOURCE_P05_unit_boxes](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_EPI_SOURCE_P05_unit_boxes.pdf)
- [NICHE_MP10_SPOT_LOCAL_BH05_heatmap](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_heatmap.pdf)
- [NICHE_MP10_SPOT_LOCAL_BH05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_k_sensitivity.pdf)
- [NICHE_MP10_SPOT_LOCAL_BH05_unit_boxes](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_LOCAL_BH05_unit_boxes.pdf)
- [NICHE_MP10_SPOT_SOURCE_P05_heatmap](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_heatmap.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_heatmap.pdf)
- [NICHE_MP10_SPOT_SOURCE_P05_k_sensitivity](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_k_sensitivity.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_k_sensitivity.pdf)
- [NICHE_MP10_SPOT_SOURCE_P05_unit_boxes](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_unit_boxes.png) · [PDF](../06_figures/ALL_DATASETS/summaries/NICHE_MP10_SPOT_SOURCE_P05_unit_boxes.pdf)

## 七dataset技术选择的空间overview

- [CanDis__UNVERIFIED：CanDis_2026_Sample9](../06_figures/CanDis__UNVERIFIED/sections/CanDis_2026_Sample9/MP10_overview.png) · [PDF](../06_figures/CanDis__UNVERIFIED/sections/CanDis_2026_Sample9/MP10_overview.pdf)
- [GSE171351：GSE171351_B1](../06_figures/GSE171351/sections/GSE171351_B1/MP10_overview.png) · [PDF](../06_figures/GSE171351/sections/GSE171351_B1/MP10_overview.pdf)
- [GSE238145：GSE238145_B14](../06_figures/GSE238145/sections/GSE238145_B14/MP10_overview.png) · [PDF](../06_figures/GSE238145/sections/GSE238145_B14/MP10_overview.pdf)
- [GSE285715：GSE285715_T2](../06_figures/GSE285715/sections/GSE285715_T2/MP10_overview.png) · [PDF](../06_figures/GSE285715/sections/GSE285715_T2/MP10_overview.pdf)
- [GSE319536：GSE319536_BC01](../06_figures/GSE319536/sections/GSE319536_BC01/MP10_overview.png) · [PDF](../06_figures/GSE319536/sections/GSE319536_BC01/MP10_overview.pdf)
- [Inhousedata__UNVERIFIED：Inhousedata_T4469479](../06_figures/Inhousedata__UNVERIFIED/sections/Inhousedata_T4469479/MP10_overview.png) · [PDF](../06_figures/Inhousedata__UNVERIFIED/sections/Inhousedata_T4469479/MP10_overview.pdf)
- [NPJ__UNVERIFIED：NPJ_sample6](../06_figures/NPJ__UNVERIFIED/sections/NPJ_sample6/MP10_overview.png) · [PDF](../06_figures/NPJ__UNVERIFIED/sections/NPJ_sample6/MP10_overview.pdf)

## 各dataset完整图包

- `figures_CanDis__UNVERIFIED_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_GSE171351_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_GSE238145_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_GSE285715_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_GSE319536_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_Inhousedata__UNVERIFIED_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
- `figures_NPJ__UNVERIFIED_spottedpy_v8.zip`：解压后打开其 REVIEW_INDEX.md 或 INDEX.html；含该dataset全部切片/29目标及总体关系页。
