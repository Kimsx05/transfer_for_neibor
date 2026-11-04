# MP10逐切片Pearson分析

本目录包含全部63片实际计算，主分析为SCT data AddModuleScore + 原始逐spot Pearson；主MP10仍为冻结AUCell。

- `full_overview.pdf`: 全队列总览。
- `RESULTS.md` / `summary.md`: 中文解释、实际数值与限制。
- `full_sections.tsv`: 每片身份、输入、覆盖及模块状态。
- `reports/<section>.pdf`: 每片3页主报告；同名MD为覆盖及具体关系。
- `atlas/<dataset>.pdf`: 全部切片与同dataset样本比较。
- `figure_index.tsv`: 主图、逐片图位置；各图有300dpi PNG和有图例时的三个PDF。
- `tables/pearson_all.parquet`: 完整85×29×两评分×各人群/范围长表。
- `tables/delta_all.parquet`: full_range描述、common_support正式比较、0–1000共同支持。
- `tables/MP10_candidate_all.tsv.gz` / `closeness_all.tsv.gz`: 比较四/一。
- `tables/context85_effects.parquet` / `context85_hierarchical_summary.tsv`: 旧80+新5的UCell上下文模型，β不等于Pearson。
- `tables/unit_correlations.tsv.gz` / `dataset_heterogeneity.tsv.gz`: 真实身份层级与空间不确定性异质性。
- `tables/front_only_niche_all.tsv.gz`: 仅前沿带新增九类niche比较。
- `tables/dataset_priorities.tsv` / `section_relationship_priorities.tsv.gz`: 覆盖和效应驱动的审阅建议。
- `tables/program_gene_coverage.tsv.gz` / `program_symbol_mapping.tsv.gz`: 实际成员、别名、检测与SCT覆盖。
- `spots/<section>.parquet`: 冻结身份、三距离轴、29候选、85 UCell及85 AMS绘图源数据，无重复表达矩阵。
- `scores/`: 紧凑评分与参数；不保存大型SCT对象。
- `scripts/run_all.sh`: 可重启入口；`scripts/redraw.py`: 任意候选/程序散点。

所有普通Pearson P/q都有独立spot假设限制；空间区间/P另列。灰色NA不等于0。相邻正常与肿瘤来源分开。完整Parquet可用pandas/pyarrow读取，主程序与敏感性在program_family列明确区分。review.zip只放说明、总览、关键小表和数据集森林图；大型全组合表、逐点缓存、表达数据不放入包。
