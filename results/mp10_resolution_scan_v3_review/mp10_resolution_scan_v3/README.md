# MP10 resolution scan v3

首先打开 [REPORT.html](REPORT.html)。目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3`。

24套固定FULL_RAW/FULL_Z分区，38个固定候选，108个固定对照（105个精确模型键）。状态：{'COMPLETE': 103, 'NOT_TESTABLE': 5}。旧分区成员复现通过。两种high及候选/局部对照均在DE前冻结。

|候选|细胞|NMF P/R|UCell P/R|展示marker|
|---|---:|---|---|---|
|FULL_RAW_r0.40__C4|12585|28.02%/28.50%|34.01%/34.59%|无符合规则的基因|
|FULL_Z_r0.80__C13|6249|35.27%/17.82%|62.52%/31.58%|IFRD1;NEDD9;EZH2|

没有保证找到“漂亮群”。结果显示high定义依赖、来源集中、零召回和成员边界变化；所有marker/功能均是同批发现RNA，不是外部验证。样本不是已确认独立患者。

主要表：[全部候选](08_review/FINAL_CANDIDATES.tsv)、[24套概览](08_review/24_partition_overview.tsv)、[模型状态](04_pseudobulk_markers/model_summary.tsv)、[完整短名单](08_review/all_shortlists.tsv)、[功能摘要](08_review/function_family_summary.tsv)、[图索引](07_figures/figure_index.tsv)。完整DE和过滤在04_pseudobulk_markers/models，每个比较通过comparison_model_index.tsv追溯；功能全表按同model_key在05_function。

复现：code/01_audit_partitions.R后执行code/run_remaining.sh；出图/报告/验收脚本见code。相同配置/精确成员的缓存可复用。原矩阵/Seurat未打包，轻量审阅ZIP包含完整DE/功能表、短表、图、代码和说明。
