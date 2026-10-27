# MP10 round3：全TME总通讯

本轮读取v2的12个完整共享状态网络，完成模块0–4。主图是全部LR聚合的全细胞群网络；B/Plasma、T/NK与旧AX列表均不限制范围。v1/v2保持只读。

阅读入口：

1. [中文报告](round3_report.md)：全TME格局→上皮双向总量→Treg定位→其他群→LR贡献→NECTIN2–TIGIT补充。
2. [方法](methods.md)：W/C/R、自环、NA、排名与覆盖。
3. [逐队列图索引](FIGURE_INDEX.md)：ALL及五队列，major/fine完整输出。
4. [轻量阅读包](review.zip)：核心PDF、主要表和报告，不含PNG、旧网络或表达矩阵。

`tables/object_index.tsv`记录来源与参数；`total_node_edges.tsv.gz`是完整source×target；`global_incoming_outgoing.tsv`、`epi_allTME_total_contrasts.tsv`、`cross_dataset_direction_summary.tsv`是全群主结果；Treg及单轴补充从同一全量结果提取。`objects/aggregated_networks.qs`只包含新增的小型W/C/R矩阵。

高低指**标本内相对epi状态**，共享对象只有一个Treg节点。ALL不算独立队列。记录患者ID沿用derived/provisional状态。CellChat内部P不是状态差异P，总强度不能直接解释为免疫抑制。

旧Normal/Tumor、T1/T2和标本MP10分组背景见[round1 v2](../mp10_subtype_screen_round1_v2/README.md)，本轮没有重复推断或重画这些背景。

脚本顺序见 `scripts/run_round3.sh`；使用已有R/Python环境，不升级依赖。
