# BLCA Visium HD 与 Xenium 首轮分析

[结果报告](REPORT.md) · [方法与阈值](METHODS.md) · [逐分支状态](reports/branch_status.tsv) · [输入清单](shared/input_manifest.tsv) · [统一 taxonomy](shared/taxonomy.tsv)

| 位置 | 内容 |
|---|---|
| hd/切片 | 16 µm 单元、原始 8 µm 成员、RCTD 权重、冻结区域、50 µm 图及独立区域邻接 |
| xenium/样本 | H5 与分割核验、互斥身份、panel 和基因审计、上皮状态、50/200 µm 空间图 |
| shared | 参考成员、共同基因、统一颜色、软件版本、MSigDB 原始 JSON、留出验证 |
| figures/样本 | 独立 PNG/PDF。矩阵图色条为 Z-score，主报告不使用对角线做生物学判断 |
| reports | 完成状态、逐样本效应、计数、验证及方向对应表 |
| code、logs | 重跑代码、各轮输出、版本及故障记录 |

从项目根执行 `bash results/mp10_hd_xenium_colocalization_v1/code/rerun.sh`。已完成缓存按样本复用。原始矩阵路径由 manifest 标识。Xenium 矩阵按 GEO 的细胞标识取交集并核验逐细胞计数。稀疏图行列均按 graph_nodes.parquet 的 id 列排序；Epi64 的矩形图行顺序见 epi64_centers.parquet，列顺序见 HD graph_nodes.parquet。counts 子目录中的 CSR 由 code/base.py 的 read_export 读取，按单元×基因排列。分割和表达均保留原始 ID，跨样本连接时必须同时指定样本。

区分进度、缺测与零值。fit_complete.txt 仅在整片原生 RCTD 全部批次完成后写入。受支持标签见 labels.parquet；低信息或混合单元保留。MPBC1 和 MPBC2 仅登记。

逐项执行和不能复用作者终标签的影响见 [完成清单](reports/completion.tsv)。读者可直接打开[图件索引](reports/FIGURES.md)和[首轮判读](reports/INTERPRETATION.md)。矩阵补充文件 matrix_full_taxonomy.parquet 同时列出在当片未出现的预设类型、0 个节点及不可估计项。
