# MP10 第一轮亚群筛选

本目录是本轮独立结果目录。阅读入口：`round1_report.md`（完成后生成）、`methods.md`、`review/`。执行范围仅模块 0–3。

- `tables/samples.tsv`：共同标本表、患者映射状态、MP10 原值与 dataset 内分组。
- `tables/proportions.tsv.gz`、`proportion_models.tsv`：四种分母及全部比较，含非显著/不可估计结果。
- `tables/stacked_bar_summary.tsv`、`stacked_bar_boundaries.tsv`：图 A 的样本/dataset 权重与累计边界。
- `tables/CellChat_*`：条件网络、原始细胞贡献和背景通讯；`objects/networks/` 仅保存可复用的小型数值网络。
- `tables/candidate_*`、`selected_axes.tsv`：发现、样本支持与最终候选（网络完成后生成）。
- `tables/circle_nodes.tsv`、`circle_edges.tsv`：图 B 的真实 LR 及共享布局/标尺。
- `figures/overall/`、`figures/by_dataset/`：透明 PDF 与 300 dpi PNG。比例总体图为共同 dataset 等权；CellChat ALL 为标本×亚群封顶后的 pooled 网络，两者不是同一种平均。
- `scripts/`：可重复运行的分模块脚本；`logs/`：实际运行记录。

脚本使用现有 R 4.5/4.4 用户库。Python 运行器为 `/usr/bin/python3`，脚本指向项目既有 Python 3.12 运行库；旧虚拟环境解释器软链接在当前主机不可用。无需重建环境。

顺序：01_samples.py → 02_extract_expression.R → 03_proportion_models.R / 04_cellchat.R / 05_abundance_figures.py → 06_model_sensitivity.R → 07_network_tables.R → 08_network_overviews.py → 09_candidate_support.R → 候选证据整合 → 11_circle_figures.R → 阅读包。CellChat 可复用已完成的同参数网络缓存。首次使用精确 triMean 加速前先运行 validate_fast_trimean.R。

患者数均为已有 derived/provisional 映射 ID 计数，不能冒称已临床核实；当前没有 verified patient ID。主要 MP10 比例结果须同时阅读原始 P、BH q、分期完整子集结果和单标本/映射患者敏感性分析。
