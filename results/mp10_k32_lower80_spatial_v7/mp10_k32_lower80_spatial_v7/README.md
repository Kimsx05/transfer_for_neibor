# K32 lower80 spatial v7

打开 REPORT.html 查看中文总览。轻量上传目录：07_review；最终包：07_review/review.zip。服务器完整模型路径及续跑命令见 MODEL_INDEX.tsv。

本轮已授权：新31类NB参考、精简DEV、两张真实切片；旧v5/v6只读。本轮训练完成与数值稳定、生物学证据分开记录。

## 实际模型状态

- NB_LOWER80: 250 epochs；STABLE；1000 draws；生物学证据：部分支持。
- GSE171351_B1: 30000 epochs；PROVISIONAL_NOT_CONVERGED；1000 draws；生物学证据：不可评价。
- DEV: 9000 epochs；STABLE；1000 draws；生物学证据：部分支持。
- GSE238145_B17: 29000 epochs；STABLE；1000 draws；生物学证据：不可评价。

## 阅读入口

- [中文总览](REPORT.html)
- [实际方法](METHODS.md)
- [模型路径与命令](MODEL_INDEX.tsv)
- 核心已知真值结果：05_evaluation/DEV_recovery_metrics.tsv；逐来源、低背景及完整29类分配表在同目录。
- 空间PDF/300dpi PNG、英文图注和plot_data：06_figures/by_dataset。
- 实际耗时、稳定性与生物学结论分别见STATUS.json及中文报告；不以技术完成代表生物学通过。
