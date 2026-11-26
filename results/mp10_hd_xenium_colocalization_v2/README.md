# BLCA MP10–Treg 共定位与 niche，v2

本目录的初版为 provisional / Codex-reviewed，human_review_status=pending。HD 以 16 µm Treg 支持 bin 为目标，32 µm 上皮中心只用于尺度敏感性。Xenium 使用原分割细胞。

- [结果报告](REPORT.md)
- [分析方法](METHODS.md) · [冻结配置](config/frozen.json) · [输入清单](audit/input_manifest.tsv)
- [人工复核入口](review/index.html) · [填写说明](review/README.md)
- [逐片统计](reports/endpoints.tsv.gz) · [主检验](reports/primary.tsv) · [工具状态](reports/tools.tsv)
- [复现入口](code/rerun.sh) · [人工意见导入](code/import_review.py)

复核表导入在项目既有 Python 环境中运行，例如从 v2 目录执行：

```bash
bash code/import.sh review/Cervilla/human_review.tsv
```

初审表在 review/样本/candidates.parquet，人工表在 review/样本/human_review.tsv，各片标签、合格集合和冻结图在 hd/ 或 xenium/样本 下，niche 指派在 niche/样本 下。CD3 等 marker 逐项以实测计数或 NA 存储。人工签署前优先审阅 D1_Johnson 的混合 bin，以及目标较稀疏的 A1R4、A1R6。
