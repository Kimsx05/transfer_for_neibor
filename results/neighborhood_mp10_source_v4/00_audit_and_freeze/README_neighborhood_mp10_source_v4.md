# 启动说明：病理检索 + 九个 niche 的中央 MP10 / 上皮归属 MP10

## 交给 Codex

将同目录的 `neighborhood_mp10_pathology_epi_codex_prompt_v4.md` 发给服务器 Codex，并附上：

```text
请完整读取并执行附件 neighborhood_mp10_pathology_epi_codex_prompt_v4.md。

本轮采用已经冻结的 SAMPLE_R12_LEIDEN，命名 niche0–niche8。
先实际检索可复用的病理区域标注，再完成中央 spot 的冻结 Top20/AUCell
和基于现有 cell2location 模型的上皮归属 MP10。病理未解决不阻塞评分。

不要重跑 niche、Harmony、CopyKAT 或全量 cell2location；不要只交付预检。
先在小批次验证表达归属实现，通过后继续全部模型可用切片。

完成九群结果、来源证据、具体病理待办和轻量回传包后停止。
```

## 执行顺序

M00 定位已有资产并冻结分析规格 → M01 本地与原始公开来源查找病理 → M02 中央整体 MP10 → M03 模型估计 Epi 表达归属和同背景评分 → M04 九群比较 → M05 可复用来源评估及有限混合对照 → M06 图表与回传。

M02/M03 可以在第一轮病理检索完成后并行；不等待全部 ROI 审阅完成。

## 回传后先看

优先查看 `RUN_SUMMARY.md`、`niche_MP10_summary_all9.tsv`、整体/上皮配对主图和 `source_evidence_all9.tsv`。病理问题另看 `WHAT_USER_NEEDS_TO_DO.md`，只处理 Codex 查找资料后仍未解决的项。

## 必须保留的区分

- `spot_MP10_frozen`：原冻结整体 spot 分数，不被覆盖。
- `spot_MP10_model_universe`：在模型覆盖的同一基因背景上重新计算的整体分数，用作配对对照。
- `epi_MP10_model_universe`：在上述同一背景上的模型估计上皮归属分数，不是实测纯上皮表达。

缺基因、缺模型、未知病理以及无法判断的来源均应报告。不能乘/除 Epi 比例来构造上皮分数，也不能只在 Top20 的20个基因内排名。

本轮结束后不自动启动通讯、机制分析或重新聚类。
