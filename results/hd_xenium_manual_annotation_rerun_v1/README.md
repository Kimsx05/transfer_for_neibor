# HD 与 Xenium 表达聚类人工注释包

优先打开 [人工注释工作簿](manual_annotation.xlsx) 和 [结果报告](REPORT.md)。工作簿的 HD、Xenium 两页分别使用各自平台共同表达图的 Leiden 0.5。先在样本分布、联合 UMAP 和代表性 marker 图上定位簇，再结合全基因表和样本内效应填写绿色列。正式大类和细类均由复核者填写。

| 内容 | 文件 |
| --- | --- |
| 原始 counts、行号定义、坐标、QC | [输入清单](input_manifest.tsv)、[过滤汇总](tables/qc_summary.tsv)、每平台的 `<sample>_qc.parquet` |
| 主分辨率 TSV 工作表 | [HD](hd/manual_review.tsv)、[Xenium](xenium/manual_review.tsv) |
| 优先复核线索 | [HD](hd/review_priorities.tsv)、[Xenium](xenium/review_priorities.tsv) |
| 对象级 UMAP、三个分辨率和来源 | [HD](hd/objects.parquet)、[Xenium](xenium/objects.parquet) |
| UMAP、空间图、marker 单图 PNG/PDF | [HD](figures/hd/)、[Xenium](figures/xenium/) |
| 全基因结果、样本内效应与 Top20 | 平台目录下各 `r0.2/`、`r0.5/`、`r1.0/`，见 [表行数与状态](tables/deg_inventory.tsv) |
| 重跑配置、环境、资源与校验 | [方法](METHODS.md)、[配置](tables/config.json)、[版本](tables/versions.json)、[进程树实测](logs/resources.tsv)、[验收表](tables/validation.tsv) |

对象关联请使用 `sample + original_cell_id` 或唯一 `cell_key`。零起始 `source_row` 只用于回查输入。按平台、分辨率、cluster 三列连接簇表；各平台和各分辨率分别按对象数从大到小编号。编号没有细胞类型含义，与旧 niche 编号也无映射。图中同一平台、同一分辨率的编号使用对应 `colors.tsv`。

每个分辨率保留对象级即时检查点 `labels.parquet` 及编号表 `cluster_key.tsv`。每平台只保存一套新模型、latent 和表达图；图中的行顺序在 `graph_rows.tsv.gz`。原始 counts 保持在输入位置。新 label 回填和后续子聚类可以通过对象表的键定位 counts。

在本项目相同目录层级的新空结果目录重跑时，复制本轮 `code/`，创建 `hd xenium tables figures logs` 子目录，执行 `bash code/start.sh model:hd graph:hd deg:hd plot:hd model:xenium graph:xenium deg:xenium plot:xenium finish:both`。这条命令提交带任务级硬限制的 user systemd service，需在同一服务器、原环境和数据挂载下运行。检查各 stages 的退出码和验证表。已完成模型受写入保护。

