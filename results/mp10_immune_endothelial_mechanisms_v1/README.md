# MP10 immune and endothelial mechanisms v1

本目录属于膀胱癌 **MP10** 项目。分析使用单细胞全转录组 RNA，不使用 LLPS/C5-CAF 分组，也不把既有 pooled CellChat 的内部 P 当作患者差异检验。

代码与结果分离：

- 代码：`../../mp10_immune_endothelial_mechanisms_v1/scripts/`
- 项目独立 R 库及 Python 环境：`../../mp10_immune_endothelial_mechanisms_v1/{Rlib,venv}/`
- 共用冻结参考：`../../resources/mp10_regulatory_priors/`
- 本目录：表、紧凑 pseudobulk、必要推断缓存、图、配置和日志。

## 输入与定义

原始对象在 `subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs`，就地读取，不复制。原 RNA counts 和 normalized data 层及 barcode/feature 对齐分别核验。输入来源、大小、时间、哈希见 `tables/input_manifest.tsv`；完整注释/细胞成员、评分及伪批量成员均有独立表。

使用最终 c2l 30 类的 exact reference membership、基线 run `c2lfinal_20261108T123140Z_0b5fd64f` 的 40 个原始合格 Tumor specimen、冻结 MP10 Top20/MP10_cell。每 specimen 全部 Epi 的 R type-7 分位数定义 H50/H10/H05 和共同 L50，ties 不拆分。原始 eligible Epi 为 116,448 个。高分尾部不自动等于真实活跃细胞。

主 pseudobulk 门槛为 **>=50**，敏感性 **>=30**；使用全部细胞而非 CellChat capped membership。T4463209a/b 是已确认同患者不同取材，主模型固定 a，敏感性替换 b。其他记录/派生 patient ID 的独立性不视为已经全部确认。

## 主要交付

本轮已完成全部可运行模块：81 个预设 DE 模型中 68 个完成、13 个明确 NE；三条谱系的 pySCENIC 全量评分完成。50 张图共 198 个图文件；651 项核心检查和 37 项 SCENIC 检查通过。实际发现及不能建立的机制见 RESULTS.md。


- `METHODS.md`：冻结设计、模型、阈值、先验及检验家族。
- `RESULTS.md`：中文实际发现与限制。
- `FIGURE_INDEX.md`：图清单、直接绘图数据及统计含义。
- `tables/patient_relationship_audit.tsv`：源 patient ID 状态与用户确认的 a/b 重复关系分列，不把其余派生 ID 当成已验证独立患者。
- `tables/model_plan.tsv`、`model_status.tsv`、`model_membership.tsv.gz`：全部模型、实际覆盖与 NE 原因。
- `tables/de/M####.tsv.gz`：全部可检验基因，不只候选。
- `tables/activity_associations.tsv.gz`、`profile_regulatory_activities.tsv.gz`、`gsea_results.tsv.gz`：功能/调控结果与各自 P/q。
- `tables/nichenet_analysis_status.tsv`、`nichenet_ligand_activities.tsv.gz`、`nichenet_targets.tsv.gz`、`nichenet_backgrounds.tsv.gz`、`nichenet_target_links.tsv.gz`：双向先验、目标、背景、完整活动排序。
- `tables/scenic_*`：共享网络、regulon、全细胞 AUCell、样本层关联、RSS、第二种子关键连接复核。
- `tables/complete_candidate_evidence.tsv.gz`：保留所有可审查 LR 候选的分项证据，不使用任意加权总分。
- `module_status.tsv`：completed / not_evaluable / technically_blocked 及范围。
- `tables/validation_checks.tsv`、`qc/`、`logs/`：一致性与渲染核查。

`not_evaluable`、feature absent、低检出、观察到的零表达分别记录；缺失统计量不是零效应。NicheNet 的 activity 是排序，不是患者层 P/FDR。ULM 内部 enrichment P 与样本层模型 P 分开。多设置、同一矩阵上的多个算法不当作独立验证。

## 重画与复现

所有命令从项目根目录运行：

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
bash mp10_immune_endothelial_mechanisms_v1/scripts/replot_all.sh
```

重画仅依赖 `tables/` 与固定绘图配置，不读取原始 Seurat、全量表达或重跑 CellChat。白底、Arial 约 12 pt、矢量 PDF、300 dpi PNG；有图例时同时输出完整图、主体 PDF、图例 PDF。直接数据在 `plot_data/`。

完整复现/断点执行入口：

```bash
bash mp10_immune_endothelial_mechanisms_v1/scripts/run_all.sh resume
# 重新建立 pseudobulk 只在明确需要从原始输入重建时执行：
bash mp10_immune_endothelial_mechanisms_v1/scripts/run_all.sh rebuild_pseudobulk
```

完整启动器按模块顺序执行；GRN/ctx 的 worker 数由冻结配置或显式运行参数记录，AUCell 与关键连接复核每谱系最多 2 个单线程 worker。并行执行谱系时，全部本任务计算 worker 合计不得超过 4；已运行的 Epi/T 谱系各固定两个 CPU。首次生成三套 GRN 的输入、seed123、软件及输出哈希在各 lineage metadata 中记录。缓存使用实际输入/参数哈希；不以同名文件存在替代验证。

本机全局 `Rprofile.site` 会重设线程数为 24。复现命令使用 `Rscript --vanilla`，`common.R` 另显式重置 BLAS/OMP/RcppParallel/data.table。早期运行发现的线程约束偏差及修正见 `logs/resource_constraint_incident.json`，不隐瞒该执行偏差。原始数据和既有 HD/CellChat 作业未改动。

不提供包含原始大对象的压缩包。共享数据库无需在每次重画时下载。安装/版本与资源来源、哈希保存在 `config/`、`logs/` 和共享参考 manifest；不自动升级环境。

The complete excluded-feature list is `tables/de_gene_filter_exclusions.tsv.gz`: observed zero UMI across a model differs from nonzero expression failing filterByExpr. Full per-profile positive fractions are exactly recoverable from `objects/pseudobulk.qs` as detected-cell counts divided by profile cell counts. `tables/nichenet_ligand_activities_annotated.tsv.gz` adds actual, model-specific sender sources; `nichenet_target_mapping_audit.tsv` separates no DE targets from targets absent in the prior.
