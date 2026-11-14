# mp10_hypothesis_followup_v1

已执行的定向内部检验，属于BLCA MP10项目；没有使用LLPS/C5分组。先读RESULTS.md、HYPOTHESIS_UPDATE.md和FIGURE_INDEX.md。原结果只读保留。

结果表：model_results.tsv.gz（主与敏感性，所有NE）、model_registry.tsv（1473项，含36项条件临床分支和10项精确质量匹配对照）、model_membership.tsv.gz（实际输入和协变量）、sample_level_scores.tsv.gz（raw score与固定reference尺度）、coverage.tsv.gz、treg_cycle_composition.tsv。

原冻结计划见config/analysis_plan.yaml；临床条件分支、cluster字段修复及数值精度说明独立留档。旧统计见evidence_matrix.tsv和reused_*；不能混淆旧limma、GSEA与新HC3。输入原对象/pseudobulk/资源路径和hash见input_reuse_audit.tsv与original_source_input_manifest.tsv。原始数据未打包。

从项目根运行（R --vanilla，内部线程1，最多4总核）：

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_hypothesis_followup_v1/scripts/resume.py
```

只从结果表重画，无需原始Seurat、pseudobulk或重跑统计：

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 220 Rscript --vanilla mp10_hypothesis_followup_v1/scripts/09_replot.R
```

移至其他位置可设置`MP10_RESULTS_DIR=/absolute/path/to/results/mp10_hypothesis_followup_v1`，调用交付目录scripts/09_replot.R。依赖版本见logs/sessionInfo_*及config/environment_manifest.tsv。Arial、Cairo、ragg、ggplot2、cowplot用于原验收风格。

resume仅按阶段输入/成员/参数/代码/输出hash复用，任何不匹配会重算依赖步骤；不是仅看文件存在。M0冻结入口不在resume中自动重写。新增必要小缓存位于objects；大型原始对象与既有完整分析对象未复制。

每个图有完整PDF、PNG及必要的panel/legend PDF；plot_data为直接数据。qc中有数值/成员/FDR验证和渲染检查。module_status.tsv区分completed、not_evaluable、technically_blocked；NE不是阴性发现。
