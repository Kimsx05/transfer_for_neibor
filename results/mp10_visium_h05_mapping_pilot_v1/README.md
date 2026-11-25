# 普通 Visium 四片 H05 映射

[打开图文报告](REPORT.html) · [纯文本报告](REPORT.md) · [方法](METHODS.md)

在项目的 results 目录中独立运行，绝对路径为 `/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_visium_h05_mapping_pilot_v1`。从项目根目录运行 `bash results/mp10_visium_h05_mapping_pilot_v1/run.sh`。仅使用已有环境。运行入口会刷新本目录的分析文件；上游只读。独立单图索引在 [figure_index.tsv](figure_index.tsv)。

| 文件 | 内容 |
|---|---|
| input_manifest.tsv、barcode_audit.tsv、coordinate_score_audit.tsv | 原始输入位置、大小与时间戳，两次逐键核验和冻结 query 数 |
| frozen_definitions.tsv、frozen_specimen_thresholds.tsv、membership_version_difference.tsv | 原标签定义、40 个 specimen 的原 Q95 及一枚边界差异 |
| reference_sampling_ids.tsv.gz、reference_sampling_proportions.tsv | 参考细胞和来源、原大类、状态和抽样前后比例 |
| config.json、versions_R.txt、versions_python.txt、features/ | 固定配置、版本、各切片及分支真正用于 PCA 和 anchors 的特征 |
| spots/切片名.tsv.gz | 全片主表，以 section_id 和完整 barcode 为键；所有大类支持度、H05 原始/条件支持度、Epi 分母、旧评分、位置、marker 及 Top20 原始 UMI |
| spots/*.noTop20.tsv.gz | 每个 spot 消融后所有类别的支持度，同一主键 |
| section_metrics.tsv、common_support_comparison.tsv、composition_depth_diagnostics.tsv | 片内相关、缺失、低分母、共同覆盖与组成/深度调整 |
| feature_ablation.tsv | 相同可用点上的去 Top20 对比 |
| denominator_diagnostics.tsv、non_epi_diagnostics.tsv、label_context.tsv | 低分母、c2l 优势谱系与 marker 回看 |
| display_thresholds.tsv、discordant_regions.tsv、discordant_expression.tsv | 展示阈值、六邻接连通域和原表达；原始支持度最高 5% 仅作描述 |
| pseudo_spot_validation.tsv、pseudo_paired_depth.tsv、pseudo_spot_predictions.tsv.gz、pseudo_truth.tsv、pseudo_members.tsv.gz | 旧 TEST 的重标真值、每个细胞身份、每个 mixture 的真实含量、配对深度及结果 |
| figures/、figure_index.tsv | 每片 11 张独立图；正式 PNG/PDF，图例可分离 |
| status.tsv、logs/、code/ | 交付核验、分阶段运行日志与 GNU time 资源记录、实际运行代码 |

空间图和浏览器报告均可在本地打开。完整 report 预览使用相对路径，转移时保留 figures/ 目录。原始支持度与条件支持度需要结合 Epi_support 及低分母列阅读。

