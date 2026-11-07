# 诊断字段与复跑说明

`<section>/diagnosis.tsv.gz` 是唯一逐对象诊断表；键为 `section + cell_id`。固定原分析集合，不含本来未进入 qc_keep_analysis 的对象。

- `final_major`, `annotation_status`, `annotation_reason`, `raw_prediction`, `reference_score_max/margin`：原流程结果，未修订。原九类详细参考支持仍位于 `../hd_major_annotation_v1/<section>/tables/annotations.tsv.gz`。
- `diag_original_reason`：原始历史拒绝/接纳理由；`diag_first_rejection_reconstructed`：按原代码书写的主规则顺序重建，仅用于解释最终Unknown。不是历史短路日志。
- `diag_fail_*`, `diag_rule_overlap`, `diag_n_main_rule_failures`：可重叠。`reference_band` 主要影响置信状态及单marker补充，不是所有对象的硬拒绝规则。
- `diag_all_marker_strength_zero`, `diag_marker_top_tied`：原规则argmax证据缺失或并列诊断；此时argmax的首个类别不能视作成组谱系证据。
- `singler_{blca,hpca}_labels`：fine-tuning后原始候选；`pruned_labels`：全片统一默认pruning后候选，NA包含被prune或不可评估者，具体看status/pruned。
- `singler_*_status`：evaluated 或 not_assessable_constant_expression_on_model_genes。后者仍保留完整ID及原始数据。
- `singler_*_delta_next`：SingleR原生fine-tuning诊断；`delta_median`：最终所选类别初始分数相对各类别分数中位数的差。
- `singler_*_initial_top_label`, `selected_initial_score`：fine-tuning前分数诊断；最高分列可能不同于最终labels。分数不是概率，不与Seurat支持分数直接数值比较。
- `diag_candidate_pattern`：两参考候选一致、不同、单路或均无；纯描述，不是新的接纳门槛。
- `diag_*_candidate_marker_count`, `diag_*_marker_evidence`：原完整RNA panel的检出/强度描述。`other_group_supported`要求另一组≥2枚marker且top3平均log强度不低于候选；这不是自动拒绝规则，也不是已证实混合。`no_candidate_marker`是缺支持，不是反证。
- `diag_hpca_*_coverage`：HPCA六类覆盖限制；Mast/Plasma/Mural未覆盖。原Unknown的真身份不已知，因此基于BLCA候选的覆盖描述也是条件性的。

`singler_*_scores.rds` 是 `list(scores=matrix, meaning=...)`，行名cell_id、列名训练大类，包含全部固定分母。不可评估行NA。不是表达矩阵。

`candidate_marker_statistics.tsv` 在完整分析集合上按原Unknown/已知及SingleR候选分组，报告原RNA log-normalized平均值与原counts>0比例；没有插补、缩放或群体差异检验。

`rerun.sh` 使用现有服务器环境和缓存。模型与预测存在时复用；报告和review.zip最后更新。要有意重新绘图，删除对应片的 `roi_comparison.tsv` 完成标记后复跑；要更改核心算法或输入，须先另行确认重建相应缓存，不可直接绕过cache_manifest.json。该manifest只对相关输入做大小/mtime记录、对少量核心脚本做签名，不哈希大数据。

HPCA来自隔离安装的celldex 1.12.0；其官方源包存档版本用于复用已安装ExperimentHub依赖，HPCA固定1.0.0参考来自官方资源，不使用其他标签库代替。原分析R库没有被替换。训练缓存和ExperimentHub文件留服务器，不进入review.zip。
