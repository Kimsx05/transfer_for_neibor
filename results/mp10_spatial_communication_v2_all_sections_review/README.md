# MP10 spatial communication v2 — all 63 ordinary Visium sections

**全量任务已终态；成功任务均已完成表图**。生成时间 2026-10-29T06:37:10+08:00。请首先阅读 [中文报告](中文报告.md)。状态基于明确的任务/后处理完成标记；运行中不是全量完成。

- [逐片输入、身份、组织、评分及状态](comparison/section_input_identity_coverage.tsv)
- [126任务运行与缓存位置](comparison/section_method_status.tsv)；实时权威表 [logs/task_status.tsv](logs/task_status.tsv)
- [全部逐片H/L、flow、双份额和支持](comparison/all_section_contrasts.tsv)
- [患者与未核验记录](comparison/patient_and_unverified_record_summary.tsv)、[dataset](comparison/dataset_stratified_summary.tsv)、[dataset等权](comparison/dataset_equal_summary.tsv)、[LODO](comparison/leave_one_dataset_out.tsv)
- [方法一致性](comparison/method_concordance_summary.tsv)、[共同LR敏感性](comparison/common_LR_sensitivity.tsv)、[质量/几何诊断](comparison/batch_geometry_diagnostics.tsv)
- [全图索引](figure_index.tsv)、[方法](methods.md)、[实际参数](parameters.tsv)、[脚本](scripts/)、[轻量review包](review.zip)

两种方法单位不同；先片内H/L，再同组织/治疗患者内汇总，再dataset。unknown不当tumor、未核验记录不当核验患者，19/20单列。0、NE与failed/blocked严格区分。所有谱系相同规则，候选判断不继承试点排序。

新目录只保存v2结果；v1缓存只读以source_dir引用。恢复入口 `bash scripts/resume.sh`；全量自动最终链 common_lr.py → summarize.py → report.py。输出中有queued/running时不得称为完成。
