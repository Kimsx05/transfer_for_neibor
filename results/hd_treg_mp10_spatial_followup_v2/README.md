# hd_treg_mp10_spatial_followup_v2

先读 [REPORT.md](REPORT.md)。三分支均已执行；未采用候选新证据门槛，五片v1主结果保留。

```r
source("/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_followup_v2/code/load_result.R")
r <- load_followup("/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_followup_v2", "Cervilla", "regions")
p <- load_followup("/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_followup_v2", "Cervilla", "paired")
```

```python
import sys
sys.path.insert(0, "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_followup_v2/code")
from load_result import load_followup, attach_region_evidence
regions = load_followup("MPBC2")
members = load_followup("Cervilla", "membership")
```

区域分数不复制为独立bin观测。加载器仅给原T32区域表追加fu2_证据字段，或返回bin→T32成员ID；不会覆盖旧labels/counts/QC。通用AnnData/Seurat对象仍通过上游loader读取。

- `tables/Cervilla_paired_analysis.parquet` + `design_r200_common.parquet`：共同247中心、两来源结局与同一数值X。
- `paired_effects.tsv`、`new_test_registry.tsv`：主/辅助检验来源与家族。
- `common_support_*`：有限counts重建、整数总量与成员，非全片新细胞对象。
- `expanded_mixtures.parquet`/`mixture_rule_results.parquet`、`mixture_parent_membership.parquet`：混合设计、真值、单规则/联合结果与参考成员。
- `*_evidence.parquet`：五片全部T32旧/新证据及秩窗口解释。
- `ROI_*`、`ROI_REVIEW.md`：24区域原始支持、完整50/200邻居及图册。
- `v1_five_section_main.tsv`：仅转录v1全五片主结果，不是本轮重新估计。

依赖根为同级 `hd_treg_mp10_spatial_all_v1`、`hd_mp10_treg_endo_precommunication_v1`、`segmentation_assignment_v2` 与manifest的真实counts/图像/参考。移动根目录时保持analysis同级结构；原始绝对路径可通过manifest检查，包内没有假装自包含的表达矩阵。

执行 `bash rerun.sh` 续跑；`bash rerun.sh --force` 重算本轮各步；`bash replot.sh` 只重绘。默认沿用项目23:00–07:00窗口；明确授权即时计算时使用 `HD_ALLOW_DAYTIME=1`。配置变化时不静默使用旧结果，应另建版本并检查输入。review.zip不包含counts/参考/大型对象。
