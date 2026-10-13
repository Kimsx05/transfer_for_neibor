# MP10 v5：启动与审阅说明

## 本轮做什么

固定现有 niche0–niche8，回答：Top20是否能代表原MP10，以及cell2location上皮归属是否能保留真实状态差异。

这轮同时完成：完整上皮单细胞评分保真度、已保存模拟的真实上皮恢复、有限合成模型校准、三个同样本直接对比（niche1−2、niche1−6、niche7−6），以及全部九群/七dataset的证据汇总。

不重跑niche；不以预期正向更换signature；不启动通讯、预后或新机制实验。病理未知不阻塞本轮，取消的代表切片人工复核不再加入。

## 发给 Codex

把 `neighborhood_mp10_calibration_codex_prompt_v5.md` 发给服务器 Codex，然后发送：

```text
请完整读取并执行附件 neighborhood_mp10_calibration_codex_prompt_v5.md。

以 neighborhood_mp10_source_v4 为只读来源，输出到新的
/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_mp10_calibration_v5。

本轮重点是实际比较模型估计上皮与模拟中已知真实上皮，
不要只重复counts重构验证，或只报告原模型权重回放就称端到端验证。

同时完成 niche1−2、niche1−6、niche7−6 的同切片/分析单位直接比较，
按dataset汇总并保留全部九群结果。不能用两个vs-rest效应相减代替。

允许按附件上限运行最多两个小规模合成数据c2l拟合；
不重训真实切片模型，不修改niche、Top20或已有分数。

不要只做预检，不等待全量病理审定。完成校准、真实对比、证据表
和轻量回传包后停止，返回实际完成情况及可执行续跑命令。
```

## 首轮应该拿到什么

先看 `reports/RUN_SUMMARY.md` 与 `completion_matrix.tsv`，确认恢复验证究竟做到了哪一层：

- `FROZEN_WEIGHT_REPLAY`：只回放既有真实模型权重；是压力测试。
- `KNOWN_MIXTURE_FIXED_SIGNATURE`：已知混合组成下检查固定reference状态分配。
- `SYNTHETIC_REFIT`：由合成观测counts重新估计丰度，再进行归属；仍是技术benchmark。

重点审阅：

1. `03_attribution_calibration/calibration_report.md`、`paired_state_recovery.tsv`：真值能否恢复，是否出现背景相关反转。
2. `01_score_fidelity/fidelity_by_patient_dataset.tsv`：是否只在合并数据相关，患者内是否成立。
3. `04_fixed_niche_contrasts/paired_dataset_effects.tsv`：三个直接对比是否存在共同样本支持。
4. `05_evidence_synthesis/evidence_matrix_all9.tsv`、`priority_decision.md`：技术证据和真实关联是否支持同一个方向。

## 什么时候停下来审阅

完成上述结果和报告后停止，不自动进入CellChat/NicheNet、预后、药物或机制分析。

若来源校准显示系统偏差，先报告偏差发生的背景/深度/基因和失败范围，不能为了恢复期望正向而改模型。

若模拟尚未收敛或关键reference/scRNA输入不完整，报告技术未解决，继续其余已授权分析；不给出生物阴性或方法已无效的结论。

负向、不显著、队列异质或来源仍不清楚，都是允许返回的完整科学结果。
