# 复核索引

先看[运行总结](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/reports/RUN_SUMMARY.md)和[六方案表](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/05_comparison/six_scheme_comparison.tsv)，再看[优先ROI表](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/ROI_priority_review.tsv)。

完整复核表8257个证据连通ROI；优先队列133个（每切片最多3个）。小岛不删除。优先队列是审阅顺序，不是把其他ROI默认接受。

| roi_id | candidate_label | n_common_anchors | question | impact |
|---|---|---|---|---|
| CanDis_2026_Sample8__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 10161 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| GSE319536_BC11__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 4669 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| CanDis_2026_Sample1__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 4291 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| CanDis_2026_Sample3__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 4269 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| GSE285715_T1__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 3913 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| CanDis_2026_Sample2__tumor_supported_epithelium__0001 | tumor_supported_epithelium | 3871 | 上皮与非整倍体证据相符；形态是否支持肿瘤区域，边缘是否受缺口影响？ | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| GSE285715_T2__uncertain_or_conflicting__0001 | uncertain_or_conflicting | 3496 | 是否正常上皮、肿瘤或其他组织？diploid/缺失CNV不能单独裁定正常。 | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |
| GSE238145_B11__tumor_supported_epithelium__0001 | tumor_supported_epithelium | 3436 | 上皮与非整倍体证据相符；形态是否支持肿瘤区域，边缘是否受缺口影响？ | 改变区域分层/候选界面；不改变共同anchor或六方案标签 |

## 每dataset证据包

- [CanDis__UNVERIFIED](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/CanDis__UNVERIFIED/REVIEW.md)

- [GSE171351](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/GSE171351/REVIEW.md)

- [GSE238145](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/GSE238145/REVIEW.md)

- [GSE285715](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/GSE285715/REVIEW.md)

- [GSE319536](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/GSE319536/REVIEW.md)

- [Inhousedata__UNVERIFIED](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/Inhousedata__UNVERIFIED/REVIEW.md)

- [NPJ__UNVERIFIED](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/01_regions/NPJ__UNVERIFIED/REVIEW.md)

## 关键结果

- [环境组成字典](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/05_comparison/environment_dictionary.tsv)

- [患者与dataset支持](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/05_comparison/cluster_patient_dataset_support.tsv)

- [来源审计图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2/07_source_validation/single_cell_source.png)

大特征、模型、稀疏图留在项目目录；review.zip仅报告、清单、复核表和必要图。
