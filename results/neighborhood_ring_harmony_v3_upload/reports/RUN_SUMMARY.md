# Niche v3 实际完成情况

18/18 主配置已完成：NONE、DATASET、SAMPLE × R1、R2、R12 × KM、Leiden。另保留全部 R12 五表示 pilot。共同 anchors=129632，完整双环=108345，reference=20000；63 张目标切片中 61 张有 anchors。

## 主要结果

|配置|群数|稳定性 ARI 中位数|标签接受率|状态|

|---|---:|---:|---:|---|

|DATASET_R12_KM|6|0.988|91.3%|niche_working_candidate|

|DATASET_R12_LEIDEN|8|0.893|97.3%|niche_working_candidate|

|DATASET_R1_KM|6|0.993|90.0%|niche_working_candidate|

|DATASET_R1_LEIDEN|8|0.966|97.8%|niche_working_candidate|

|DATASET_R2_KM|6|0.983|91.1%|niche_working_candidate|

|DATASET_R2_LEIDEN|8|0.817|97.7%|niche_working_candidate|

|NONE_R12_KM|6|0.970|90.7%|niche_working_candidate|

|NONE_R12_LEIDEN|27|0.813|97.0%|niche_working_candidate|

|NONE_R1_KM|6|0.973|90.2%|niche_working_candidate|

|NONE_R1_LEIDEN|6|0.999|99.0%|niche_working_candidate|

|NONE_R2_KM|6|0.966|90.9%|niche_working_candidate|

|NONE_R2_LEIDEN|25|0.723|96.8%|exploratory_with_limitations|

|SAMPLE_R12_KM|6|0.976|92.0%|niche_working_candidate|

|SAMPLE_R12_LEIDEN|9|0.971|97.3%|niche_working_candidate|

|SAMPLE_R1_KM|6|0.988|92.2%|niche_working_candidate|

|SAMPLE_R1_LEIDEN|8|0.838|97.3%|exploratory_with_limitations|

|SAMPLE_R2_KM|6|0.989|93.8%|niche_working_candidate|

|SAMPLE_R2_LEIDEN|5|0.954|97.9%|niche_working_candidate|

## 三种整合模式的收益与代价

DATASET: selected WEAK; strong mixing improvement=0.04066, raw Hellinger cost=8.269%, stability preserved=False, coverage preserved=True, composition preserved=False. Conservative WEAK when tradeoffs conflict. No MP10 or UMAP used.

SAMPLE: selected WEAK; strong mixing improvement=0.12318, raw Hellinger cost=28.968%, stability preserved=False, coverage preserved=False, composition preserved=False. Conservative WEAK when tradeoffs conflict. No MP10 or UMAP used.

Final finite-grid candidate reconciliation checks positive original-composition separation where supported. This changes some Leiden working granularities but does not change either frozen Harmony strength; the raw neighbor-distance cost alone exceeds the strong-selection bound in both modes. Initial choices remain archived.


|R12 表示|同 dataset 邻居|同 sample 邻居|原始 Hellinger 邻居距离|

|---|---:|---:|---:|

|DATASET_STRONG_R12|61.0%|59.7%|0.2070|

|DATASET_WEAK_R12|65.1%|63.6%|0.1912|

|NONE_R12|94.1%|90.3%|0.1422|

|SAMPLE_STRONG_R12|39.5%|26.8%|0.3014|

|SAMPLE_WEAK_R12|51.4%|39.1%|0.2337|

来源混合不能单独证明批次效应被正确消除；patient/clinical 差异也可能是真实组成差异。请结合原始组成距离、来源支持、稳定性和拒绝率评估，暂不宣称唯一最佳方案。

## 工作粒度和组成解释

v2 的 K2 是粗粒度结构，原 Leiden 38–45 群偏细。本轮 K2/4 仅作粗诊断，工作 K-means 从 K6 起评估；Leiden 从更低分辨率评估。负轮廓系数或强来源依赖的配置被保留为探索结果，不命名为已验证的病理 niche。详细逐配置选择依据见比较目录。

[niche_dictionary.tsv](../04_comparison/niche_dictionary.tsv) 给出每个 niche 的 R1/R2 相对富集、来源共享、真实患者证据与来源主导标记。SPP1 TAM、CAF、T/NK 等只按原始组成描述，不推断纯细胞群或机制。

## 结果入口

- [18方案比较](../04_comparison/eighteen_scheme_comparison.tsv)
- [图形目录](../04_comparison/figures/)
- [逐数据集空间图](../04_comparison/by_dataset/)
- [ROI 操作说明](../05_user_review/WHAT_USER_NEEDS_TO_DO.md)
- [轻量回传包](review.zip)

## 你现在需要做什么

先查看 ROI 操作说明中的五个重点区域，填写 ROI_review_to_fill.tsv 的 reviewed_label、review_comment、reviewer、review_date。允许 uncertain 和 mixed_needs_split；混合区域可返回 polygon/spot成员表，不强迫整个区域单一标签。人工判断不影响本轮 niche 完成。

## 停止范围与续跑

本轮止于 niche。MP10 分析、解盲、来源验证全部 DEFERRED_BY_USER，留待后续授权。未运行旧 v2 resume.sh。

```bash
bash "/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/code/resume_niche_only.sh"
```

续跑命令复用已完成 checkpoint，仅生成/补全本轮 niche 产物，不触发 MP10。

## 建议优先复核的工作结果

K-means 可先以 DATASET_WEAK 的 R12 K6 作为较保守的整合工作视图，同时保留 NONE K6 作为未校正基线；SAMPLE_WEAK 的 R12 Leiden K9 可作为细分环境的补充视图。它们是不同用途的工作候选，不宣称唯一最佳；具体稳定性、覆盖率和组成代价见上表。

NONE R2/R12 Leiden 的低分辨率中等群数解在原始组成上分离不足，因此工作结果采用有限网格中的25/27群，R1保留有正轮廓证据的6群。它们仍有来源集中限制；v2的38–45群不被直接当作生物学真值。

|SAMPLE_WEAK R12 Leiden ID|占共同anchors|R1 相对富集|R2 相对富集|支持数据集|最大单dataset占比|

|---|---:|---|---|---:|---:|

|0|34.2%|Epithelial 82.1% (+0.34 log2)|Epithelial 80.5% (+0.34 log2)|6|31.3%|

|1|21.9%|GC B 15.0% (+1.32 log2), Naive/Memory B 3.7% (+0.46 log2), Epithelial 59.1% (-0.13 log2)|GC B 14.7% (+1.30 log2), Naive/Memory B 3.8% (+0.41 log2), Epithelial 58.7% (-0.11 log2)|7|35.6%|

|2|15.6%|Plasma 4.2% (+1.32 log2), Naive/Memory B 6.1% (+1.20 log2), Epithelial 51.3% (-0.33 log2)|Plasma 4.6% (+1.31 log2), Naive/Memory B 6.4% (+1.18 log2), GC B 4.8% (-0.33 log2)|7|27.5%|

|3|11.6%|Matrix_CAF 7.7% (+1.78 log2), Lymphatic EC 3.1% (+0.94 log2), vCAF 3.6% (+0.91 log2)|Matrix_CAF 7.9% (+1.71 log2), Lymphatic EC 3.0% (+0.94 log2), vCAF 3.6% (+0.88 log2)|7|35.2%|

|4|3.9%|vCAF 8.5% (+2.14 log2), Lymphatic EC 5.5% (+1.77 log2), Epithelial 59.2% (-0.13 log2)|vCAF 8.3% (+2.08 log2), Lymphatic EC 5.4% (+1.77 log2), Epithelial 59.0% (-0.11 log2)|7|53.0%|

|5|3.5%|SMC 11.0% (+3.63 log2), Matrix_CAF 4.4% (+0.97 log2), vCAF 3.2% (+0.71 log2)|SMC 11.9% (+3.54 log2), Matrix_CAF 4.5% (+0.92 log2), vCAF 3.2% (+0.71 log2)|7|27.7%|

|6|4.4%|SPP1 TAM 8.0% (+2.55 log2), Lymphatic EC 3.3% (+1.06 log2), Neutrophil 4.0% (+0.99 log2)|SPP1 TAM 7.9% (+2.50 log2), Lymphatic EC 3.3% (+1.06 log2), Neutrophil 4.2% (+1.01 log2)|7|50.9%|

|7|3.1%|Neutrophil 15.4% (+2.94 log2), SPP1 TAM 3.6% (+1.41 log2), Epithelial 56.1% (-0.21 log2)|Neutrophil 16.5% (+3.00 log2), SPP1 TAM 3.7% (+1.41 log2), Epithelial 53.4% (-0.25 log2)|6|46.8%|

|8|1.8%|Capillary EC 7.9% (+2.85 log2), Pericyte 3.6% (+2.31 log2), Naive/Memory B 4.4% (+0.73 log2)|Capillary EC 6.6% (+2.50 log2), Naive/Memory B 4.5% (+0.65 log2), GC B 6.3% (+0.07 log2)|6|32.4%|

这些是原始 q05 双环组成的相对富集描述。多数群仍有大量 Epithelial 成分；不能把8% SPP1 TAM 或8% Matrix_CAF称为纯髓系/纯CAF niche。R1/R2总体相近，部分群如中性粒细胞富集群在R2更高，尚不能据此声称组织界面或作用方向。

来源主导群详见 niche_dictionary.tsv；明确共享需结合多来源支持而非要求每群包含全部7个dataset。已知患者数与匿名/假设分析单位分别记录。

SAMPLE R1 Leiden 的工作候选在原始 Hellinger 空间仍为负轮廓系数，属于探索性结果；高条件稳定性并不能消除这一限制。DATASET R1 最终取8群，SAMPLE R2取5群，均是既定有限网格中的组成分离受支持候选。
