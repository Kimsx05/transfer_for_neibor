# Module 04 — PASS_WITH_WARNING

完成 7 个冻结 dataset、63 张切片、265,002 个 spot；有效组成 265,002，无效组成 0。
每个 dataset 完成 3 个 k × 15 个 resolution = 45 个候选，共 315 个候选；选出 dataset-scoped niches 共 48 个。

| dataset | n_spots | n_niches | selected_candidate | n_evaluation | evaluation_mode | status |
| --- | --- | --- | --- | --- | --- | --- |
| GSE171351 | 4086 | 5 | k30_r0.4 | 4086 | all_valid_spots | PASS_WITH_WARNING |
| GSE285715 | 17241 | 3 | k30_r0.1 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |
| NPJ__UNVERIFIED | 23800 | 4 | k30_r0.1 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |
| Inhousedata__UNVERIFIED | 31213 | 5 | k30_r0.2 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |
| CanDis__UNVERIFIED | 56357 | 7 | k20_r0.3 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |
| GSE319536 | 58041 | 6 | k10_r0.1 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |
| GSE238145 | 74264 | 18 | k30_r0.8 | 5000 | fixed_section_stratified_sample_approximation | PASS_WITH_WARNING |

输入来自 Module 00 冻结的原始完整 30-factor posterior mean abundance CSV；全 factor 行和归一化为比例。按真实 (section_id, barcode) 对齐，原始 abundance 文件 SHA 与上游冻结 SHA 均未改变。所有有效组成进入发现，无 Epi/CNV 预筛；MP、H&E、空间坐标和患者标签不进入 feature matrix。未做 log、z、ILR、PCA、Harmony 或 c2l 重训。

文献核心来自 [Wang et al. 2025](https://link.springer.com/article/10.1186/s12943-025-02377-9)：cell2location compositions、scran SNN 的 k=10/20/30、Louvain、silhouette 选参和 niche-vs-rest Wilcoxon/FDR。
公开作者实现未恢复；posterior mean proportions、dataset 分池、resolution=0.1–1.5、rank weighting、精确 Euclidean NN、spot-weighted silhouette、固定 section 分层评估集和 BH family 都是透明项目适配，见 methods_parity.tsv。并非逐参数原样复刻。

silhouette 只在不超过 5,000 个固定评估 spot 的原始 30 维组成空间计算。N≤5,000 的 GSE171351 全量评估，其余六组是 sample-based approximation；同 dataset 的所有候选使用完全相同的评估集。预先设定 1e-8 tie tolerance，先低 resolution 再小 k。cluster size、section dominance、centroid redundancy 仅记录诊断，未据此改解或粗化。

CPU 预算8：dataset 逐组串行，组内4个 candidate workers，graph/BLAS/OpenMP均1线程，无嵌套并行或GPU。前两组的单 worker 和4-worker全部候选 memberships、silhouette逐值一致，见 execution_parallelism_parity.json。内存为逐进程 VmHWM；父进程+4倍最大worker峰值是保守上界，会重复计入fork共享内存，不宣称是真实合计RSS。

Wilcoxon 为独立两组的双侧 rank-sum，使用 ties 与 continuity correction；BH 在各 dataset 的所有 niche×factor tests 内校正。输出 U、rank-biserial effect、组内/组外 median fraction、差值、mean fraction、绝对 mean abundance 和 FDR。优化实现与 scipy 独立算法核对通过。注释描述组成特征，不是独立验证，也不能将空间 spot 当独立生物重复。

niche 描述使用实际的前三个 mean fractions；ID 仅在 dataset 内定义，跨 dataset 同编号没有对应意义。未自动合并标签；major 聚合只用于解释图。低可辨识性原始 subtype 标记随表提供，所有 T/NK factors 保留。

空间图 63 张，组成/major/选参/coverage 图 28 张，共 91 views，均 PDF + 600 dpi PNG，另有缩略图。优先复用原 R save_both 函数，调整 PNG 至600 dpi；H&E alpha=0.25，spot alpha=1。零值黑色，NA 灰色，使用冻结 lowres 坐标，无再次缩放。视觉抽查覆盖七个 dataset。

警告表示完成的分析存在解释或复刻边界，不表示运行报错：
- Posterior mean proportions, dataset pools, resolution grid, SNN rank weighting, evaluation sampling and detailed selection rules are transparent project adaptations; author code not recovered.
- Six datasets use a fixed 5000-spot silhouette approximation; GSE171351 uses all valid spots.
- Frozen patient_uid absent; patient coverage is SKIPPED_NOT_ESTIMABLE; section coverage is available.
- Low-identifiability subtype factors retained and flagged; defining Wilcoxon/BH annotations are not independent validation.
- CanDis, Inhousedata and NPJ dataset identities retain frozen UNVERIFIED labels.
- Installed scran 1.38.1 differs from paper 1.32.0; actual interfaces and versions recorded.

patient_uid 全缺失，因此患者层覆盖标 SKIPPED_NOT_ESTIMABLE；切片层覆盖正常输出。未推断患者编号。

交付：input_composition_manifest.tsv、selected_parameters.tsv、silhouette_grid.tsv、spot_niche_assignments.parquet、niche_centroids.tsv、niche_annotation_wilcoxon.tsv、niche_sample_coverage.tsv、methods_parity.tsv；另有图、代码、版本、校验、逐 dataset 日志、output_manifest.tsv、handoff.json 和 review_bundle.zip。

**停止位置：Module 04 完成。未运行 MP–niche 检验、CNV 关联或下一模块。原始数据、历史 MP 图与 Module 00/01 输出均未覆盖。**
