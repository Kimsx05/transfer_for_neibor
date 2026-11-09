# 固定 8×8 μm bin 基线与细胞对象比较

两片全量基线、两路 SingleR、固定 ROI 原始信号对照、Cervilla 作者方法子集对照和两片 ROI RCTD 均已完成。原有细胞正式标签、分割、归属、QC 和表达对象未覆盖。8×8 μm 方格与“核边界外扩最多 8 μm”是不同单位。

## 完成状态与输入

| section | before_bins | before_UMI | QC_bins | QC_UMI | genes |
| --- | --- | --- | --- | --- | --- |
| Cervilla | 297675 | 64528669 | 268924 | 64443819 | 18085 |
| D1_Johnson | 400147 | 131539416 | 364805 | 131440767 | 18085 |

两片均使用原生 square_008um。Cervilla 原生文件从已定位的作者 Zenodo 原始 ZIP 中按成员范围恢复，仅下载矩阵与坐标；不是由细胞矩阵反推。新恢复的原始文件保存在 Cervilla/native8。QC 为 in_tissue 且 nCount > 5；Cervilla 阈值来自公开作者代码，D1 为本项目共同对照设置。全部分析 bins 均保留 18,085 个输入基因，不再以 nFeature 过滤。原生输入全量 QC/覆盖字段在每片 bin_metadata.tsv.gz。

本轮按片独立 LogNormalize(10,000)、3,000 HVG、50,000 LeverageScore sketch、50 PCs / 前 30 维、Louvain resolution 0.8、UMAP，再投影到全体 QC bins。seed=123。SingleR 直接预测所有可评价 bins；不是 sketch 标签传播。两片主 RDS 已实际重读，counts、ID、全量 PCA/UMAP、FOV 与物理坐标核对通过，并由重读对象生成 readback_UMI 图。

方法来源、项目控制参数与不能确认的原文细节见 [METHODS.md](METHODS.md)。软件版本见每片 sessionInfo.txt；SingleR 2.12.0，Seurat 5.5.1，SeuratObject 5.4.0；原版 spacexr 2.2.1。

## 共同覆盖、表达回收与边界差异

| section | common_supplied2_UMI | previous_unassigned_in_QC8_UMI | previous_unassigned_fraction | analysis_cell_UMI_in_QC8 | native8_minus_common2_in_QC |
| --- | --- | --- | --- | --- | --- |
| Cervilla | 64348326 | 3707867 | 0.0576 | 60626009 | 95493 |
| D1_Johnson | 131364518 | 3669471 | 0.0279 | 127664092 | 76249 |

表中分母是 QC 8 μm 方格所覆盖的 supplied 2 μm 原始 UMI。previous_unassigned 由原始 2 μm bin→cell 映射直接统计，不通过最近质心配对。native8_minus_common2 是原生 filtered 输入/边界掩膜差异，不能算作细胞分配改善。未进入细胞分析 QC 的已分配表达另有贡献；完整账目保存在 bin_metadata.tsv.gz 和 gene_count_accounting.tsv。原生 8 μm QC 前比 supplied 2 μm 多出的 UMI 分别为 96,873 和 78,796。所有 supplied 2 μm 都有对应的原生 8 μm 网格父格；对这些相同输入，逐基因聚合守恒。

Cervilla 4,937、D1 5,045 个原生 bins 仅有部分 supplied 2 μm 覆盖；coverage_fraction_2um 记录实际提供的小格数/16，未知小格不当作已测零。Cervilla 图像中无 supplied bins 的区域继续用旧覆盖遮罩区分，不产生“Unknown 细胞”。图像与网格配准沿用验证变换；D1 裁剪原点 (23764,35851) 保留。原生 bin 中心与 2 μm 仿射外推的最大差异均 <0.05 μm，最终使用原生空间坐标。

## SingleR 与信息量对照

复用冻结 BLCA 44,896 个参考细胞的九大类模型，以及 HPCA 十个原生 main 类合并成六大类的模型。BLCA 包含 Mural；HPCA 主对照不能公平区分 Mast、Plasma、Mural。基因/标签映射和训练对象路径见 METHODS.md。按每片、每参考合并所有 5,000-bin 分块后，仅统一 pruning 一次。scores 是相关性；delta.next 为官方输出，不是自行计算的普通 top1−top2。

| section | unit | n | UMI | median_nFeature | covered_disagree_fraction | epi_fib_ge2_fraction | reference_disagree_UMI_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | cell8 | 93593 | 60670091 | 300 | 0.0980 | 0.2867 | 0.0775 |
| Cervilla | bin8 | 268924 | 64443819 | 114 | 0.1365 | 0.1436 | 0.1072 |
| D1_Johnson | cell8 | 188237 | 127719635 | 247 | 0.3115 | 0.3368 | 0.3249 |
| D1_Johnson | bin8 | 364805 | 131440767 | 179 | 0.3204 | 0.3141 | 0.2921 |

covered_disagree_fraction 分母为两路均有候选且 BLCA 类型在 HPCA 六类覆盖内的对象；epi_fib_ge2_fraction 为两组各检出至少两枚既有 marker 的对象比例。最后一列为两参考分歧对象所承载的 UMI / 该单位全部分析 UMI。细胞与 bins 的单位、QC 和分母不同，这些数值不能解释为准确率变化。

### Cervilla

BLCA 可评价 266,501/268,924，默认 pruning 拒判 33；HPCA 可评价 261,586/268,924，拒判 0。其余保留 ID，并标为模型选定基因表达恒定、不可评价。未被 pruning 不等于高置信度。

nFeature≥500 层：细胞→bins 的上皮/胶原双组共检出为 48.06%→30.07%；覆盖类别内参考分歧为 3.65%→1.83%。该层内部的 nFeature 分布仍不相同，因此也不是严格匹配检测量的纯度实验。完整分层见 information_strata.tsv。

### D1_Johnson

BLCA 可评价 360,348/364,805，默认 pruning 拒判 3；HPCA 可评价 352,131/364,805，拒判 0。其余保留 ID，并标为模型选定基因表达恒定、不可评价。未被 pruning 不等于高置信度。

nFeature≥500 层：细胞→bins 的上皮/胶原双组共检出为 71.56%→68.80%；覆盖类别内参考分歧为 27.09%→28.98%。该层内部的 nFeature 分布仍不相同，因此也不是严格匹配检测量的纯度实验。完整分层见 information_strata.tsv。

## 固定 ROI 与 RCTD

| section | ROI_bins | ROI_UMI | fitted_bins | fitted_UMI | fitted_UMI_fraction | below_UMI_min100 | doublet_certain | doublet_uncertain | not_in_bin_QC | reject | singlet |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 4569 | 1059406 | 2271 | 981476 | 0.9264 | 2262 | 382 | 228 | 36 | 99 | 1562 |
| D1_Johnson | 4460 | 786929 | 2532 | 730470 | 0.9283 | 1718 | 376 | 275 | 210 | 291 | 1590 |

同片五个旧 ROI 合并拟合；输入为未转换 counts，BLCA 参考也为原始 counts。Reference n_max_cells=50,000，实际保留同一批 44,896 细胞。doublet 模式只拟合一至两种主要参考成分，不代表真实 bin 最多含两类细胞。默认 UMI_min=100、UMI_min_sigma=300、counts_MIN=10；其他参数逐片写入 rctd_config.json。ROI 外为 not_tested，无投影。

结果表保留 first_type、second_type、weight_1/2、spot_class、拟合分数及收敛字段。weight_1/2 是候选两成分拟合权重，singlet 的第二权重也可能非零；不能脱离 spot_class 将其解释成已确认 doublet。reject 的候选与权重仍保留供追溯，不当作通过拟合的身份。低于算法门槛/未过 bin QC 的全部 ROI bins 也在 rctd_ROI_results.tsv.gz。

- Cervilla：2,271/4,569 个原生 ROI bins 进入拟合，覆盖 ROI UMI 的 92.64%。far_UMI_hotspot 的上皮—成纤维参考分歧 bins 中，81/81 个拟合；44 个为这两类的 doublet_certain，另 14 个为 doublet_uncertain。
- D1_Johnson：2,532/4,460 个原生 ROI bins 进入拟合，覆盖 ROI UMI 的 92.83%。stroma 的上皮—成纤维参考分歧 bins 中，67/120 个拟合；20 个为这两类的 doublet_certain，另 7 个为 doublet_uncertain。

D1 stroma：实际查看 H&E、原始 2 μm、4/8 μm 细胞归属及固定 bins。上皮/胶原信号在原始尺度已交错，固定 bins 仍有共检出；RCTD 在部分分歧 bins 中支持两成分解释。BLCA/HPCA 差异在同一固定 bin 上存在，因而不能全归因于扩张。仍不能凭这些图分清相邻细胞、胞质延伸、局部混合与背景信号各占多少。dense/small_nuclei 作为免疫相关区域对照，far_UMI_hotspot 两参考主要支持上皮信号，但部分强信号位于 H&E 腔隙/缺乏核支持处，不能作为明确细胞身份对照；该 ROI 包含 95,896 个此前未归属的 supplied 2 μm UMI（28.53%）。这是覆盖收益，也提示胞外/背景信号或局部配准支持不足的可能，未在本轮强行归入细胞；edge 极低信息，D1 edge 无 bins 达到 RCTD 门槛，不能对它给出解卷积改善结论。

同一批对象的 4→8 μm 归属对照：Cervilla 热点的双组各≥2 marker 对象为 202/236→224/236；D1 stroma 为 97/302→118/302。扩大归属确实增加共同检出，但 4 μm 版本已存在明显信号竞争，不能将问题全部归咎于 8 μm 扩张。详见 assignment4_vs8_ROI.tsv。Cervilla far_UMI_hotspot：固定 bins 仍有明显上皮/胶原共同信号；dense、stroma、small_nuclei、edge 的结果分别保留，不把“细胞对象更多表达”视为已分离。免疫/内皮原始 panel 信号另见 ROI_immune_endothelial_small_nuclei。marker 采用既有面板，未包含基因与包含但零检出在 marker_coverage.tsv 区分；KRT18、TPSAB1、TPSB2 未在输入面板中，不将其缺失解释为阴性生物学。

## Cervilla 作者方法的有限敏感性分析

固定区域/区域外 × nFeature 分层抽取 4,968 个相同 bins；作者 HPCA 保留 T/B fine 标签再映射到大类，独立于主对照的六大类 HPCA 模型。B_cell:Plasma_cell 明确映射为 Plasma，其余 B fine 标签为 B；这少量 fine 参考不等于主 HPCA 覆盖全部 Plasma。完整映射见 author_HPCA_label_mapping.tsv。

sctransform 0.4.3，v2、nb_offset（未安装 glmGamPoi，使用原生拟合），ncells=4968、conserve.memory=TRUE。输入 SingleR 为 SCT data，即 log1p 校正 counts，不是 scale.data/Pearson residuals。共用 11,144 个基因训练作者标签模型，两路分别预测。SCT 留下的表达基因比临时 symbol counts 少 6,204 个；这是该子集建模保留范围，未修改主对象的完整基因。

两路均可评价 4,924 个；其中 297 个大类发生变化（6.03%）。作者 fine-label 变化数、分层结果与 UMI 见 author_information_strata.tsv、author_major_crosstab.tsv。部分 LogNormalize 不可评价的低信息 bins 在 SCT 校正后得到候选，不能将此称为新增实测信号。

该子集有 2,113 个 bins 检出≥2 枚上皮 marker、1,443 个检出≥2 枚成纤维 marker、155 个检出≥2 枚髓系 marker；集合可重叠，不是细胞类型计数。T/NK ≥2 marker 仅 5 个，B 为 0；稀有 T/B 对照代表性不足。SCT 丢失的输入中可用 marker 包括 KRT6A, KRT6B, KRT4, KRTDAP, DSG3, CD3D, CD3G, CD247, LCK, NKG7, KLRD1, GNLY, PRF1, EMCN, MS4A1, CD79A, CD79B, CD22, CD19, TNFRSF17。两路主要比较使用相同的交集基因，仍不能推广到全片免疫细分。原生 NB 拟合出现 iteration limit / alternation limit 警告，已保留在 author_warnings.txt；本次结果仅作有这一限制的敏感性对照。没有据此改写主线或宣称精确复现原论文。归一化会影响局部候选，尚不足以改变“换 bin 不能整体消除混合与参考差异”的判断；不自动扩展到全片 SCT。

## 后续单位建议

| 切片 | 建议 | 证据与限制 |
|---|---|---|
| Cervilla | 保留细胞主线，以固定 bins 补充空间表达覆盖与争议区域复核 | 细胞信息量较高、整体参考一致性较好；bins 回收未归属信号，并在部分高信息区域减少共同检出，但增加低信息单位，热点仍混合。类型/空间位置使用有表达支持的可靠细胞；连续组织表达和漏分配热点同时查看 bins。 |
| D1_Johnson | 固定 bins 作为主要组织空间表达基线，与可靠细胞群配合 | bins 提供更完整、独立于核外扩张的信号覆盖；高信息上皮—胶原共检出和参考分歧仍明显。适合研究组织表达分布，不适合将每个 bin 强制视作单细胞。保留明确免疫/上皮等细胞证据作定位辅助。 |

现有证据不能证明换 bin 提高注释准确率。D1 更支持混合表达与参考差异共同存在，局部归属影响尚不能排除；如下一轮需要改细胞标签，应优先限定在 stroma/热点中有跨边界证据的对象，逐一结合原 2 μm 归属检查，不重跑全片核分割。没有启动 16 μm 全片/ROI 扩展：两个 8 μm 基线整体具备可用信息，主要争议并非普遍检测量不足。

## 交付与读取

每片只保存一份 seurat.rds：原始稀疏 RNA counts、LogNormalize data、全输入基因、HVG、sketch 与 full reductions、空间中心、QC、直接 SingleR 候选和 marker 字段。移除不必要的 sketch scale.data；没有重复 h5ad/全量 SCT 对象。原始输入/QC 未通过 bins 在 bin_metadata.tsv.gz，矩阵仍从原生输入路径读取。大 H&E、核 mask、细胞归属结果通过原目录引用。

```r
library(data.table)  # project installation: load before Seurat
library(Seurat)
root <- '/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_bin8_comparison_v1'
source(file.path(root, 'code/load_result.R'))
cervilla <- load_bin8(file.path(root, 'Cervilla'))
# d1 <- load_bin8(file.path(root, 'D1_Johnson'))
counts <- LayerData(cervilla, assay='RNA', layer='counts')
umap <- Embeddings(cervilla, 'umap.full')
# singler_* fields are candidates; rctd_status outside fixed ROIs = not_tested
```

load_bin8 默认按 ID 附加 ROI RCTD 表，核/图像依赖按结果目录相对路径解析。坐标 metadata 为原 fullres pixels 与物理 μm；FOV/绘图使用负 y。每片 gene_mapping.tsv、transform.json、fixed_rois.tsv、bin_metadata.tsv.gz、predictions.tsv.gz、两份 scores.rds 和 ROI 结果可以独立复核。summary tables 位于根 tables/ 与每片目录。

复跑：`bash rerun.sh`。已完成的输入、结构、预测和拟合缓存按阶段复用；固定配置改变时请使用独立的新版本输出目录。入口最后重新汇总报告与 review.zip。使用最多 8 个主数值线程、BLAS=1、RCTD 4 workers，本次并发总线程控制在 48 以下，峰值实际内存远低于 200 GB；没有替换原 R 环境。spacexr、mapproj、pals 的分析加载路径在本轮 Rlib。依赖安装曾意外新增用户库 pals/mapproj；自动审批拒绝迁移已有库目录，随后复制 pals 到隔离库并显式加载，未删除/替换用户库内容。

review.zip 仅含报告、方法、主要脚本、关键小表和精选图；不含大对象、表达矩阵、图像副本、参考/训练缓存或分类 scores 大表。全部 PDF/300 dpi PNG 留在每片 figures/；有图例的 PDF 另含 legend_only / panel_only。精选空间图与 ROI 的栅格组织底图不声称为全矢量。

## 图说明

- `umap_*`：全量 QC bins 的投影坐标，分别绘制 cluster 与两个参考候选；两片 UMAP 独立。
- `spatial_candidate_*`：原图配准下的全量候选；遮罩表示无 supplied 2 μm 输入范围。
- `ROI_candidates_*`：Cell8 BLCA、Bin8 BLCA、Bin8 HPCA、RCTD 第一成分。同一区域、同色板。RCTD 灰色为未拟合/拒判，不等同原细胞 Unknown。
- `ROI_raw_signals_*`：原始面板 UMI>0 的空间位置；并非细胞类型。灰轮廓始终是 v2 8 μm 归属范围，青色为原核。
- `stratified_comparison`：参考覆盖内分歧与双组 marker 共检出；区间只作描述。
- `author_normalization_comparison`：同一子集、同一交集基因、作者标签设置的归一化对照。
- `RCTD_statuses`：含未过 bin QC 和算法门槛者的完整 ROI 分母。
- `input_QC_coverage` / `unassigned_UMI_contribution`：原生输入过滤位置及此前未分配信号的贡献。

本轮结果目录新增约 2.69 GiB（含隔离 R 包和 review.zip）；review.zip 约 26.6 MiB。用户 R 库另外新增依赖约 2.7 MiB。原始/历史结果和冻结参考未复制进新目录。
