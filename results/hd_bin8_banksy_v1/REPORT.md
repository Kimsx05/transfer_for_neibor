# 8 μm bins BANKSY 匹配对照

| 切片 | 全量 QC bins | 原始 UMI | 原大类分歧 (双方可评价) | 四粗类分歧 (同分母) | 下降 (百分点) | λ0 / 0.2 / 0.8 簇数 | 推荐空间层 | 状态 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 268924 | 64443819 | 27.13% | 16.26% | 10.87 | 101 / 29 / 41 | λ0.2 主用，λ0.8 粗域对照 | 完成并实际加载验证 |
| D1_Johnson | 364805 | 131440767 | 40.26% | 28.81% | 11.45 | 101 / 31 / 40 | λ0.8 组织域，λ0.2 局部对照 | 完成并实际加载验证 |

两片均保留原有全部 QC bins、18,085 个输入基因及原 RNA counts/data。三分支为全量直接 PCA/聚类，非 sketch 投影。旧 cluster_projected、逐 bin SingleR、分割、归属和 RCTD 均未重跑或覆盖。区域解释与候选类型分层保存，未据此正式修订细胞身份。

## 主要结论

主表分歧率的双方可评价分母：Cervilla 260,214 bins，D1_Johnson 349,659 bins。按同一集合的 UMI 加权，粗分类前后分歧承载比例分别为 10.74%→6.78% 和 29.25%→26.38%；尤其 D1 的主要高表达分歧仍然保留。

**已有表达簇本来就有组织结构。** Cervilla 的旧簇可辨认上皮样巢及间质穿插；D1 可辨认大片组织区域、细胞密集带和局部聚集。旧图的空间结构不是本轮新产生的。四粗类使免疫内部、Fibroblast/Mural 等粒度差异减少，但上皮—间质分歧仍然存在；两参考仍分别保留。

**Cervilla：以 λ=0.2 作为类型相关表达结构的主要辅助层，λ=0.8 保留为粗区域对照。** 匹配 λ0/0.2/0.8 分别得到 101/29/41 簇。邻接同簇率由 51.91% 到 77.38% 再到 97.39%，只是预期平滑程度。λ0 的许多小簇是低信息量表达模式，均保留；λ0.8 不是简单“合并成更少簇”，而是产生更连续的区域划分。原始 UPK2/WFDC2/PHGR1/PSCA 支持部分上皮富集簇；LUM/SFRP2/BGN/TAGLN/POSTN 支持间质成分；局部 LYZ/CD74/C1QA/免疫球蛋白及 JCHAIN 与间质信号共同出现，宜描述为间质—免疫复合区域，不一概命名为纯免疫区。

Cervilla 的 far_UMI_hotspot 中 1,022/1,024 bins 落入 λ0.8 的簇 7。H&E 可见细胞岛与纤维材料交错；该簇仍有 55.1% bins 同时检出两枚以上上皮及胶原面板基因。两个聚合参考都预测 Fibroblast，也不能据此把这些 bins 全部视为纯间质细胞。该 ROI 既有 RCTD 中，201/1,024 个 QC bins 被标为 Epithelial/Fibroblast 的 doublet_certain（不分先后），可作为混合成分的相关证据；RCTD 与 BLCA 共享参考，不是独立真值。λ0.8 的边缘簇 21 覆盖 edge ROI 的 813/816 bins，小核 ROI 的 625/650 bins 归入簇 27；两者全簇 nFeature 中位数分别为 17、15。连续区域容易辨认，但细小结构与信息量差异被压缩，不能据此宣布低信息量身份得到解决。stroma ROI 中 λ0.8 簇 15 的局部 nFeature 中位数仅 39，低于该簇全片 208.5，提示“域平均解释”不能替代逐 bin 证据。

Cervilla 类型层保留 BLCA 的完整九类候选及 HPCA 对照：BLCA 对 Mural/Mast/Plasma 有额外参考覆盖，但本轮没有证明其在所有区域更正确；HPCA 也不能独立排除这些类型。上皮—胶原热点保留混合信号说明。Cervilla 的 T/NK、内皮面板整体较稀疏，不仅凭空间连续图或参考标签宣称出现了可靠的大片相应类型区域。没有 supplied bins 的 H&E 区域继续单列，不归入任何 domain。

**D1：以 λ=0.8 作组织区域层，配合 λ=0.2 检查局部类型相关结构；HPCA、BLCA 继续独立保留。** 三分支为 101/31/40 簇，邻接同簇率为 58.94%/78.87%/96.23%。空间结构更连续，但没有证明类型准确率提升。原始 KRT14/KRT6A/KRT6B/KRT17/COL17A1 等支持主要分歧域中的上皮成分，COL1A1/COL1A2/FN1/POSTN 等与之共存。λ0.8 簇 1、2 的上皮—胶原双面板检出率分别约 84.0%、84.6%，其聚合 BLCA/HPCA 仍分别为 Fibroblast/Epithelial。BANKSY 帮助把争议落到可复核的区域，并未“纠正”为纯上皮或纯成纤维。

D1 主要方向的冲突共有 **37,972 bins（全 QC 的 10.41%），承载 23,967,631 UMI（18.23%）**。其中 33,421 bins（88.02%）检出至少 100 个基因；至少 500 基因的 20,116 bins 中，98.11% 同时检出两枚以上上皮及胶原面板基因。只有 113/37,972（0.30%）位于 supplied 输入边界 16 μm 内，这排除了“主要由输入边缘/低检测量造成”的解释，但输入边缘不是上皮—间质解剖边界。相反方向仅 1,429 bins，其中 1,186（83.00%）nFeature<100，信息量构成明显不同，不能对称理解为同一种错误。

主要冲突广泛位于高信息量的上皮相关/上皮—间质混合表达域，而不是仅在组织外缘。新增 conflict_focus 位于原图物理坐标 x≈7,040–7,296 μm、y≈11,968–12,224 μm；H&E 可见连续上皮样细胞结构和小腔隙，不能解释为单纯薄层间质边界。该 ROI 的 495/1,024 bins 为主要冲突，其中 98.79% 有双面板多基因检出。它们在匹配 λ0 中已主要分属两个表达簇（645、329 bins），λ0.2 主要汇入簇 1（930），λ0.8 主要汇入簇 1（1,018）。因此不是 BANKSY 首次发现表达结构，而是区域表示更连续。既有 stroma ROI 确实有上皮样结构与纤维带交界；目前没有全片人工病理分区，不能定量声称所有冲突都在上皮主体，更不能将某一参考判为真值。

综合而言，D1 更支持**高信息量混合表达加参考响应差异**，不支持把全部冲突归咎于旧核外扩归属、随机低 UMI 或组织边缘。固定原生 bins 中仍有相同分歧，也不能证明更换单位已解决来源问题。仍无法区分 bin 内多种细胞/结构贡献、局部 RNA 混入与单个细胞自身的 ECM 相关表达；本轮不推断 EMT、恶性或 CNV。HPCA 的 epithelial 候选在重点 ROI 有形态及角蛋白支持，但不能据此覆盖 BLCA，也不能把 BLCA Fibroblast 直接解释为 CAF。

D1 的小核 ROI 在 λ0.8 中保留簇 15 的免疫相关成分（339/1,002 bins；MS4A1/CD37/IGHM/LTB、两聚合参考 B），同时有簇 11 的胶原富集成分（554 bins）。这是区域层的复合组织结构，不是每个 bin 都能可靠细分。该 ROI 原 RCTD 有 305 个 QC bins 未达到 100 UMI 门槛，已有拟合可见 T/NK、B、Plasma 等成分；它与 BLCA 共享参考。stroma ROI 则有 662/1,007 个 QC bins 因同一 RCTD 门槛未拟合，不能只拿剩余高信息量部分证明整个 ROI 已解决。新增 conflict_focus 未拟合 RCTD，明确为 not_tested。far_UMI_hotspot 的 840/840 bins 归入 λ0.8 簇 17，原始上皮相关信号与聚合参考一致，但 H&E 仍含碎片及腔隙，不能把所有信号点视为完整细胞。边缘 ROI 主要为低信息量簇 20/3，强空间连续性也没有修复这一限制。

**下一步建议与边界。** 可基于本轮保留的类型候选与区域层，设计经过核验的区域表达比较和精细亚群分析；先复核 D1 高信息量上皮—胶原域、Cervilla 胶原热点，以及低信息量边缘/薄结构。正式细胞身份修订、通路/代谢评分、邻域与通讯应等待对这些证据的评估，不把混合 bins 或域内所有 bins 当成一种细胞。两片均保留 λ0 基线和全部原始 counts/候选，便于后续区分表达与空间效应。

**剩余实现限制：** k=50 在少数稀疏边缘仍可能跨块内孔洞，距离尾部已标记；分组仅阻断独立 supplied 组织块。D1 UMAP 近似邻居检索 recall 约 58–69%，因此细小 UMAP 分隔不作为主结论。没有人工真值，不宣称准确率、纯度或 Unknown 得到改善。各域解释是依据表与重点 ROI 复核的暂定解释，不是全片人工病理真值。

**未完成分支：无。** 两片三种 λ、四粗类对照、λ0.2/0.8 的 cluster-level SingleR、原始 marker 汇总、固定 ROI/既有 RCTD 对照、增量加载和实际重读验证均已完成。未重跑逐 bin SingleR、RCTD、分割或本轮范围外分析。


## 分母与粗分类

主表分母是同片两套参考均有候选的 bins。全 QC 和原 HPCA 覆盖类别限定的分母另存 `tables/coarsening_all_denominators.tsv`，不可混用。原生大类→四粗类只是解释粒度改变；Mural+Fibroblast 为 Mesenchymal，Mast/Plasma 等归 Immune 但原字段保留。HPCA 缺少 Mural/Mast/Plasma 的独立覆盖，不构成对这些类型不存在的证据。Unknown 是缺失/不可评价的候选；低信息量区域、未测试 RCTD 和没有输入覆盖另行区分。SingleR pruning 保留不等于高置信度。

### Cervilla

独立原始/剪枝后粗类的 bin 数、比例及 UMI 占比：`Cervilla/coarse_composition.tsv`；交叉表：`Cervilla/coarse_crosstabs.tsv`。

| branch | clusters | adjacent_pairs | same_cluster_pairs | adjacent_agreement | connected_fragments |
| --- | --- | --- | --- | --- | --- |
| historical | 20 | 517600 | 287881 | 55.62% | 50023 |
| l00 | 101 | 517600 | 268688 | 51.91% | 92036 |
| l02 | 29 | 517600 | 400495 | 77.38% | 28919 |
| l08 | 41 | 517600 | 504097 | 97.39% | 2111 |

相邻网格同簇率和碎片数是描述性的空间连续性指标；BANKSY 本来就鼓励连续，不能把这些值解释为准确率。旧 sketch 聚类与新匹配 λ0 的实现/规模/resolution 均不同，旧→新差异不能全部归因于空间信息。

空间邻域：24 个 supplied in-tissue 连通块；跨块边为 0；211 个 bins 因所在块较小使用少于 50 个邻居。第 50/实际最远邻居的中位距离 32.98 μm，99% 分位 58.23 μm，最大 201.89 μm。块内孔洞/稀疏边缘仍可能被跨过，逐 bin 距离已交付；该尾部不是新 QC 删除标准。分组只限制邻接，不改变全片 scaling。

### D1_Johnson

独立原始/剪枝后粗类的 bin 数、比例及 UMI 占比：`D1_Johnson/coarse_composition.tsv`；交叉表：`D1_Johnson/coarse_crosstabs.tsv`。

| branch | clusters | adjacent_pairs | same_cluster_pairs | adjacent_agreement | connected_fragments |
| --- | --- | --- | --- | --- | --- |
| historical | 16 | 707325 | 435616 | 61.59% | 51527 |
| l00 | 101 | 707325 | 416879 | 58.94% | 95760 |
| l02 | 31 | 707325 | 557888 | 78.87% | 30618 |
| l08 | 40 | 707325 | 680631 | 96.23% | 2461 |

相邻网格同簇率和碎片数是描述性的空间连续性指标；BANKSY 本来就鼓励连续，不能把这些值解释为准确率。旧 sketch 聚类与新匹配 λ0 的实现/规模/resolution 均不同，旧→新差异不能全部归因于空间信息。

空间邻域：11 个 supplied in-tissue 连通块；跨块边为 0；23 个 bins 因所在块较小使用少于 50 个邻居。第 50/实际最远邻居的中位距离 33.00 μm，99% 分位 56.02 μm，最大 351.59 μm。块内孔洞/稀疏边缘仍可能被跨过，逐 bin 距离已交付；该尾部不是新 QC 删除标准。分组只限制邻接，不改变全片 scaling。

## 证据与文件

- `cluster_summary.tsv`、`cluster_reference_composition.tsv`：所有历史/三路簇的 counts、信息量、marker 检出、参考组成及边缘关系。
- `cluster_markers.tsv.gz`：全原始基因表达上的描述性效应/检出差，非只用 HVG；`panel_marker_statistics.tsv.gz` 保留成组 marker 均值和比例。
- `cluster_SingleR.tsv`、`cluster_SingleR_scores.tsv`：λ0.2/0.8 的原始表达聚合对照，保留两参考分歧；没有给整簇强行补同一种 bin 身份。
- `cluster_domain_interpretation.tsv`、`reviewed_region_overrides.json`：原 cluster ID、证据、暂定区域名及是否经过重点图像复核。未重点复核的小簇仍标为摘要级解释。
- `conflict_information_summary.tsv`、`conflict_cluster_distribution.tsv`、`conflict_ROI_summary.tsv`、`ROI_cluster_summary.tsv`：D1 正反方向分歧及一致区的固定分母对照。
- `fixed_rois.tsv`：沿用 ROI；D1 新增冲突热点的选择依据和坐标在 `additional_roi_selection.json`。新增热点未重跑 RCTD；原 RCTD 与 BLCA 共享参考，不是独立真值。
- `metadata.tsv.gz` 与 `l00/l02/l08_embeddings.rds`：增量结果，三路 ID 与坐标严格对齐；不含复制的表达矩阵。原详细候选、scores、marker 支持继续从源 Seurat 获取。

主要图：`figures/spatial_clusters_historical`、`spatial_coarse_references`、`spatial_clusters_l00/l02/l08`、`spatial_interpreted_three_lambda`、`spatial_recommended_regions_HE`、`ROI_types_*`、`ROI_clusters_*`、`ROI_raw_markers_*`；D1 另有 `D1_reference_conflict_location`、`D1_conflict_region_composition` 与 `D1_conflict_raw_marker_evidence`。图内无大标题；H&E/高密度点层为栅格，文字图例为矢量；白底 PNG 300 dpi，含图例图均有全图/主体/图例 PDF。

## 加载与复跑

依赖原对象：`../hd_bin8_comparison_v1/Cervilla/seurat.rds`、`../hd_bin8_comparison_v1/D1_Johnson/seurat.rds`；外部图像继续由原加载器解析。加载不会保存重复的完整 Seurat。

```r
root <- '/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_bin8_banksy_v1'
source(file.path(root, 'code/load_result.R'))
x <- load_banksy(file.path(root, 'D1_Johnson'))
table(x$banksy_l08_region)
table(x$coarse_blca, x$coarse_hpca)
SeuratObject::Embeddings(x, 'banksy_l02_umap')
SeuratObject::LayerData(x, assay='RNA', layer='counts')
```

`bash rerun.sh` 顺序复用已完成分支，检查固定参数，更新证据汇总和最后报告。软件版本在各片 `model_versions.json`、`sessionInfo.txt`；完整设置及与官方示例的区别见 `METHODS.md`、`parameters.json`。只进行一次最终加载验证，结果在各片 `validation.json`，并从实际重读对象绘制 `figures/reloaded_attachment`。

本轮未扩展 MP10、精细亚群、正式邻域生态、通讯、整合或细胞身份修订。后续先评估本轮组织区域与候选证据；争议区域仍需病理/细胞边界层面的针对性核验。

review.zip 精简收录关键 PNG、必要表格和代码；完整 PDF（含主体/图例版本）及其余图件保留在服务器结果目录。


本轮目录实际新增文件总量（含 review.zip）：0.851 GiB；review.zip：37.8 MiB。未复制 counts、完整 Seurat/h5ad、参考及全片图像。
