# BLCA Visium 肿瘤与组织区域计算注释 v2

已完成 **63/63 张切片**的逐 spot 区域注释和核心图，覆盖 **265,002 个原始位置、7 个 dataset**。SpaCET 完成 63 张，BANKSY 完成 63 张，实际 CNA profile 支持完成 63 张。

**61 张**存在可评估肿瘤—非肿瘤界面；**55 张**达到内侧、外侧、内部各至少 20 个 spot 的几何可用条件且未因高不确定性/坐标精度降级。其中 **19 张**同时满足较高的计算证据一致性条件，可优先进入下一轮比较。它们仍是计算推断，不是病理确诊或已确认的侵袭前沿。

| dataset_id | 完成切片 | 推断肿瘤spot | 可评估界面切片 | 不确定spot | 原Epi门槛下新增肿瘤spot |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | 13 | 34974 | 13 | 6221 | 2849 |
| GSE171351 | 4 | 1747 | 4 | 316 | 404 |
| GSE238145 | 10 | 32358 | 10 | 2604 | 12186 |
| GSE285715 | 6 | 4841 | 4 | 7789 | 9 |
| GSE319536 | 16 | 18667 | 16 | 5486 | 3590 |
| Inhousedata__UNVERIFIED | 7 | 13132 | 7 | 1457 | 2773 |
| NPJ__UNVERIFIED | 7 | 13408 | 7 | 1497 | 112 |

## 实际区域与主要差异

| 组织身份 | spot 数 |
| --- | --- |
| Tumor | 119127 |
| Other tissue | 69556 |
| Uncertain/conflicting | 25370 |
| Immune-rich region | 19881 |
| Non-tumor stromal region | 16822 |
| Smooth-muscle-rich region | 14010 |
| Normal urothelium | 236 |

新主 mask 包含 **119,127** 个推断肿瘤 spot。其中 **21,923** 个在旧 Epi≥0.25 mask 之外；旧 Epi mask 中 **46,257** 个位置未进入新肿瘤 mask。具体到每片的 Jaccard、新增/移除数量见 [切片汇总](tables/section_region_summary.tsv)。这不是调整上皮门槛，而是改用恶性模型、实际 CNA 模式和组织上下文。

正常尿路上皮仅在有来源和分化表达/形态上下文的局部作推断注释；不会由 diploid 或低 CNV 自动推出。N1 的表面上皮带保留了与 SpaCET 高恶性估计的冲突标记，不能用这次上下文辅助重标记声称特异性验证成功。平滑肌标签基于 SMC 组成和 MYH11、ACTG2/CNN1 等多基因证据；没有自动称为肌层。未强行输出坏死或 TLS。

## 优先查看的图与具体冲突

- [Inhousedata_T4413900](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4413900/01_core_regions.pdf)：大片肿瘤区与 CNA 支持较一致，邻近免疫丰富组织和局部混合区可对照；内部小缺口未通过填孔强行消除。
- [GSE238145_B13](by_dataset/GSE238145/GSE238145_B13/01_core_regions.pdf)：展示不以旧 Epi gate 为入口后恢复的大片推断肿瘤区；同时看按染色体排序的真实 CNA profile。
- [NPJ_sample1](by_dataset/NPJ__UNVERIFIED/NPJ_sample1/01_core_regions.pdf)：较直观地展示肿瘤与免疫丰富、平滑肌成分丰富及混合非肿瘤组织的空间关系。
- [CanDis_2026_Sample8](by_dataset/CanDis__UNVERIFIED/CanDis_2026_Sample8/01_core_regions.pdf)：原生 CopyKAT 全 diploid，但 SpaCET 和部分连续 CNA 模式提示肿瘤。保留推断肿瘤，同时明确跨方法冲突；不能称为 DNA 验证。
- [GSE285715_T3](by_dataset/GSE285715/GSE285715_T3/01_core_regions.pdf)：原生低置信 diploid、基线问题及表达模型/CNA 分歧集中；局部不确定较多，边界解释优先复核。
- [GSE285715_N2](by_dataset/GSE285715/GSE285715_N2/01_core_regions.pdf)：来源为癌旁，大片 H&E 肌束同时表达 MYH11/ACTG2/CNN1 且 SMC 占优；SpaCET/CopyKAT 却提示恶性。主组织身份经上下文修正为平滑肌成分丰富区，并保留模型冲突。残余肿瘤候选不列入优先边界分析。
- [GSE285715_N3](by_dataset/GSE285715/GSE285715_N3/01_core_regions.pdf)：癌旁来源与大片上皮区的恶性模型估计冲突，缺少足够证据判定正常或肿瘤，局部保留 uncertain；没有可评估主肿瘤界面。
- [GSE285715_N1](by_dataset/GSE285715/GSE285715_N1/01_core_regions.pdf)：癌旁来源的肌束与表面上皮带；正常尿路上皮是上下文支持的计算推断，仍有 SpaCET 冲突。

每张片的 `01_core_regions` 同时展示原始 H&E/位置、H&E 叠加、组织分类、SpaCET 连续成分、CNA 支持和主 mask 分带。`02_distance_native_domains` 展示有符号距离、SpaCET 原生 0.5 对照和 BANKSY 域；`03_CNA_chromosome_profiles` 是实际 CNA 热图；`04_final_tissue_regions` 是独立分类图。所有正式图均有 PDF 和 300 dpi PNG。

## 下一轮可用性

以下名单优先用于计算定义下的肿瘤内部/内侧边缘/外侧邻域比较；目标 MP10、Treg、EC、CAF 等信号没有参与参数或阈值选择。组成和 marker 已用于区域命名，之后不能把同一特征的富集称为独立验证。

GSE238145_B10, NPJ_sample1, CanDis_2026_Sample10, CanDis_2026_Sample12, CanDis_2026_Sample2, GSE238145_B11, GSE238145_B13, GSE238145_B14, GSE238145_B9, GSE319536_BC12, GSE319536_BC15, GSE319536_BC19, Inhousedata_T4413900, Inhousedata_T4463209, Inhousedata_T4471438, NPJ_sample2, NPJ_sample3, NPJ_sample5, NPJ_sample7

其余存在界面的切片仍可作为探索性或敏感性分析对象，但需按具体冲突分层；没有内部核心、有效外侧或完整证据的切片不应强行进入同一主比较。

| section_id | n_tumor | n_inner | n_outer | n_interior | uncertain_fraction | boundary_analysis_status |
| --- | --- | --- | --- | --- | --- | --- |
| GSE285715_N3 | 0 | 0 | 0 | 0 | 0.674 | no_assessable_tumor_interface |
| GSE285715_T3 | 1203 | 78 | 38 | 851 | 0.616 | review_required_high_uncertainty |
| GSE285715_T2 | 1518 | 82 | 51 | 1308 | 0.529 | review_required_high_uncertainty |
| GSE285715_T1 | 2089 | 272 | 184 | 1626 | 0.44 | review_required_high_uncertainty |
| CanDis_2026_Sample3 | 2726 | 56 | 46 | 2500 | 0.353 | eligible_computational |
| GSE319536_BC17 | 1655 | 172 | 124 | 1408 | 0.283 | eligible_computational |
| CanDis_2026_Sample4 | 2036 | 55 | 45 | 1852 | 0.253 | eligible_computational |
| GSE319536_BC02 | 1072 | 118 | 89 | 858 | 0.223 | eligible_computational |
| CanDis_2026_Sample6 | 992 | 71 | 68 | 860 | 0.216 | eligible_computational |
| CanDis_2026_Sample7 | 1253 | 184 | 162 | 1053 | 0.21 | eligible_computational |
| GSE171351_C1 | 490 | 287 | 280 | 203 | 0.184 | eligible_computational |
| GSE319536_BC11 | 2696 | 456 | 364 | 2055 | 0.183 | eligible_computational |

## 必须保留的解释边界

- SpaCET malignant fraction 是模型组分比例，不是恶性概率；实际签名来源逐片保存。来源计数：CNA:BLCA: 59 张, expr:BLCA: 4 张。
- CopyKAT 自动参考未外部验证与基线未收敛分别保留；连续模式、全局偏移和原生标签冲突分开记录。原生 `c1:diploid:low.conf` 共 3341 个位置，未并入无结果。
- BANKSY 分辨率 0.4/0.6 的试运行 ARI 为 0.484–0.841，T3 的细域更敏感；固定规则下仅 2/42,087 个试运行 spot 的主肿瘤身份改变。主 mask 稳定不意味着细域就是病理真值。
- 主分带来自最终主 mask 和真实六邻接。组织外缘、QC 缺口和 uncertain 邻区不构成生物学界面。两圈敏感性标签已保存。距离单位为对齐图像像素，另提供邻居间距标准化值；没有虚构微米尺度。
- 全部 `review_status` 均为无已记录的人工病理复核。助手图像对照与真实病理医师确认不同；不能把本结果称为已确认侵袭前沿。
- 未证明有充分坏死/无细胞区形态证据时，不按低 UMI 强行命名。未找到可用肿瘤界面的切片仍保留组织分类和原始位置。

## 文件入口

- [全队列逐 spot 核心表](tables/all_spot_region_labels.tsv.gz)；完整附加证据在同目录 parquet 和各切片目录。
- [每片区域/界面汇总](tables/section_region_summary.tsv)、[区域与分带计数](tables/region_band_counts.tsv)、[一圈/两圈比较](tables/one_vs_two_ring_counts.tsv)。
- [定义检查](tables/definition_checks.tsv)、[原生 CopyKAT 标签计数](tables/native_CopyKAT_label_counts.tsv)、[试运行 mask 敏感性](tables/pilot_mask_domain_sensitivity.tsv)。
- [方法、规则及原始文献](METHODS.md)。轻量 `review.zip` 包含核心 PDF、小表和方法，不含原始表达、完整 CNA 矩阵和软件库。

## Dataset 总览

- [CanDis__UNVERIFIED 分类总览](figures/CanDis__UNVERIFIED_overview.pdf)
- [GSE171351 分类总览](figures/GSE171351_overview.pdf)
- [GSE238145 分类总览](figures/GSE238145_overview.pdf)
- [GSE285715 分类总览](figures/GSE285715_overview.pdf)
- [GSE319536 分类总览](figures/GSE319536_overview.pdf)
- [Inhousedata__UNVERIFIED 分类总览](figures/Inhousedata__UNVERIFIED_overview.pdf)
- [NPJ__UNVERIFIED 分类总览](figures/NPJ__UNVERIFIED_overview.pdf)
