# 两张 Visium HD 全谱系细胞注释 v2.1 审阅材料

生成时间：2026-10-23T18:48:46（服务器本地时间）。本轮计算止于细胞注释。所有标签均为 **computational_proposal_pending_user_review**，不是人工确认或最终冻结。

## 交付与范围

- `cell_annotations.parquet`：全体 308,549 个 v2 保留细胞，一行一个 ID，含原 QC/候选/旧校准、SingleR 大类和细类完整精简诊断、新建议、状态、具体 Unknown 原因与审阅索引。
- `cell_annotations.compact.tsv.gz`：便于一般表格工具读取的精简副本。
- `ALL_CLUSTER_REVIEW.tsv`：所有全局簇和全部谱系局部簇的联合及逐切片证据/命名审阅总表。包含多数标签同时也包含完整组成，不应将多数标签赋给所有成员。
- `metadata/*_markers.parquet`：实测 marker 表，含平均 log1p10k、检出率、与同分析集合其余细胞的差、描述性 z/FDR；`nearest_cluster_marker_comparisons.parquet` 提供局部 PCA 最近竞争簇的表达和检出率差。
- `figures/`：两切片共享标签配色、透明背景；所有大类与各谱系已生成 known/data-driven marker 图及空间/表达图。全局 UMAP 只画既有 core；空间图包含全部保留细胞。局部图标明 informative 数量。
- `code/attach_annotations.py`：按 cell_id 对齐，默认只验证；指定新的输出路径才生成附加 v21 字段的对象，拒绝覆盖原文件。

原始 counts、v1/v2 对象、分割和坐标均未更改。本轮没有重新训练 scVI、分割、恢复图像、下载新的生物数据，也未做 MP10/胆固醇、空间邻域、通讯、CNV/恶性推断。

## 各切片大类建议

| proposed_major | Cervilla | D1_Johnson |
| --- | --- | --- |
| B | 23 | 1500 |
| Endothelial | 360 | 738 |
| Epithelial | 55782 | 46991 |
| Erythroid | 0 | 2 |
| Fibroblast | 5435 | 32859 |
| Mast | 77 | 240 |
| Mural | 477 | 925 |
| Myeloid | 1246 | 7465 |
| Plasma | 631 | 736 |
| T/NK | 150 | 4048 |
| Unknown | 33746 | 115118 |

## 全谱系分辨率和状态

| lineage | section | n | fine_supported | subtype_unknown | subtypes | states |
| --- | --- | --- | --- | --- | --- | --- |
| B | Cervilla | 23 | 1 | 22 | Memory B (1) | none |
| B | D1_Johnson | 1500 | 129 | 1371 | Naive B (96); Memory B (33) | Stress (10); Interferon (3) |
| Endothelial | Cervilla | 360 | 27 | 333 | Capillary EC (11); Arterial EC (10); Venous EC (5); Lymphatic EC (1) | Stress (4) |
| Endothelial | D1_Johnson | 738 | 74 | 664 | Arterial EC (41); Venous EC (26); Capillary EC (6); Lymphatic EC (1) | Stress (55); Interferon (7) |
| Epithelial | Cervilla | 55782 | 0 | 55782 | none | Cycling (904); Stress (9994); Interferon (121); Epithelial_basal (5822); Epithelial_luminal (20710); Epithelial_squamous (722) |
| Epithelial | D1_Johnson | 46991 | 0 | 46991 | none | Cycling (6169); Stress (15727); Interferon (6872); Epithelial_basal (22582); Epithelial_luminal (608); Epithelial_squamous (17352) |
| Erythroid | Cervilla | 0 | 0 | 0 | none | none |
| Erythroid | D1_Johnson | 2 | 0 | 2 | none | none |
| Fibroblast | Cervilla | 5435 | 2389 | 3046 | Matrix_CAF (2162); iCAF (222); rCAF (5) | Cycling (19); Stress (148); Interferon (9) |
| Fibroblast | D1_Johnson | 32859 | 18424 | 14435 | Matrix_CAF (18069); iCAF (313); rCAF (26); vCAF (16) | Cycling (307); Stress (4729); Interferon (982) |
| Mast | Cervilla | 77 | 0 | 77 | none | Stress (3) |
| Mast | D1_Johnson | 240 | 0 | 240 | none | Stress (6); Mast_activation (1) |
| Mural | Cervilla | 477 | 261 | 216 | Pericyte (176); SMC (85) | Stress (8) |
| Mural | D1_Johnson | 925 | 590 | 335 | SMC (308); Pericyte (282) | Stress (81); Interferon (7) |
| Myeloid | Cervilla | 1246 | 417 | 829 | C1QC TAM (379); SPP1 TAM (22); Classical Mono (6); FOLR2 TAM (3); Neutrophil (3); cDC2 (3); LAMP3 DC (1) | Cycling (5); Stress (57); Interferon (7); C1Q_program (207); SPP1_program (46); FOLR2_program (82); Myeloid_S100A8_A9_program (24) |
| Myeloid | D1_Johnson | 7465 | 3820 | 3645 | C1QC TAM (3299); SPP1 TAM (361); Classical Mono (131); FOLR2 TAM (15); LAMP3 DC (10); cDC2 (4) | Cycling (297); Stress (1412); Interferon (632); C1Q_program (2538); SPP1_program (817); FOLR2_program (151); Myeloid_S100A8_A9_program (2605) |
| Plasma | Cervilla | 631 | 0 | 631 | none | Stress (1); Plasma_IgA (318); Plasma_IgG (175) |
| Plasma | D1_Johnson | 736 | 0 | 736 | none | Stress (17); Interferon (1); Plasma_IgA (271); Plasma_IgG (96) |
| T/NK | Cervilla | 150 | 1 | 149 | Treg (1) | Stress (3) |
| T/NK | D1_Johnson | 4048 | 142 | 3906 | CD4 Naive/Memory-like (71); Treg (30); Tfh/CXCL13-like (23); NK (11); CD8 GZMK/Effector-memory-like (5); CD4 Helper/Conventional (2) | Cycling (2); Stress (88); Interferon (4) |

Epithelial、Plasma、Mast 也重新检查了谱系内表达结构及状态；参考对这些群体本来只有 broad 标签。没有将上皮 basal/luminal/squamous、浆细胞 IgA/IgG 或 mast 活化程序自动升格为新身份。保留大类而细类 Unknown 是可用的父类结论，不代表未分析。D1 的2个 Erythroid 建议属于参考之外的实测血红蛋白程序对象，必须结合 artifact flag 审查，不据此宣称发现成熟红细胞核或新群体。小群体不强制聚类。

## 实际表达发现与需要区分的解释

- 上皮：Cervilla 的 luminal 程序较多，D1 的 basal/squamous 程序较多；程序可以重叠。Epithelial::0 富集 KRT5/KRT14/KRT6A，::5 有 OLFM4/HPGD/MUC1/PIGR，::2 主要由低深度/MT 差异区分。保留 Epithelial 身份并另列状态，不推断恶性或机制。
- Fibroblast/CAF：Matrix_CAF 有 COL11A1、MMP11、POSTN、CTHRC1/FAP 等支持，而不是仅用 COL1A1。局部 ::1/::3/::5 的数据驱动差异混有 KRT/上皮程序，::6 带 VWF/RGS5/血管信号，不能把这些局部簇当成纯 CAF 新类型。逐细胞竞争转录本另列 secondary_lineage_signal。
- 髓系：Myeloid::5 是 D1 的43个细胞，CHIT1/CHI3L1/APOE/SPP1/GPNMB 同向支持 SPP1 相关巨噬程序；更多 C1Q 细胞分布于其他簇。::2/::3 有 S100A8/A9 和上皮差异，S100A8/A9 本身不决定 Classical Mono，需 FCN1/CCR2/S100A12 等更具体支持；炎症程序另存 state。C1Q/SPP1 并存保留，不把缺乏某个 marker 当互斥证据。
- T/NK：::0 有 CXCL13、SELL/GIMAP7，::1 有 CD3D/E，但并非整簇均可定 Tfh 或同一 T 细类；::2 受 Ig 信号、::3 受炎症/上皮信号影响。多数保留 T/NK，具体 Treg/Tfh/naive-memory/NK 等建议仅对应有证据的成员。
- B：两个局部群主要是 BANK1 等 B 程序与 Ig 丰富背景差异，少部分支持 Naive/Memory；本轮没有足够证据提出 GC B。极少 Cervilla B 和个位数细类不构成稳定群体的证明。
- 内皮：局部两簇并未自然分成四种血管身份，背景差异明显；只对有区分 marker 的成员提出 Arterial/Venous/Capillary/Lymphatic EC，其余保留 Endothelial。
- Mural：RGS5/CSPG4/PDGFRB 与 MYH11/CNN1/MYOCD/LMOD1 用于区分 Pericyte/SMC；局部两簇仍有上皮/Ig 背景差异，不能凭簇号等同于这两个身份。原 Fibroblast 到 Mural 的纠正与剩余 CAF 均分别审阅。
- Plasma：两簇有 IgG/IgA 相关表达差异，但只命名 Plasma 并单列免疫球蛋白状态；游离 Ig 表达不能单独建立 Plasma 身份。
- Mast：两个小表达簇以一般表达差异为主，没有可靠细身份分离；保留 Mast。活化程序支持极少，不能由 marker 缺失推断其他细胞静息。
- 参考外红系信号：仅2个 D1 对象满足协同血红蛋白程序，均带 artifact_review_flag；这是待核查对象，不是已确认红系群体。没有支持独立 Neural_glial 群体；不能把稀疏信号缺失写成生物学不存在。

这些命名是可修订的计算建议。主要身份有支持但另有竞争转录本的细胞，可能受环境 RNA 或重建边界影响；本轮不能确定其来源。局部 HVG/PCA 揭示了这些差异，但表达聚类本身并不证明新的细胞身份。

## 与 v2 的变化

| section | n | old_major_unknown | new_major_unknown | old_display_unknown | old_fine_supported | new_fine_supported | major_changed | new_mixed_flag | artifact_review |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 97927 | 28210 | 33746 | 28129 | 1822 | 3096 | 26314 | 1137 | 38 |
| D1_Johnson | 210622 | 67680 | 115118 | 67513 | 11635 | 23179 | 88698 | 17811 | 163 |

| section | old_predicted_Myeloid_to_T | old_predicted_Myeloid_to_B | old_predicted_Fibroblast_to_Mural | C1QC | SPP1 | C1Q_and_SPP1_measured_programs | fine_supported_old_calibration_failed |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 16 | 0 | 265 | 379 | 22 | 17 | 1210 |
| D1_Johnson | 574 | 104 | 560 | 3299 | 361 | 578 | 11312 |

完整转移见 `metadata/major_transition_v2_to_v21.tsv` 和 `label_transition_v2_to_v21.tsv`。old_predicted 的回归检查与 old_major_label（v2 已接受标签）是不同基线，勿混用。Mural 从旧 Fibroblast 父类拆出是层级修正的一部分。C1Q 与 SPP1 程序允许在同一对象共现，两个状态不互斥；提出的 TAM 标签不能被解释成互斥稳定机制。旧校准失败可被实际 marker/表达支持覆盖，但这并不消除不确定性。

## Unknown 与质量的具体原因

| unknown_reason | Cervilla | D1_Johnson |
| --- | --- | --- |
| broad_identity_supported_finer_identity_not_established | 56490 | 47969 |
| competing_lineage_programs | 1137 | 17811 |
| competing_related_subtype_programs | 1 | 12 |
| low_information_identity_unresolved | 12494 | 29867 |
| major_supported_subtype_information_insufficient | 4594 | 24344 |
| no_identity_anchor_detected | 9260 | 22402 |
| not_applicable | 3096 | 23179 |
| possible_segmentation_issue | 6 | 2 |
| reference_external_signal_insufficient | 526 | 1940 |
| sparse_or_competing_expression | 10323 | 43096 |

`major_unknown_reason` 对大类未定单列原因；大类清楚但细类不清保留父类。low_information、v2 mixed、当前竞争程序和 artifact flag 均独立保存，彼此可重叠，不可相加当互斥细胞数。当前 mixed 仅表示同对象可比的跨谱系表达，不等于已证明 doublet。血红蛋白背景、宽泛 ECM 和可能分割问题需要区分；没有为减少 Unknown 改变分母。低信息总数仍为 Cervilla 12,573、D1 30,342。

## 全簇审查与需集中判断的问题

共 72 个联合 review 条目：全局 34 个（33 个表达簇 + 一个明确标出的 noncore 集合），谱系 38 个。low_information_unclustered 是审阅集合而非表达簇。所有条目均有逐切片支持；细胞很少或某切片没有成员时不制造空群。

| cluster | lineage | n | suggested_dominant_label | review_question |
| --- | --- | --- | --- | --- |
| Myeloid::3 | Myeloid | 1689 | C1QC TAM | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Fibroblast::1 | Fibroblast | 10351 | Matrix_CAF | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Fibroblast::6 | Fibroblast | 1036 | Matrix_CAF | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份 |
| Myeloid::4 | Myeloid | 213 | Myeloid | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份；大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类 |
| Epithelial::0 | Epithelial | 24959 | Epithelial | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Mast::1 | Mast | 98 | Mast | 主要身份外仍有成套竞争转录本；局部簇可能由背景/边界混入驱动，勿将其差异直接视为新身份 |
| Epithelial::2 | Epithelial | 20303 | Epithelial | 低深度/线粒体主导；保留可支持的大类，勿作为新身份 |
| Erythroid::small_group | Erythroid | 2 | Erythroid | 无满足描述性筛选的特异富集基因或样本太少；不由此创造新身份；仅2个参考外血红蛋白程序对象；优先核查非核/边界污染，不能认定稳定红系群体；大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| T/NK::low_information_unclustered | T/NK | 2 | T/NK | 低深度/线粒体主导；保留可支持的大类，勿作为新身份；无满足描述性筛选的特异富集基因或样本太少；不由此创造新身份；大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| B::0 | B | 929 | B | 大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Endothelial::0 | Endothelial | 552 | Endothelial | 大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类 |
| Mural::0 | Mural | 1032 | Mural | 核对 RGS5/PDGFRB 与 MYH11/CNN1 的 pericyte/SMC 分界及 fibroblast 反证；大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类 |
| B::1 | B | 594 | B | 大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Mast::0 | Mast | 219 | Mast | 核对建议命名与所列 marker；当前无明显专门阻塞 |
| Plasma::0 | Plasma | 792 | Plasma | 核对建议命名与所列 marker；当前无明显专门阻塞 |
| Endothelial::1 | Endothelial | 546 | Endothelial | 大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类 |
| T/NK::3 | T/NK | 303 | T/NK | 大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类；来源偏倚；目前每来源一张，不能拆分技术与生物差异 |
| Plasma::1 | Plasma | 575 | Plasma | 核对建议命名与所列 marker；当前无明显专门阻塞 |
| Mural::1 | Mural | 370 | Mural | 核对 RGS5/PDGFRB 与 MYH11/CNN1 的 pericyte/SMC 分界及 fibroblast 反证；大类可支持，细亚群 marker 稀疏或相近程序竞争；保留父类 |

全局 1/11/12/14/18/21/25 等含不同谱系的成员，不能凭一个多数 marker 整簇覆盖。多个小全局簇的 top genes 主要为 MT 基因且深度低，已明确标出，不命名为新身份。`ALL_CLUSTER_REVIEW.tsv` 的 competing_lineage_markers 必须结合细胞级 mixed 标记读取：整簇出现不同谱系不等于每个成员混合。上表仅为集中查看入口，全部谱系和有把握簇也在总表中。

## 方法及解释边界

1. 使用原 v2 重构细胞 raw integer X。marker 使用全部输入基因，映射符号并合并重复后 18,082 个；参考不包含的实测基因仍可用于 marker 检查。实测检出率和 log1p(10,000/细胞总 counts) 的差为描述性证据，不以 scVI 去噪/插补表达发现 marker。
2. 原共享 scVI 涉及 **2 sections**，core 252,601（81,789 + 170,812）；本轮只复用其 latent20、全局簇和 core UMAP。大类首先来自身份 marker/竞争证据，允许跨旧标签纠正；弱成员可通过 marker 支持种子的表达近邻投射，具体阈值、距离范围和证据路径见 `major_rules.json`、`major_projection.json`。没有物理邻居标签传播。
3. 每个有足够细胞的纠正后谱系重新选 section-aware HVG（至多 2,000，去 MT/RPL/RPS，补充实测身份 marker），稀疏归一化和方差缩放，PCA30、表达 kNN15、Leiden0.65，seed123；小群体仅做 marker 审阅。未强制来源混合。仅每来源一张切片，来源/技术与生物差异无法完全分离。局部 marker 以同谱系为对照，并记录最近局部 PCA 竞争簇。
4. SingleR **确实运行**：R4.5.3、SingleR2.12.0，原冻结 `reference_celltype_v2` 32 类，原参考300,033细胞按标签及 dataset/specimen 轮转选至多400/类，12775参考细胞（vCAF全375）。specimen不假定独立患者。共享14,327去重基因、实测 counts→log1p10k；wilcox/de.n20，aggr.ref=TRUE，seed123，8线程；大类9类（SMC/Pericyte归Mural）和细类32类分别训练一次，所有 query 按5000分块逐细胞预测，按整切片统一prune。保存 candidate/pruned/selected score/delta/最初第一第二名；score/delta 是相似性诊断，不是正确概率。within-major initial 候选来自同一细类初始分数，不称为额外分层微调运行。安装缺失 scrapper1.4.0 及其依赖；服务器 R 启动路径优先级使最终依赖位于现有用户库，未重装项目环境，日志和 sessionInfo 保存。
5. 复用 CellTypist，不重训全参考。两工具同参考一致只作为算法间支持；任何单工具或旧校准均不是真值。总表的 major agreement 是原始标签精确匹配，旧 CellTypist 的 Fibroblast 与新 SingleR 的 Mural 层级差异也会计为不一致，因此不能将该比例直接理解为生物学冲突。直接 marker、特异 anchor、参考兼容和表达簇支持是不同路线，没有全细胞通用的双工具一致/旧校准通过/固定 marker 数/旧候选纯度 AND 门槛。实现仍有显式可审查的描述性阈值，未声称通过独立 HD 真值校准。
6. 亚群证据要求相近身份的区别；近似竞争时保留父类。state 是实测程序，不是机制或恶性结论。单 marker 未检出不是不存在证据。GC B 不预设存在；单独 BCL6/cycling 不足以定 GC。
7. 大量细胞不能替代生物学重复，z/FDR 仅用于排序和描述；nearest 比较主要看效应与检出率。分割核+局部扩张不是精确膜边界；本轮复用 v2 已有形态审阅，8 个原形态疑问仍留表，没有新人工确认。原 H&E/分割例图位于 `../cell_annotation_v2/{Cervilla,D1_Johnson}/figures/representative_large_nuclei.png`。

官方实现依据：[SingleR 单细胞参考](https://www.bioconductor.org/books/release/SingleRBook/sc-mode.html)、[参考聚合与复用训练](https://www.bioconductor.org/books/release/SingleRBook/advanced-options.html)、[诊断和剪枝](https://bioconductor.org/books/release/SingleRBook/annotation-diagnostics.html)。

## 复现和必要检查

Python：`../cell_annotation_freeze_v1/.venv/bin/python`，版本见 `metadata/python_method_versions.json`；R版本/依赖见 `metadata/SingleR_session.txt`。依次运行 code 下 prepare_inputs.py → run_singleR.py；expression_evidence.py → major_proposals.py → lineage_structure.py → subtype_proposals.py → review_tables.py → make_review_figures.py → validate_and_report.py。SingleR 使用现有模型时固定训练不变；全部路径从脚本目录解析。

必要一致性检查 22 项通过，见 `metadata/delivery_checks.tsv`：ID 唯一与完整、按 ID 对齐、两张分母正确、low_info/core 保留、空间坐标未改、源文件大小/修改时间未变、标签层级、Unknown 原因、SingleR 实际全覆盖、全簇/全谱系审阅覆盖。没有额外全文件哈希或多份完整表达矩阵副本。

本轮停在集中审阅材料；建议标签尚未经用户确认，不进入下游分析。

本轮三个可再生 SingleR H5 桥接文件已清理（约330 MB），重算时先运行 prepare_inputs.py。保留模型、抽样 ID、所有最终标签及 marker/低维证据。图件索引：metadata/figure_index.tsv，共 81 张透明 PNG。
