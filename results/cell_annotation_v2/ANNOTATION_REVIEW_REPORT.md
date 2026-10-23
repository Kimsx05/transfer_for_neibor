# 两张 Visium HD 细胞重建与亚群注释 v2 审阅版

2026-10-23。已分别完成 Cervilla、D1_Johnson 的新图像核分割、真实2µm原始计数聚合、宽松QC、冻结参考候选预测与本轮深度校准，再完成两张共享scVI、表达聚类和marker审查。旧v1及原始输入保留。本目录是可修改的审阅版，没有新设冻结文件；本轮计算止于注释。

## 实际细胞数量与处理范围

| section | 图像核 | 有非零表达的重建细胞 | QC保留 | 低信息 | 低信息占比 | 原始counts分配率 | 训练与UMAP核心 | 未映射但保留 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 173057 | 107071 | 97927 | 12573 | 12.84% | 70.57% | 81789 | 16138 |
| D1_Johnson | 229392 | 227272 | 210622 | 30342 | 14.41% | 81.61% | 170812 | 39810 |

“重建细胞”仅统计分配后counts>0的对象。另有Cervilla 232、D1 1316个核虽被分配到提供的bins，但合计counts=0，未称为有表达细胞；它们在排除元数据中。未分配任何bin的核只保存几何和原因。QC排除对象分别为 9376、17966，包括上述零表达对象。

| section | 输入bins | 图像未覆盖bins | 分配bins | 未分配bins | 输入counts | 分配counts | 未分配counts | QC保留counts |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 4722911 | 1274 | 2996804 | 1726107 | 64431796 | 45469172 | 18962624 | 45417647 |
| D1_Johnson | 6362240 | 0 | 4892217 | 1470023 | 131460620 | 107287385 | 24173235 | 107203097 |

Cervilla使用恢复的11418×8964同片RGB，原生0.8636122748 µm/px；D1使用27318×26836捕获区RGB，原生0.2478437293 µm/px，原fullres裁剪 `[23764,35851,51082,62687)`。D1局部坐标仅减一次裁剪原点。模型mpp约0.5，重采样不会创造原图没有的细节，也不扩大表达测量视野。模型像素中心变换为 `(fullres−cropOrigin+0.5)×scale−0.5`；物理坐标为原fullres×原生mpp。四套坐标分别保存于 `spatial`、`spatial_image`、`spatial_model`、`spatial_um`；详见每张 `metadata/transform.json`。没有counts的图像区域不能解释为零表达组织。

## 分割、扩张与QC

StarDist `2D_versatile_he`，prob=0.05、nms=0.3；先查看致密、间质/混合与组织边缘ROI，再分块2048推断。D1初选edge实际仍位于组织内部，目视检查后改为跨组织边界的窗口。各张 `pilot_review.json`、`segmentation_pilot.png` 保留实际审查。Bin2cell最大扩张2个2µm array bins（4µm），种子123、PCA30用于冲突处理；没有GEX-only补细胞。H&E核及扩张只是细胞近似，并非实测细胞膜。

正式X来自原始Space Ranger整数counts逐bin求和。destripe只影响临时冲突处理视图，汇总前恢复原始矩阵；没有将归一化值取整为counts。完整18,085个输入features均保留，HVG只用于模型视图。逐基因分配计数守恒；未分配counts保持未分配。

保留条件为n_counts≥10且n_genes≥10；10–30 genes单独标记但保留。未套200/500 genes或统一线粒体删除阈值。5MAD counts/genes/面积、触图像边缘等仅标记并影响表达核心选择，不单项删除。面积与counts未作为自动删除条件；实际代表图中确认的4个D1非单核血样聚集区假分割被排除（均原本不在scVI核心），另8个大而不规则/多叶轮廓保留但标为形态歧义。Cervilla其中1个在审阅前已进入表达核心，历史training标记保留，其最终身份不作为确定标签。代表形态审查：两张代表窗口的核与bin分配没有明显全局坐标偏移，扩张局部且保留未分配空隙。Cervilla仍受0.864µm原生像素及拥挤/细长轮廓限制；D1边缘血样聚集/碎片区域出现明确非核假分割。本次4个逐一确认对象排除，另8个大而不规则轮廓保留为歧义；这是代表对象审查，并非已穷尽全片血液/碎片区伪分割。

## 参考、证据与建议标签

继续使用 `cell2location_v2_frozen_20260916__scRNA32`，字段 `reference_celltype_v2` 的32类词表；参考文件名含30type不改变实际词表。三套已训练CellTypist模型复用，没有重训全参考。两张输入features顺序与映射已验证：14,329个映射features合并为14,327参考基因。按这些实测counts做10k归一化、log1p与训练时相同scaler及截断；没有把scVI latent或估计表达传给分类器。

冻结留出参考4745细胞复用于深度校准，按每张新细胞共享基因的真实中位counts稀释：Cervilla普通/低信息={'ordinary': 355, 'low_information': 21}；D1={'ordinary': 307, 'low_information': 20}。类别至少20个接受预测、经验精度≥0.8才获得可用margin门槛，结合full-depth是否可校准；这不是独立HD真值准确率。保存候选、sigmoid模型支持分数、margin及校准接受状态；支持分数不等同于经HD真值校准的身份概率。

全局表达簇检查主要谱系；在预测谱系内重建表达图和分群，再按切片分别审阅实测marker频率。模型候选仍逐细胞保留。细胞自己的marker组合或相对一致表达群的成组marker均可提供证据，不再以每个细胞检出两个低表达特异marker作为统一否决。群体证据只支持模型自身候选一致的细胞，不将标签传播给簇中所有细胞；详见 `metadata/annotation_rules.json` 与每张 `cluster_annotation_review.tsv`、`cluster_marker_evidence.tsv`。参考和cluster证据不是独立真值。

主要谱系计数（包括保留的mixed对象原有主要谱系建议；mixed身份应结合状态排除于确定类型统计）：

| major | Cervilla | D1_Johnson |
| --- | --- | --- |
| B | 4 | 2956 |
| Endothelial | 1123 | 2230 |
| Epithelial | 55136 | 75850 |
| Fibroblast | 8332 | 32507 |
| Mast | 56 | 208 |
| Myeloid | 1524 | 12601 |
| Plasma | 3449 | 6230 |
| T/NK | 93 | 10360 |
| Unknown | 28210 | 67680 |

最终注释状态，分母为该张所有QC保留细胞：

| status | Cervilla | D1_Johnson |
| --- | --- | --- |
| lineage_unresolved | 28129 | 67513 |
| mixed_or_ambiguous | 3323 | 8914 |
| reference_broad_identity | 55484 | 73944 |
| subtype_supported | 1822 | 11635 |
| subtype_unresolved | 9169 | 48616 |

获得亚群级支持的建议标签：

| subtype | Cervilla | D1_Johnson |
| --- | --- | --- |
| Arterial EC | 11 | 323 |
| C1QC TAM | 40 | 0 |
| CD8 GZMK/Effector-memory-like | 0 | 17 |
| Capillary EC | 0 | 2 |
| Classical Mono | 1 | 9 |
| FOLR2 TAM | 3 | 19 |
| LAMP3 DC | 1 | 64 |
| Lymphatic EC | 10 | 1 |
| Matrix_CAF | 1517 | 10469 |
| Memory B | 0 | 20 |
| NK | 0 | 16 |
| Naive B | 0 | 2 |
| Neutrophil | 2 | 0 |
| Pericyte | 53 | 340 |
| SMC | 123 | 172 |
| SPP1 TAM | 2 | 16 |
| Tfh/CXCL13-like | 0 | 25 |
| Treg | 0 | 13 |
| Venous EC | 2 | 10 |
| cDC1 | 2 | 2 |
| cDC2 | 1 | 13 |
| iCAF | 51 | 74 |
| pDC | 2 | 20 |
| rCAF | 1 | 8 |

`predicted_subtype`为模型候选，`final_label`为建议标签，`annotation_status/evidence`保留支持路径及不确定原因。主要谱系明确但亚群不足时保留对应父类unresolved。Epithelial、Plasma、Mast在本参考粒度下保留broad identity；上皮不推断恶性状态。GC B需B谱系及成组GC证据，不能仅凭BCL6/cycling。`low_information`与身份不确定分开；混合嫌疑来自实际跨谱系marker共存，不等同于已证实doublet。

实际重点谱系/亚群审阅结论：Matrix_CAF在Cervilla三个谱系内簇、D1六个簇有成组COL11A1/COL10A1/MMP11/POSTN/CTHRC1/FAP（D1亦LRRC15）支持，最终分别1517/10469。D1 Arterial EC有GJA5/EFNB2/HEY1簇证据（最终323；Cervilla仅11个逐细胞建议）。上皮保持broad；Plasma及多数免疫谱系仍有大量候选未获marker支持。Cervilla B候选仅约2.0%检出任一B marker，T/NK候选约3.3%，不能将模型候选直接当作细胞组成，也不能将未支持解释为生物学不存在。SPP1 TAM仅2/16个逐细胞组合支持，没有形成满足本轮一致性条件的SPP1群体；GC B两张均0，模型候选190/1087仍保留。其他极小亚群的逐细胞建议不等于已验证稳定群体。COL1A1及部分基质信号横跨多个标签，可能含环境/跨细胞污染或重建边界影响，不据此宣称所有类型均产生胶原。

低置信结果按表达簇与谱系集中审阅，16行谱系概览见 `metadata/uncertain_lineage_review.tsv`，每张一行一个表达群的概览见 `metadata/uncertain_group_review.tsv`；详细证据见每张 `cluster_annotation_review.tsv`、`lineage_cluster_top_markers.tsv` 与marker dotplot；没有创建逐个Unknown人工任务。描述性Welch排序只用于浏览实测log表达差异，不是独立确认性检验。

## 两张共同表达空间

scVI实际输入 **2 sections**：{'Cervilla': 81789, 'D1_Johnson': 170812}，合计252,601核心细胞、2000个HVG；section_id为batch，NB，latent20、hidden128、2层，seed123，batch1024，最多50epochs、KL预热10epochs、early stopping patience10，实际50epochs。GPU为Quadro RTX 6000，PyTorch 2.8.0+cu128 / CUDA 12.8。训练约2.5分钟。HVG按两张各自变异排名等权选择并优先共同HVG；每张保留全基因主对象。

有HVG计数的保留细胞保存latent；零HVG计数的Cervilla 1,948 / D1 4,330个细胞latent明确NA，均保留在全基因主对象及空间分母。初始核心中676个零HVG对象在输入复核后转为仅预测，初次短训练被替代；最终模型按修正核心训练。表达近邻、Leiden、UMAP仅用核心，非核心UMAP明确NA但空间坐标和总数分母保留。另保存未整合log-count PCA基线。最终模型验证ELBO从112.328降至93.319（最小值约93.312），最后5轮范围约0.030，KL权重为1；达到50轮上限，未再次因未成熟的预热早停。这里的表达近邻不是物理邻域，未建立跨切片物理邻接。每来源仅一张，section、队列与真实组织差异混杂，不能把充分混合作为验收标准。没有为UMAP外观反复训练。一次默认400epochs的KL预热与50epochs上限不匹配，导致早停时KL权重仅0.025；查看实际ELBO与重建损失后，将预热设为10epochs重新训练，相关轨迹与原因保留于scvi_warmup_correction.json和日志。

## Cervilla v1/v2共同覆盖比较

| version | reconstructed_full | retained_full | low_information_full | common_reconstructed_fully_covered | common_retained_fully_covered | common_assigned_fraction | common_median_counts | common_median_genes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| v1 | 106722 | 97847 | 12798 | 106722 | 97847 | 0.709007 | 271 | 227 |
| v2 | 107071 | 97927 | 12573 | 107071 | 97927 | 0.705695 | 268 | 224 |

共同覆盖按同一原始bin表两版image_covered交集定义；counts只统计交集bins，细胞比较保留所有已分配bins均在交集的对象，不按cell_id或行号配对。v1旧报告重建数106,934含212个零表达对象；本表统一按非零表达计算106,722。完整组成/状态对照见 `metadata/Cervilla_v1_v2_composition.tsv`。D1旧结果为NA，不记作0。

本轮图像、模型图像归一化视图和注释证据使用方式均有变化，不能把全部标签变化归因于图像恢复。分割/计数变化与注释决策分别列出；不以细胞数增加、Unknown下降或genes中位数提高单项判断成功。相同物理ROI的分割对照见 `Cervilla_v1_v2_segmentation.png`。实际比较审查：相同物理ROI相交核数：dense 1707→1709，stroma 1260→1241，edge607→612；局部拆分/边界变化存在，但没有可据以声明全片分割准确率大幅改善的真值证据。共同覆盖counts分配70.90%→70.57%，保留细胞97847→97927，中位counts271→268、genes227→224。Unknown下降包含群体marker证据规则改变，不能都归因于新图像。

## 交付与验证

主对象：
- [Cervilla/cells.h5ad](</home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/cell_annotation_v2/Cervilla/cells.h5ad>)
- [D1_Johnson/cells.h5ad](</home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/cell_annotation_v2/D1_Johnson/cells.h5ad>)

每张同时保存 `cell_labels.parquet`、`excluded_cell_metadata.parquet`、压缩核mask、2µm array cell mask、bin_to_cell映射、未获bins核几何及实际变换。mask与bin-to-cell保存重建时全部数值标签（包含QC排除对象），通过source_cell_id连接最终保留表，不把所有mask标签都当作保留细胞。没有将原图或重复全基因矩阵嵌入对象。最终共享模型为 `models/scvi_two_sections/`；参考模型只读复用v1。GPU同版本软件覆盖层复用原环境其余依赖，没有修改v1环境。

每张 `figures/` 包含原图/核/扩张ROI、QC、major/subtype/不确定状态UMAP、全捕获空间图、重点谱系图、marker dotplot和代表大核。跨切片同一配色、透明背景；绘图分母见标题。共享表达基线和新旧比较位于 `metadata/`。

验证通过102项，重点为原始整数X、完整features、行counts/genes、映射/标签对齐、原始/局部/模型/µm坐标及低信息分母；每张随机25个细胞加最高counts细胞独立从原2µm H5重汇总比对全部基因，非重复全量复算。结果见 `metadata/delivery_validation.tsv`。图像与marker是计算审查，不替代病理专家或独立细胞真值验证。

## 补充数据支线

复用既有65行候选审计，更新于 `metadata/candidate_datasets_v2.tsv`。本轮没有新增就绪HD切片、没有新增平台全流程：
- GSE326226另外四张Xenium5K仅各读131,072字节H5 features，实际均5001个gene targets、固定Top20覆盖9/20；与此前A2R2一起五张现均有实际面板核对。只登记覆盖，不替补、不计算MP10分数。相同研究/跨平台不重复计作独立患者。
- UNC LCCC_553_Xenium读取嵌套ZIP目录与少量Zarr元数据（cells约0.78MB、transcripts约3.17MB）：作者细胞数258,581，cell/nucleus centroid、area及边界schema、µm单位已确认；转录本gene/quality/location字段已确认。未找到独立原始cell count矩阵，转录本到细胞分配尚未建立，不能称为已就绪细胞counts。未下载9.1GB整包。
- MPBC1/2仍需作者细胞原始counts+空间映射，或同片图像/positions/scales；UC_SCC仍缺现成2µm及形态输入，未下载22.69GB FASTQ。STT0000103与Oh/Jung仍属受限访问。其他已核候选沿用前轮证据，没有把未复查项写成新近复验。

人工操作集中于 [MANUAL_INPUT_ACCESS.md](</home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/cell_annotation_v2/metadata/MANUAL_INPUT_ACCESS.md>)，列明完整链接、所需文件、阻碍、操作、落盘位置与可解锁切片。未联系作者/发邮件/索取凭据，未让这些缺项阻塞本次两张。

## 后续是否具备物理邻域分析基础

已有两张独立cell_id、真实原始counts、经检查物理坐标、分割映射、建议谱系/亚群和明确不确定性字段，具备从本审阅版设计后续分析的技术基础。血液/碎片区和广泛基质背景仍需针对性审查，不能把本版等同于已证明准确的单细胞真值。应按切片分别构建物理图，预先约定低信息、未定类及mixed的处理，对分割/扩张和标签不确定性做敏感性检查；稀少或不支持的亚群不可强行解释，两个来源各一张不足以分离患者/批次/组织效应。

本轮没有运行MP10/胆固醇评分、neighborhood/niche、spottedpy、空间阳性检验、通讯或CNV。交付后停止，不自动启动下游。
