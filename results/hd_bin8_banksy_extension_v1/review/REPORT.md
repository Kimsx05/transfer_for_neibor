# 8 μm BANKSY extension: status, evidence and limits

Updated: 2026-11-11T10:24:54+08:00 (Asia/Shanghai)

All executable branches are complete and reread-validated. Original two-section results are reused; three new sections have full-bin results. Missing tissue/image inputs remain explicit limitations.

## 本轮结论

本轮可执行分支均已完成：D1补图、UC_SCC原生8μm分析、MPBC1/2相对网格恢复和8μm聚合、三张新片的全量逐bin两参考SingleR及BANKSY三路、原始表达审阅、五片统一汇总和加载验证。缺失的MPBC组织掩膜/H&E与UC原始高分辨率H&E没有被伪造或当作已恢复。

|切片|分析bins|原始UMI|中位nFeature|nFeature<20|建议区域主分支|
|---|---:|---:|---:|---:|---|
|Cervilla|268,924|64,443,819|114|18.10%|既有λ0.2；λ0.8作粗结构对照|
|D1_Johnson|364,805|131,440,767|179|17.93%|既有λ0.8区域，辅以λ0.2；增加角蛋白/胶原审阅层|
|UC_SCC|692,681|126,389,649|127|4.71%|λ0.2；λ0.8仅作粗区域解释|
|MPBC1|675,380|107,299,054|124|9.08%|λ0.2，保留小型免疫及血管相关成分|
|MPBC2|422,741|39,606,795|35|37.67%|λ0.2保守探索；λ0.8不能填补低信息身份|

前三片分母为supplied、组织内、nCount>5；MPBC仅为supplied且nCount>5，组织范围未知。全部保留18,085输入基因。MPBC聚合前后逐基因及全局守恒，QC去除UMI另有账目。旧片没有重跑SingleR/BANKSY/RCTD。

**D1解释有所改变，但没有改写正式标签。** λ0.8簇1/2的KRT14、KRT6A、KRT17、KRT5检出约77%–85%，并且相对其他区域富集；COL1A1/2仍广泛检出。更合适的补充描述是“上皮信号明显，同时伴胶原信号”。D1正向参考冲突37,972 bins中86.97%有≥2枚扩展角蛋白、67.77%同时有≥2枚胶原基因。这不支持将整块简单当成纯成纤维区，但也不能证明真实混合细胞比例、EMT或恶性。原规则区域名、两参考候选均保留。

**UC_SCC有清楚原始结构，但没有解决参考冲突。** 两路均有候选的四粗类分歧为Cervilla16.26%、D1 28.81%、UC27.66%、MPBC1 17.88%、MPBC2 25.78%；相应UMI加权为6.78%、26.38%、25.50%、10.10%、12.96%。同一信息量分层加权的全部分析bins分歧为19.51%、30.36%、30.19%（C/D/UC），只比较这三个有原组织范围的切片。因此UC在相似检测量下并未显示比D1更稳定的参考标签。加权指标与“两参考均有候选”原始比例分母不同，不能直接相减。

UC的KRT5/6A/6B/14/17、COL17A1空间巢与H&E粗结构相符；typical/high ROI在上皮样巢及纤维样区之间呈现角蛋白/胶原分布差别，局部有共同信号。其鳞状分化背景与D1的扩展角蛋白模式，比Cervilla/MPBC的EPCAM/KRT7/KRT19/UPK2模式更接近。单一“尿路上皮marker必须检出”规则会误伤前者，单一“扩展角蛋白必须检出”规则又会误伤后者。

**有直接证据支持参考适用性/特征选择参与分歧，同时存在样本与输入差异，不能简化成旧片质量差。** 六个新增上皮基因在两路query输入均可匹配，但冻结BLCA模型只选择了KRT17作为判别marker，HPCA选用了全部六个；两路都选用了四个胶原基因。冻结BLCA上皮参考中KRT14/6A/5/COL17A1检出9.2%/13.0%/19.0%/10.3%。这说明有效参考特征覆盖与D1/UC的表型不同；不是ID错配，也不是分类器准确率验证。没有做特征消融，不能定量声称这一因素解释了全部冲突。胶原的具体来源（相邻细胞贡献、组织复合、其他生物/技术来源）仍未完全区分。

MPBC1的原始尿路上皮对照基因和胶原表现出较清楚的上/下区域差异，正反向参考冲突各约1.9% bins；这种不同在描述层面支持组织构成/输入数据差别，却不能排除缺掩膜的分母影响。MPBC2大量低信息位置与输入缺口、部分条带重合，但局部高信息上皮/间质/免疫相关区域仍有原始表达支持。不能把其2μm输入统计与8μm结果混比，也不能将缺少组织掩膜造成的未知背景归属问题说成已解决。

**空间信息的增益与限制。** 三片λ0即有部分表达结构；λ0.2通常减少散点并保留区域内表达差别，λ0.8形成更连续的较大区域。UC免疫/间质复合簇16、MPBC1髓系簇16和血管相关簇7、MPBC2免疫相关簇4有成组原始表达支持，不能仅按最大marker面板将它们看成纯间质。MPBC1 λ0.8没有独立保留较小的免疫解释区；MPBC2低信息ROI在λ0.8高度统一却没有新增原始marker证据。连续性不是准确率；本轮没有据此减少Unknown或做邻居投票。

两套候选标签均继续保留，不选一个参考作为全片唯一真值。D1/UC扩展角蛋白明确的区域可参考HPCA上皮候选，但需保留BLCA分歧及胶原证据；MPBC与Cervilla应继续用其适用的尿路上皮面板核验。`ext_recommended_region`是独立区域审阅层，`ext_region_review_note`记录少量明确修订的证据，原类型层/区域层未覆盖。

下一轮可在Cervilla/D1/UC和MPBC1原始表达支持的区域做预先定义的区域层面比较；MPBC2限定信息充分局部并保留低信息对照。正式细胞类型特异、精细亚群或病理亚型归因仍需参考适用性复核、MPBC组织范围/H&E补齐及局部共表达来源核验。UC/MPBC不是与旧片病理构成匹配的技术重复；本轮不启动MP10、通讯、恶性判定或邻域生态。

未恢复事项仅为输入与解释限制：MPBC组织掩膜/图像配准、UC原始高分辨率H&E，以及胶原共表达的细胞来源。计算分支没有以试跑或sketch代替全片；所有新片主结果都是完整有效bins的直接计算。

现有本地pipeline摘要另显示Cervilla约477.4百万reads、81.3% saturation，D1约221.5百万reads、30.05% saturation；confident probe mapping分别90.2%与97.19%。这些是各自原始文库/覆盖范围的指标，不是本轮相同组织面积的深度；新三片没有定位到可比的测序/比对摘要，不能从表达矩阵反推。因此目前既不能将所有差异归因于测序深度，也不能据这些数字对五片作统一质量排名。源路径与原始单位保存在`available_pipeline_metrics.tsv`。


## 完成状态

| section | input_complete | normalization_complete | annotation_complete | lambda0_complete | lambda02_complete | lambda08_complete | validation_complete | visual_review_complete | tissue_mask |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | True | True | True | True | True | True | True | True | original |
| D1_Johnson | True | True | True | True | True | True | True | True | original |
| UC_SCC | True | True | True | True | True | True | True | True | original |
| MPBC1 | True | True | True | True | True | True | True | True | unknown |
| MPBC2 | True | True | True | True | True | True | True | True | unknown |

## Input and QC

| section | input_bins | input_resolution_um | bins8 | input_UMI | analysis_bins | analysis_UMI | excluded_UMI | tissue_mask | denominator |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | nan | 8 | nan | nan | 268924 | 64443819 | nan | True | supplied in-tissue nCount>5 |
| D1_Johnson | nan | 8 | nan | nan | 364805 | 131440767 | nan | True | supplied in-tissue nCount>5 |
| UC_SCC | 695646.0 | 8 | 695646.0 | 126400557.0 | 692681 | 126389649 | 10908.0 | True | supplied in-tissue nCount>5 |
| MPBC1 | 11222500.0 | 2 | 702244.0 | 107385169.0 | 675380 | 107299054 | 86115.0 | False | supplied nCount>5; tissue unknown |
| MPBC2 | 9781598.0 | 2 | 615386.0 | 40043738.0 | 422741 | 39606795 | 436943.0 | False | supplied nCount>5; tissue unknown |

MPBC uses supplied-expression QC only; original tissue membership remains NA. It cannot be pooled as a tissue-verified denominator with Cervilla/D1/UC_SCC. Original matrices and old official results are preserved.

## D1 supplementary evidence

| group | symbol | pct_positive | mean_log | log2_mean_effect |
| --- | --- | --- | --- | --- |
| 2 | KRT14 | 78.95% | 2.986 | 1.006 |
| 2 | KRT6A | 77.49% | 2.946 | 1.153 |
| 2 | KRT6B | 46.08% | 1.481 | 1.564 |
| 2 | KRT17 | 83.11% | 3.398 | 1.488 |
| 2 | KRT5 | 82.40% | 3.298 | 1.28 |
| 2 | COL17A1 | 44.51% | 1.349 | 0.742 |
| 2 | EPCAM | 3.58% | 0.092 | -0.09 |
| 2 | KRT7 | 0.46% | 0.013 | -1.294 |
| 2 | KRT19 | 20.92% | 0.574 | 0.189 |
| 2 | UPK2 | 0.02% | 0.001 | -0.12 |
| 2 | COL1A1 | 89.35% | 3.909 | 0.533 |
| 2 | COL1A2 | 75.11% | 2.956 | 0.548 |
| 2 | DCN | 43.45% | 1.468 | -0.086 |
| 2 | LUM | 39.11% | 1.298 | 0.098 |
| 1 | KRT14 | 85.32% | 3.37 | 1.44 |
| 1 | KRT6A | 82.17% | 3.177 | 1.392 |
| 1 | KRT6B | 38.34% | 1.125 | 0.859 |
| 1 | KRT17 | 79.42% | 2.931 | 0.728 |
| 1 | KRT5 | 83.50% | 3.249 | 1.13 |
| 1 | COL17A1 | 60.25% | 1.872 | 1.418 |
| 1 | EPCAM | 4.43% | 0.105 | -0.085 |
| 1 | KRT7 | 0.81% | 0.022 | -1.244 |
| 1 | KRT19 | 11.90% | 0.297 | -0.668 |
| 1 | UPK2 | 0.01% | 0.0 | -0.126 |
| 1 | COL1A1 | 84.75% | 3.481 | -0.001 |
| 1 | COL1A2 | 71.10% | 2.646 | 0.146 |
| 1 | DCN | 47.80% | 1.562 | -0.012 |
| 1 | LUM | 45.41% | 1.457 | 0.252 |

These raw-expression measurements supplement the old epithelial/collagen conflict interpretation. Keratin detection and relative enrichment are separated; collagen detection does not by itself prove fibroblast identity, cell admixture, EMT or malignancy. Old region labels remain unchanged.

## Cervilla

实际打开本轮统一14基因面板的全片原始表达图，并沿用上一轮已复核的固定ROI/H&E和域解释。Cervilla 的 KRT14/KRT6A/KRT6B/KRT17/KRT5/COL17A1 扩展组整体稀疏，而 EPCAM/KRT7/KRT19/UPK2 对照组有较明确的区域信号；胶原/间质组呈穿插分布。不能把针对D1补齐的角蛋白组直接改成所有切片的统一上皮判定门槛，否则会引入新的组织构成偏差。

本轮全分析集合中，扩展角蛋白>=2检出约1.03%，胶原组>=2约17.72%，两者共检出约0.31%；这些是指定窄面板的检出，不是上皮/间质细胞比例。旧的完整marker面板和两参考候选保留，区域仍建议lambda0.2为主、lambda0.8作粗区域对照。没有supplied bins的H&E区域不参与统计，图中无点不等于实测零表达或Unknown。

新增字段已经按ID挂回旧对象并实际重读验证；没有重跑旧片模型或修改旧区域标签。结合已完成的跨片比较，不能据此判定Cervilla质量较差。


## D1_Johnson

已实际打开本轮全片角蛋白图、conflict_focus 单基因对照图，并结合上一轮已复核的同坐标 H&E/域图。14个指定基因全部存在，没有以近似基因替换。

补图让主要冲突区域的上皮证据更明确：lambda0.8簇1/2的 KRT14、KRT6A、KRT17、KRT5 检出率分别约79–85%、77–82%、79–83%、82–84%，相对其他 bins 的平均线性标准化表达均为正效应。COL17A1也有支持。相反，EPCAM约3.6–4.4%、KRT7约0.46–0.81%、UPK2<0.03%；只看后一组会低估这里的上皮相关表达。

conflict_focus 中 KRT14/KRT5/KRT6A/KRT17 呈连续而非单个孤立点信号，COL1A1/COL1A2 也广泛检出。需要区分“检出”和“富集”：簇1 COL1A1检出84.75%，但相对其余 bins 的 log2(1+均值)效应约−0.001，几乎没有区域特异富集；簇2 COL1A1为89.35%、效应约0.533。不能把胶原阳性直接等同为成纤维细胞身份。

建议保留原区域标签，新增补充描述“上皮信号明显，同时伴胶原信号”。这更支持原参考冲突存在面板/参考响应差异，不能据此消除实际共检出的胶原，也不能区分 bin 内多细胞贡献与细胞自身 ECM 表达；不推断 EMT、双细胞比例、CAF或恶性。双方一致上皮候选区域作为比较，不视为纯度真值。

原来的两套类型候选均保留；lambda0.8组织域配合lambda0.2局部结构仍适合继续核验。是否可进行细胞类型特异归因仍未解决。此次没有重跑旧片模型、逐binSingleR或RCTD。

冻结模型核验提供了具体参考机制线索：六个扩展上皮相关基因在BLCA test.genes中全部可匹配，但实际选用的732个判别marker中仅KRT17在内，KRT14/KRT6A/KRT6B/KRT5/COL17A1均不在；HPCA则选用了全部六个。四个胶原/间质基因被两模型共同选用。这是有效判别特征选择差异，不是query缺基因或Ensembl未对齐。它可以解释部分角蛋白明确区域仍被BLCA偏向间质参考的现象，但未进行因果消融，不能量化各因素贡献或据此把HPCA全部判为真值。完整检查在`tables/reference_expanded_panel.tsv`；冻结BLCA参考的对应基因表达分布另存，不重新训练。

The frozen BLCA epithelial reference itself detects KRT14/KRT6A/KRT5/COL17A1 in only9.2%/13.0%/19.0%/10.3% of its2,807 epithelial cells, whereas D1 conflict clusters detect these much more broadly. This describes reference/query phenotype mismatch; per-cell and per-bin detection fractions are different units and cannot be treated as an accuracy or direct purity comparison. The same reference also has nonzero keratin detection in mesenchymal classes (e.g. KRT17:31.0% fibroblast,26.0% mural); this warrants caution about marker exclusivity, not relabelling the frozen reference here.


### UC_SCC input review

实际查看已有hires图的1600px检查缩略图（仅为显示，不用于推断mpp）以及原生8um输入/QC概览。下采样H&E可辨认大块上皮样巢、纤维间隔、较深染细胞密集区域和裂隙，支持后续粗结构对照；不据此做核分割或微细细胞身份判定。物理坐标来自原始scalefactors，非图像大小估算。692,681个bins通过原组织内+nCount>5条件，中位nFeature127；该片并非简单拥有高于D1的总体检测量。这是建模前的输入观察；最终SingleR、BANKSY及跨片结论见下文完成后的审阅。

The completed, unsmoothed per-bin coarse candidate map was actually opened before BANKSY domain interpretation. Both references already show large epithelial-like islands and interleaved mesenchymal bands. HPCA calls broader epithelial areas; BLCA adds mesenchymal/immune assignments within parts of these islands, especially the central/lower and right-central capture area. Thus UC_SCC has organized candidate structure before spatial clustering, but reference disagreement is not confined to isolated low-signal pixels. Candidate structure alone is not original-marker or histological validation; the completed raw-panel/domain review is recorded separately in visual_review.md.


## UC_SCC

已完成全量692,681 bins、126,389,649原始UMI的两路逐bin SingleR与λ=0/0.2/0.8（45/24/40簇），实际重读验证通过。所有14个补充基因均可匹配；并非EPCAM/UPK2 ID缺失造成低检出。

实际查看了全片候选图、三个空间分支、原始三组marker、high/typical ROI单基因与cluster图、固定四ROI的下采样H&E，以及整片H&E配准叠图。typical ROI (3584–3840,2048–2304 μm)和high ROI (6912–7168,1536–1792 μm)中，左侧/斜行纤维样区对应COL1A1/COL1A2，右侧上皮样巢对应KRT5/KRT6A/KRT6B/KRT14/KRT17/COL17A1；有局部共检出，不能据此推断单个细胞混合或EMT。EPCAM/KRT7/KRT19/UPK2在这些ROI很稀疏，与鳞状分化背景相容，但不在本轮作恶性或亚型判定。low-information ROI仍有可见组织结构，低信息不等于无组织；sparse_input_edge是捕获范围边缘，图像可见不等于表达覆盖。

λ=0已经具有组织性；λ=0.2更清楚分开角蛋白巢、胶原带及部分免疫相关小区，适合作为主区域分支。λ=0.8适合较粗区域浏览，不能把被连成大片的间质/复合区域解释成均一细胞类型。相邻同簇率0.583→0.882→0.962仅描述连续性。尤其λ0.2簇16同时具有JCHAIN/MZB1 (26.7%/29.1%检出，表达效应2.665/2.279)与COL1A1/COL1A2 (84.3%/65.7%)；原规则“Mesenchymal-enriched”遗漏免疫成分，本轮另加“Immune-stromal composite”审阅层，不覆盖原字段。簇17有MS4A1/CD79A/TRAC成组支持，簇22有VWF/PECAM1/ENG及胶原/收缩相关表达，需保留局部复合解释。

64,352个BLCA Fibroblast / HPCA Epithelial bins占全部QC的9.29%、UMI的11.57%；其中69.47%检出至少2枚扩展角蛋白，35.59%同时检出至少2枚胶原面板基因。冲突不能只归因于低深度。两参考均有候选时四粗类分歧27.66%，与D1接近；按相同信息量加权后也无明显优势。保留两路候选；扩展角蛋白和H&E支持的巢区可作为上皮相关区域继续分析，但不能把HPCA候选普遍升级为细胞身份真值。


### MPBC1 input review

实际查看聚合8um全supplied网格的nCount/nFeature图。输入几乎覆盖一个矩形阵列；检测量本身存在连续结构、中央低量带及细条带，但没有H&E或原组织掩膜，不能把这些形状直接称为解剖结构。675,380个bins通过表达QC，中位nFeature124。全阵列提供的零表达2um位置与未提供位置已分开；floor(row2/4)、floor(col2/4)固定原点聚合通过逐基因守恒，并用R加载原始稀疏8um矩阵核对成功。后续domain需要检查是否主要追随检测量或条带。这里记录建模前输入观察，最终参考和域结论见下文。

The full, unsmoothed per-bin coarse candidate map was actually opened before BANKSY completion. Both references show broadly corresponding lower-half epithelial candidate territories, upper-half interleaved epithelial/mesenchymal patterns, and a central oblique immune/mesenchymal-rich band. Fine labels are not thereby validated; merging Fibroblast/Mural into Mesenchymal removes a reference-coverage distinction. White holes in this QC-only map are absent/excluded analysis positions, not confirmed anatomical lumina. No H&E is available to independently verify the apparent structures or their relation to micropapillary versus urothelial histology.


## MPBC1

已完成实际supplied2μm整数counts的固定原点8μm聚合（逐基因和全局守恒），675,380个表达QC bins、107,299,054UMI，两路逐bin SingleR、三路BANKSY（114/29/32簇）及实际加载验证。in_tissue_original始终NA。没有H&E、原始组织掩膜和核边界，不能进行病理结构或细胞身份的直接确认。

实际查看输入信息图、两参考四粗类图、三路区域图、全片原始marker、high/typical/sparse-edge ROI单基因和typical ROI聚类。下半部及部分弯曲带具有EPCAM/KRT7/KRT19/UPK2信号，上半部COL1A1/COL1A2/DCN/LUM较强；扩展KRT14/6A/6B/17/5/COL17A1面板整体很弱，不能因它们阴性把上皮候选否定。typical ROI (1536–1792,5376–5632 μm)表现尿路上皮对照基因与较少胶原，high ROI局部兼有COL1A1/2。sparse_input_edge位于最右侧仅部分supplied的窗口；未提供的部分保持空白，不是实测零或已知组织外。

λ=0已有上/下部表达结构。推荐λ=0.2：保留区域内部差异和一个1,371 bins的髓系相关簇16（APOE/LYZ/CD68/CTSS/C1QA/B/C），其CD68检出37.1%、相对表达效应3.811；λ=0.8把上半部合并成更连续的间质区域，局部小型免疫相关表达不再作为独立解释区，不宜仅凭更平滑优先选择。λ0.2簇7有RGS5、PLVAP、COL4A1/2、ENG、VWF、PECAM1成组相对富集，分别RGS5/VWF/PECAM1检出29.5%/5.5%/7.3%，原规则为Unresolved；另加“Vascular-associated composite”审阅层，不能将整个簇认作纯内皮或纯Mural。

两路均有候选四粗类分歧17.88%，UMI加权10.10%。正反上皮/成纤维冲突各约1.9% QC bins，明显不同于D1/UC的单向大块冲突。然而缺组织筛选、组织构成及测序来源不同，这不是更高质量或更准确的证明。保留BLCA/HPCA独立候选，λ0.2用于相对网格区域探索；MPBC病例为额外MIBC混合微乳头/尿路上皮标本，不能由本文题目误归NMIBC，也不能将本图区域直接命名为微乳头病理区。


### MPBC2 input review

实际查看8um聚合的全supplied网格nCount/nFeature图。存在大块未提供位置、多个分离的输入支持块、局部高信息区，以及相对阵列y约6300–6700um的明显低检测量横带。低量带不能在没有组织/图像证据时直接删除或命名为空白背景。表达QC保留422,741个bins，中位nFeature35，约37.7%落在nFeature<20；与前三片的组织内分母不同。仅凭全片均值会掩盖局部高信息区域，后续采用既定分层和同尺度ROI。

输入连通块将阻止跨独立supplied区域连边；这不证明背景已排除。无组织掩膜、配准图像或微乳头/尿路上皮病理分区，不能按相对坐标自行复现作者病理标注。聚合计数守恒及R稀疏counts加载已通过。这里记录建模前输入观察；两路SingleR与三路BANKSY现已完成，结论见下文。


## MPBC2

已完成2μm→8μm固定原点原始counts聚合、逐基因及全局守恒、422,741个表达QC bins（39,606,795UMI）的全量两参考SingleR、BANKSY三路（253/40/41簇）、空间绘图和实际加载验证。组织掩膜/H&E仍缺失；in_tissue_original=NA，没有用表达阈值或预测标签伪造组织掩膜。

实际查看全阵列输入信息、三路区域、全片原始面板、high/typical/sparse-edge单基因ROI及low-information ROI三路clusters。高信号局部存在KRT7/KRT19/EPCAM/UPK2与胶原分布；扩展角蛋白整体极弱。typical窗口处于低信息的上部阵列区域，14基因大多未检出；sparse-edge只有零散输入点，不能把空白解释成组织腔隙。low ROI在λ0.8变成连续同簇，但原始信号仍弱，属于空间连续性增强而非身份得到证实。

37.67%的QC bins nFeature<20，但只承载4.62% UMI；中位nFeature35、UMI37。因此不能仅看全部bin比例否定局部信息充足区域，也不能只看UMI加权指标忽视广泛低信息位置。上部、部分中心带和输入缺口附近的区域解释要谨慎。38个输入支持连通块中，258 bins邻居少于50、2个singleton仅采用自身kernel；无跨连通块边。实际第50/最大邻居距离99%分位80μm、最大299.87μm，长边位置已保存。它们是输入几何限制，未证实为独立组织块。

推荐λ0.2作为保守探索主分支，λ0.8仅用于粗结构对照。λ0.2簇3/8及若干局部有原始上皮相关表达；簇4有MS4A1/CD79A/TRAC成组相对富集（检出13.0%/8.2%/17.0%），支持免疫相关区域而非逐bin均为免疫细胞。λ0.8相邻同簇率97.2%不能证明背景被排除。原始counts支持的局部区域可继续探索；低信息大片、缺口边缘及未核验组织范围不适合细胞类型特异归因。两参考均有候选四粗类分歧25.78%，需同时报告37.67%低信息和约10%未双路可评价，不能与组织内分母简单排名。


## Same-denominator reference comparison

| section | level | n | disagree_fraction | disagree_UMI_fraction |
| --- | --- | --- | --- | --- |
| Cervilla | original_major | 260214 | 27.13% | 10.74% |
| Cervilla | four_coarse | 260214 | 16.26% | 6.78% |
| D1_Johnson | original_major | 349659 | 40.26% | 29.25% |
| D1_Johnson | four_coarse | 349659 | 28.81% | 26.38% |
| UC_SCC | original_major | 690288 | 39.82% | 32.15% |
| UC_SCC | four_coarse | 690288 | 27.66% | 25.50% |
| MPBC1 | original_major | 662017 | 29.10% | 17.22% |
| MPBC1 | four_coarse | 662017 | 17.88% | 10.10% |
| MPBC2 | original_major | 379880 | 43.25% | 21.62% |
| MPBC2 | four_coarse | 379880 | 25.78% | 12.96% |

Per-information-band bin/UMI comparisons are in `tables/reference_disagreement_strata.tsv`; expanded raw-panel support is in `tables/information_strata.tsv` and `tables/reference_conflict_strata.tsv`. `information_standardized_sensitivity.tsv` gives a single bounded common-information weighted comparison; no models/labels are changed. Histological composition and assay/input differences prevent ranking sample quality from these percentages alone.

## Cross-section interpretation

见报告开头的本轮结论；分母和分层原表随交付保留。

## Reading and rerun

```r
root <- "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_bin8_banksy_extension_v1"
source(file.path(root,"code/load_result.R"))
x <- load_extension(root,"UC_SCC") # or MPBC1, MPBC2, Cervilla, D1_Johnson
```

The loader attaches the completed incremental results by unique bin ID. Final object validation is recorded separately in each `validation.json`. New sections reference native UC_SCC H5 or their sole `counts8.h5`; `normalized.rds` is the sole normalized layer. Old sections depend on sibling `hd_bin8_comparison_v1` and `hd_bin8_banksy_v1`.

`bash rerun.sh` reuses completed models and runs light postprocessing directly. If model stages are missing, it uses the night supervisor (23:00–07:00 CST), processes UC_SCC then MPBC1/2, resumes caches, and updates comparisons/report. Active heavy child groups pause outside that window. Inspect `logs/night_status.json`, `logs/night_supervisor.log` and `logs/run_SECTION.log`. No stage-complete claim is made merely because a job is queued.

Required interpretation tables: `cluster_domain_interpretation.tsv`, `cluster_summary.tsv`, `cluster_markers.tsv.gz`, `panel_marker_statistics.tsv.gz`, `cluster_SingleR.tsv`; predictions retain both references, pruning and scoring. Supplementary 14-gene evidence is in `supplement_marker_statistics.tsv` and `supplement_gene_coverage.tsv`.

All figures are white300dpi PNG plus PDF/legend/panel variants where appropriate. Review archive selects key PNGs and evidence tables; all complete figure formats remain on server. No expression/reference/large embedding data are in review.zip.

This round does not perform MP10, subtyping, CNV/malignancy, communication, formal neighborhood ecology, new segmentation or RCTD.


Incremental files: 3.633 GiB logical, 3.635 GiB allocated, including review.zip (96.3 MiB). This value includes current temporary caches, if any; it is updated after rerun.
