# BLCA 病灶空间分区 v3：实际结果

已完成 **63/63 张、7 个 dataset、265,002 个原始 spot**。v2 保留原位。63张均有真实H&E对照、域解释、主分区、距离和组成图，逐点标签、ROI几何、两类独立边界以及100/200/300 μm敏感性结果。入口：[结果索引](index.html)。所有ROI仍为计算辅助候选，没有新增病理医师确认。

## 完成情况与可用边界

- **60 张**形成候选病灶ROI并有观测外缘；合计480个连通ROI，其中284个少于10点、合计661点。它们是保留的小候选，不应解读为已确认的480个独立肿瘤灶。
- **59 张**内、外边缘各至少20个可估计点，可作有条件的边缘比较；**57 张**有各至少20点的深部内部/远外部。标记关键解释冲突/叠图残差后，53张满足上述双侧覆盖且没有这些预设警示，仍不是病理认证名单。
- 有效主距离共 **152,739 点（57.6%）**。位置已确定但内部深度未估计72,599点、外部距离未估计34,036点，未把它们当深部core。
- N1：无支持病灶；N2：主要为肌性背景、31个局部候选冲突；N3：1,677点局部位置未定、812点显示无支持病灶。N3仍不适合以此版做肿瘤外缘分析。

| Dataset | 切片 | ROI切片 | 双侧边缘≥20点 | 可估计距离点 |
|---|---:|---:|---:|---:|
| CanDis__UNVERIFIED | 13 | 13 | 12 | 24,633 |
| GSE171351 | 4 | 4 | 4 | 3,193 |
| GSE238145 | 10 | 10 | 10 | 52,240 |
| GSE285715 | 6 | 3 | 3 | 6,744 |
| GSE319536 | 16 | 16 | 16 | 30,958 |
| Inhousedata__UNVERIFIED | 7 | 7 | 7 | 22,769 |
| NPJ__UNVERIFIED | 7 | 7 | 7 | 12,202 |

## v2到v3的实质变化

v2主要描述细胞/组织身份；v3先用连通BANKSY域组织证据，再建立病灶范围。内部免疫/间质岛可保留原身份并处于ROI内，组成不确定不再自动使空间位置不可计算。

原 Other tissue / Uncertain-conflicting 共 **94,926 点**：

- 18,732点归入病灶内；
- 74,517点归入病灶外；
- 1,677点保留局部不可定位（N3）。另有N2的31个原Tumor候选转为不可定位。

最终ROI内133,725点，其中22,781点并非v2 Tumor身份；8,152个v2 Tumor点处于v3 ROI外。内部域/包围结构归并共3,033点，其余差异来自域级范围识别。ROI内有26,178点低于旧Epi≥0.25入口；v3不依赖该gate。

这些差异是分析对象改变与计算注释调整，不能称为“纠正了多少真实肿瘤”。正常尿路上皮、免疫、基质和平滑肌等v2组织标签仍在独立字段，不把所有outside统称正常上皮。

## 主外缘和局部接触面

已分别保存 `lesion_outer` 和 `tumor_nest_contact`。主候选外缘24,459段；次级观测接触37,753段，其中18,456段两侧均处于主ROI内部。这些内部接触没有混入主外缘。两套几何、spot距离和MP10/TME摘要均单独输出。

次级接触继承v2计算恶性身份，不代表新的病理确认。缺损/空白不填补，捕获/采样截断和组织表面候选边不作为生物界面。详见METHODS中的删失定义。

## 优先审阅图

| 切片 | 推荐检查内容 | 图 |
|---|---|---|
| Inhousedata_T4463209 | 斜向连续外缘、肌性外部与上方截断 | [H&E与ROI](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4463209/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4463209/02_v2_SpaCET_v3_comparison.png) |
| Inhousedata_T4413900 | 内部混合岛归位；内部接触与外缘分离 | [H&E与ROI](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4413900/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4413900/02_v2_SpaCET_v3_comparison.png) |
| GSE319536_BC01 | 多灶/碎片、真实空白未跨越 | [H&E与ROI](by_dataset/GSE319536/GSE319536_BC01/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/GSE319536/GSE319536_BC01/02_v2_SpaCET_v3_comparison.png) |
| CanDis_2026_Sample8 | 全diploid冲突、大病灶和广泛深度删失 | [H&E与ROI](by_dataset/CanDis__UNVERIFIED/CanDis_2026_Sample8/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/CanDis__UNVERIFIED/CanDis_2026_Sample8/02_v2_SpaCET_v3_comparison.png) |
| GSE285715_N1 | 没有强制生成病灶的癌旁对照 | [H&E与ROI](by_dataset/GSE285715/GSE285715_N1/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/GSE285715/GSE285715_N1/02_v2_SpaCET_v3_comparison.png) |
| GSE285715_T3 | 低置信证据下的候选域与局部边缘 | [H&E与ROI](by_dataset/GSE285715/GSE285715_T3/01_v3_core.png) / [v2–SpaCET–v3](by_dataset/GSE285715/GSE285715_T3/02_v2_SpaCET_v3_comparison.png) |

每片5组PNG/PDF位于对应dataset目录；全队列7个总览在 `figures/*_overview.png/pdf`。主图没有Mixed图例；浅红/浅蓝代表内外位置保留而距离删失，灰色仅表示局部不可定位。

## 具体待复核区域

- **GSE285715_N3**：局部 1,677 点来源/上皮恶性证据冲突；未建立可信病灶ROI；其余812点也不是独立病理阴性证明。
- **GSE285715_N2**：31个原Tumor候选位置未定位；肌束区域的SpaCET/CNA高值不作已确认病灶。
- **GSE285715_T3**：低置信diploid与expr:BLCA依赖；左侧过渡域、10个候选ROI及截断区需复核。
- **GSE319536_BC01**：18个候选ROI，中央及周缘碎片/腔隙；小ROI和跨缺口的距离被保留/删失。
- **CanDis_2026_Sample8**：原生全diploid但模型/轮廓候选病灶；大范围内部深度删失，优先看上方观测外缘。
- **CanDis_2026_Sample11**：小且截断的133点候选ROI；没有足够深内部，内侧可估计边缘不足20点。
- **CanDis_2026_Sample1**：图像仿射残差最大约0.258 pitch；保留阵列距离，局部图像轮廓需复核。
- **CanDis_2026_Sample5**：图像仿射残差最大约0.246 pitch；保留阵列距离，局部图像轮廓需复核。

完整候选ROI位置、坐标范围与局部问题见 `tables/lesion_roi_summary.tsv`、`local_review_targets.tsv`，简短优先清单见 `priority_review.tsv`。小ROI保留但需审阅，不能把全部连通小岛当作已确认病灶。

## 带宽与工具对照

固定ROI后，主Margin点数：100 μm半宽 **23,334**；200 μm **42,768**；300 μm **64,662**。主距离/ROI完全相同，变化仅为分箱；200 μm是项目候选尺度，不是Gut公布参数。比例变化明显，因此正式结论应结合连续距离和三个尺度，不只依赖单个分带。

Cottrazm的BoundaryDefine模块实际运行三片：

| 切片 | 已赋值 / 总点数 | 未赋值 | 图 |
|---|---:|---:|---|
| Inhousedata_T4413900 | 4,385 / 4,969 | 584 | [诊断图](comparators/Inhousedata_T4413900/comparison.png) |
| GSE319536_BC01 | 1,315 / 5,278 | 3,963 | [诊断图](comparators/GSE319536_BC01/comparison.png) |
| GSE285715_T3 | 892 / 3,341 | 2,449 | [诊断图](comparators/GSE285715_T3/comparison.png) |

其输入经过明确适配，且BC01/T3未赋值较多；不是完整Cottrazm复现，也不是准确率评估。N1触发“无可接受恶性种子”保护，没有机械选择最高CNV群作肿瘤。SPATA2未运行：交付的是可加载的sf/GeoJSON等价几何实现及独立有效边界距离。

## MP10与TME：仅用于后续衔接的描述性结果

在 **57张**具备内部/外部覆盖的切片中，MP10区域中位数“内部−外部”为正的有 **46张**，负的9张，零差2张；全切片配对差中位数 **0.0200 AUCell**。同时全部57张的内部恶性估计中位数更高，所以不能把MP10总体升高直接解释成肿瘤细胞内的状态增强。

Treg、EC、CAF在内部较低的切片分别为 **45/57、38/57、45/57**。这描述的是区域组成差异，不是独立验证；域注释本身使用了部分组成和相同表达证据。

按核验/匿名分析单位合并重复切片后的方向一致性：

| Dataset | 单位数 | MP10内部更高 | Treg内部更低 | EC内部更低 | CAF内部更低 |
|---|---:|---:|---:|---:|---:|
| CanDis__UNVERIFIED | 12 | 9 | 11 | 8 | 10 |
| GSE171351 | 3 | 2 | 2 | 0 | 3 |
| GSE238145 | 9 | 8 | 5 | 2 | 5 |
| GSE285715 | 1 | 1 | 0 | 1 | 1 |
| GSE319536 | 16 | 12 | 16 | 12 | 12 |
| Inhousedata__UNVERIFIED | 7 | 7 | 5 | 6 | 5 |
| NPJ__UNVERIFIED | 7 | 5 | 5 | 7 | 7 |

可见MP10内部较高的方向在多个dataset较常见；Treg外侧较高在GSE319536最一致。EC的方向具有明显dataset差异：GSE171351和GSE238145更多为内部较高，不能泛化成统一的“血管外侧富集”。CAF多为外侧更高，但并非每片一致。GSE285715的3张肿瘤切片属于同一核验分析单位，不能算3个独立病例；未核验公开编号的来源也不冒充已核验患者。

这些内部—外部结果不等于“在Margin达到峰值”。已经提供内/外Margin独立摘要、恶性比例分层的原MP10摘要和单独的内部接触摘要；下一步可在选定可评估ROI/边段上进行患者/分析单位层面的连续距离曲线、100/200/300 μm稳健性和比例控制后的比较。缺失深度不是随机缺失，应排除或单独处理，不将浅红位置补为core。

[切片差值图](figures/MP10_TME_section_contrasts.png)；[区域摘要](tables/lesion_outer_region_MP10_TME_summary.tsv)；[内部接触摘要](tables/tumor_nest_contact_MP10_TME_summary.tsv)。

## 交付、检查与资源

- `tables/all_spots_v3.parquet`：全265,002点完整字段，包括原有证据/组成、MP10及新几何；`all_spot_labels_v3.tsv.gz`是便于审阅的核心列版。每片另有完整TSV.gz。
- `tables/domain_annotations.tsv.gz`：连通域汇总；`lesion_roi_summary.tsv`：ROI汇总；`section_summary.tsv`：切片状态。
- 每片 `lesion_rois_um.geojson`、图像像素版本、`lesion_rois_sf.rds`（有ROI的60片）、独立两类边界几何和距离表。无ROI切片保存空GeoJSON和明确状态。
- `METHODS.md`、`config.json`、`scripts/`、实际Cottrazm适配源码和运行版本：可续跑和复现，不依赖隐藏手工编辑。
- 63张坐标/符号/标签检查全部通过，原MP10完整匹配；所有核心图300 dpi PNG并有PDF，63张H&E路径可用。
- 主流程两片并发、每进程一个BLAS线程。63片逐进程实测峰值最大644.8 MiB；不应将该数当作整个Python汇总或R诊断的峰值。本轮没有接近200 GB预算，也没有重建大型表达/CNA缓存。结果目录约1.55 GiB（不含可选review.zip），以交付时磁盘实际大小为准。

复现命令见METHODS；在summarize_v3之后运行 `python3.12 results/mp10_pathology_regions_v3/scripts/write_report.py` 更新本报告和优先复核表。
