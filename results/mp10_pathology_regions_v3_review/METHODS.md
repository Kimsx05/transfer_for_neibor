# BLCA 病灶空间分区 v3：方法与复现

## 分析对象和证据层级

本版的主对象是**计算辅助病灶 ROI 及其外缘**，不是逐 spot 恶性细胞身份。ROI 可以包含内部间质、免疫浸润和混合组成。表内的 `malignant_fraction`、`CNV_support`、原 v2 组织身份、30 类组成和证据等级独立保留。图中的 Tumor interior 指病灶范围内、距可评估外缘超过指定距离的位置，不表示其中所有细胞恶性。

没有新增或虚构病理医师确认、复核人、复核日期。代表片进行了图像与计算结果的结构检查；这不等于病理专家验证。所有 ROI 标记 `computational_candidate_no_documented_pathologist_confirmation`。模型之间共享表达信息，不能把 SpaCET/CNA/c2l 一致性称为独立实验验证。

v2 整体保持原位，只读取其结果。实际清单来自 `mp10_pathology_regions_v2/tables/section_manifest.tsv`，63 张 ordinary Visium、7 个 dataset、265,002 个原 QC 位置。保留 dataset、section、原始 barcode；每次连接检查键唯一性、位置数及覆盖。v3 不新增表达过滤，也不重跑 cell2location、MP10、CopyKAT、SpaCET 或 BANKSY 全队列。

## 复用输入与软件边界

- 复用 v2 SpaCET 连续 malignant fraction、major-lineage 组成、实际 signature 来源及官方默认 0.5 的 Tumor/Interface/Stroma。原生 Interface 是紧邻肿瘤的非肿瘤一环；不是本版 200 μm margin，也不是组织学侵袭前沿。
- 复用 v2 实际 CNA profile 分析得到的 cluster 支持、染色体事件及原生 CopyKAT 标签，包括 `c1:diploid:low.conf`。不重新缩放 CNA 或套用其他论文不同 score 的阈值。v3 不再次加载大型 profile，因为本轮复用的是已完成的轮廓级证据；原 profile 仍由 v2 清单定位。
- 复用 Banksy 1.9.5 的 `clust_*res0.6` 空间域（v2 全转录组、lambda 0.2、k=30 的既有结果），不重新调参，不将 cluster 数强制成三类。识别域列使用 `clust_`，不误将 `cluster_consensus_MSQ` 当成域。
- c2l 使用原 `v2_30type` mean 表，逐 spot 除以 30 类 mean 总和得到 fraction；没有把 q05 混入 mean，没有将总上皮改名恶性。原始 mean 和派生 fraction 同时保存。
- 原 MP10 取冻结文件中的 **`EpiMP10_Top20_AUCell`**，数值不重算。只有 ROI 和边界完成后才按 barcode 接入。试点规则冻结发生在 MP10 汇总之前。
- 现有 SPATA2 3.1.4 接口文档已经核对；当前环境没有安装 SPATA2，缺少 concaveman、ggalt、keys、shinyhelper、tidytext 等 imports。为了避免升级已有分析环境，实际使用 **Shapely + sf 的局部多边形适配**，保存 GeoJSON 和 `lesion_rois_sf.rds`。没有声称调用了 SPATA2 的注释或距离函数。SPATA2 `createGroupAnnotations` / `getCoordsDfSA` 所针对的完整注释边界，与本项目剔除截断/缺损的有效边界并不相同，不能直接混用距离字段。
- Morph 和 POSST 仅用于参考分层边界、组织缺损和对象内距离的原则，未运行这两套模型。

本版是**参考 Gut 的“空间域 → 组织区域注释 → 区域内组成分析”思路的项目实现**。BANKSY v0.1.4/resolution 0.6 等文献信息不足以完整复现作者 ROI；未取得完整 domain→三类区域脚本，不声称完整复现。以下打分、ROI 归并和距离删失均为本项目规则。

## 域级组织证据和 ROI

首先按真实六邻接，将每个原 BANKSY 域拆成空间连通部分。邻居阵列偏移为 `(0,±2)、(±1,±1)`；不同切片和无捕获点之间不连边。同一表达域的分离区域有独立 `domain_component_id`，同时保留原 `spatial_domain`。

各连通部分汇总 malignant fraction 四分位数、CNA 支持、上皮/免疫/基质/平滑肌组成、尿路上皮 marker、H&E 前景、原组织解释及相邻域。域候选证据分数是可解释的启发式分数，**既不是恶性概率，也没有独立训练/验证**：

`0.65 × malignant_fraction + 0.20 × CNA_support_score + 0.15 × clipped_epithelial_program − 0.45 × contractile_context`

逐点分数在域连通部分内平均。CNA support strong=1、moderate=0.55、baseline-caution=0.35、coherent-native-diploid-conflict=0.30、weak/unavailable=0；epithelial program 为已有 marker_epithelial / 2.5 后限制到 [0,1]。contractile context 要求 c2l SMC≥0.5、至少 3 个平滑肌 marker 被检测且 epithelial program 原分数<0.8，避免 N2 类肌组织模型冲突直接成为病灶。少于 5 点的域部分用 25% 局部均值 + 75% 同组织连通区内同原域均值稳定估计；不删除小岛。

分数≥0.42 注为 Tumor-related domain，0.22–0.42 为 Interface-related domain，更低为 Non-tumor-related domain。这是固定的候选注释起点，不是文献标准阈值。不用 MP10、Treg、CAF 亚型、EC 富集或显著性挑选阈值。

初始 ROI 采用 Tumor-related 域部分；随后仅作以下有记录的结构归并：

1. Interface-related 域部分若至少 60% 跨域邻接边接触病灶域，且不接触组织外侧空白，可归入 ROI。
2. 非 ROI 的已观测连通组织若被病灶包围、不连到外侧空白，且存在 ROI 邻居、H&E 前景中位数≥0.35、nFeature 中位数≥500，可作为表达支持的内部组织归入 ROI。若图像缺失仅保留表达/几何候选依据，状态明确记录；本批实际均有图像。
3. H&E 前景为原 spot 半径附近像素 saturation>0.07 或平均 optical density>0.15 的比例，只用于区分可见组织和空白，**不是自动识别坏死或病理亚型**。该规则不能排除局部坏死混入，所有归并区域仍是可复核候选。
4. 不填补没有 QC spot 的洞，不作大凸包，不反复闭运算、删岛或填洞，不默认取最大 ROI。保留所有连通病灶及小候选，分别编号。真实腔隙/缺损仍为空白或缺损边界；未知缺损不命名坏死。

`assignment_basis` 保存初始域依据和内部组织归并理由；`domain_annotation` 保留组织域解释，不因进入 ROI 就改写成距离 Margin。

### 已知癌旁背景

N1 的表面尿路上皮与深部肌性组织、N2 的肌束结构及多 marker/c2l 背景沿用既有图像解释，因此不把模型冲突自动升级为病灶。N1 无支持病灶；N2 的 31 个原 Tumor 候选局部保留 unlocalized。N3 将既有明确冲突位置保留 unlocalized，其余位置显示无支持病灶。N3 全片病灶存在性仍未解决，不能把其图中的非灰位置当作独立病理阴性证明。

这些是结合来源/形态信息的注释决策，不是盲法特异性测试；SpaCET 原始高值没有被改成零。正常上皮的独立 v2 标签保留，但不据“在 ROI 外”重新命名正常上皮。

## 几何、物理尺度和有效外缘

ordinary Visium 官方中心间距为 100 μm、捕获点直径 55 μm。根据已核对的 array_row/array_col 定义局部物理坐标：

`x_um = array_col × 50; y_um = array_row × sqrt(3) × 50`

这是厂家标称阵列标定，不是显微镜独立实测组织缩放。图像用同一阵列到已有 lowres 图像坐标的仿射拟合，不对已对齐 x/y 再套 scalefactor。每片保存拟合矩阵、残差、像素间距与尺度来源。CanDis Sample1/5 的最大残差超过 0.15 pitch，列为局部叠图复核；物理阵列距离仍按真实网格计算。

ROI 是所纳入 spot 的 Voronoi 六边形并集。边界零点在共享六边形轮廓，**不在边缘 spot 中心**。相邻中心 100 μm，共享轮廓到各中心 50 μm。多边形表示观测网格的区域足迹，不声称达到细胞级 H&E 描边精度。面积是纳入 spot 六边形面积之和，不是组织学精确面积。

保存所有 ROI 边段并分别分类：

- `biological_lesion_outer`：ROI 与外侧连通、可定位、已观测非 ROI 组织之间的**计算候选**生物外缘。
- `internal_compartment_boundary`：内部未纳入组织的边界，不进入主外缘。
- `QC_gap`：原全组织坐标存在、当前 QC 表缺失的位置。
- `internal_tissue_defect`：内部缺损/空白，无可靠生物界面证据。
- `tissue_cut_or_surface_edge`：外侧无观测 spot 且邻近图像为空白；不区分天然表面和切片切缘。
- `capture_or_sampling_truncation`：外侧无观测 spot，但预测位置图像仍有组织。
- `image_or_capture_truncation`、`unresolved_capture_or_tissue_edge`：超出图像或类型无法判断。
- `local_contour_uncertain`：邻接局部位置不可判断或 QC 排除。

图像支持的截断细分是候选解释，不能声称已恢复仪器捕获区精确轮廓。上述非生物边段均不作为主边界。没有将组织外缘/捕获边自动称为侵袭前沿。

主距离是在同一保留组织连通区内，到有效 ROI 边段的精确点—线段欧氏距离；内部仅匹配自己的 ROI，外部保存匹配到的 `reference_lesion_roi_id`。内部负、外部正。分块计算，不构造全队列稠密距离矩阵。

若无有效外缘，保存 `boundary_not_observed`。若到未观测/缺损边的距离不大于到有效边的距离，主距离删失为 `censored_by_unobserved_or_defect_edge`。若直线路径每 25 μm 取样到同组织连通区最近中心超过 58 μm（六边形覆盖半径约57.7 μm），标记 `occluded_by_tissue_gap`，不输出跨空白距离，也不偷偷换成测地距离。保留到观测外缘的候选值在单独字段，不能用于主距离分析。没有有效距离仍可可靠保存 inside/outside。

这种删失较保守，特别是组织碎片、广泛病灶或缺口较多样本，可估计比例会下降。其意义是避免假精确核心深度，不能把删失当成随机缺失。

## 分带与次级接触面

主半宽 w=200 μm；敏感性仅100/200/300 μm，ROI不随带宽改变：d<−w 为 Tumor interior，−w≤d<0 为 Inner margin，0≤d≤w 为 Outer margin，d>w 为 Outside。主图合并两侧为 Margin；表保留两侧、连续距离和域注释。

无法估计深度的内部用 `Inside (depth unestimated)`（浅红），外部用 `Outside (distance unestimated)`（浅蓝）。局部不可定位为灰，QC问题另设深灰；主图无 Mixed/Other。无支持病灶的样本明确标注，不强造 Margin/core。

次级 `tumor_nest_contact` 使用 v2 计算 Tumor 身份的连通巢与相邻可解释非 Tumor 组织之间的实际边段，作为显式继承 v2 的局部接触比较对象；不是 v3 ROI 的真值。每巢独立编号，inside只到自己的巢边，outside只在同组织连通区并检查路径。距离是**到已观测接触面**的距离，没有声称恢复完整巢的未观测外缘。N2 残留候选接触面带 unresolved 状态，不能作已确认肿瘤界面。是否位于主 ROI 内独立保存。

次级表、几何和摘要与主外缘完全分开。后续接触摘要限于主 ROI 内观测到的接触距离，不能与主 signed_distance_um 直接合并。

## 试点与 Cottrazm 对照

8片：T4413900、T4463209、CanDis Sample8、BC01、N1、N2、T3、B13，覆盖清楚外缘、内部混合、碎片/多灶、癌旁和截断。所有片均生成 v2 / 原生SpaCET / v3 同坐标并排图。按结构检查固定规则后扩展63片；没有按目标富集筛选切片。

实际执行 Cottrazm 0.1.1 的 BoundaryDefine 模块（源码 commit `888afbb69cec007ee5e597a8722d8453fa6e2bc3`），在 T4413900、BC01、T3 三片诊断性对照。**没有运行其完整形态调整表达和 inferCNV 流程**。适配如下：

- 用真实原 counts 建立 Seurat Spatial 对象；缓存 BANKSY PCA→UMAP，原 CNA_cluster 作 CNVLabel，原 BANKSY 域作 seurat_clusters，1−SpaCET fraction 作 NormalScore。这些属于本项目替代输入。
- 显式选择 malignant≥0.65、检测≥3 epithelial markers且epithelial score≥0.8的候选；CNA cluster中此类点≥50%且平均malignant≥0.5才作种子。没有默认选择最高CNV两组。癌旁 N1 执行到种子保护步骤即停止，不作为零假阳性的性能验证。
- 原 dense Manhattan 邻接替换为真实六邻接稀疏表；修正单种子时 R 子集降维；移除 Seurat5 不兼容的绘图块，主 BoundaryDefine 迭代逻辑仍来自缓存源码。每片保存实际适配脚本、种子、标签和 sessionInfo。
- 未赋值位置保留 Unassigned，不能为提高一致率强行补齐。该比较仅诊断分区对象/输入的影响，模型一致率不称准确率，不外推全队列性能。

## MP10 和组成摘要

规则冻结后输出每切片、每分带的原 MP10 和组成均值/中位数/四分位数。Treg 为原c2l Treg；EC为四类血管/淋巴EC之和；CAF为Matrix_CAF/iCAF/rCAF/vCAF之和。使用同批组成注释区域后的组成差异不是独立发现。

配对比较要求两区各≥20个可估计距离点，只取切片区域中位数差，不将 spots 当病例。映射沿用已核对的 patient/analysis-unit crosswalk；同分析单位多切片差值先平均。匿名来源单位明确保留 `__UNVERIFIED` 或对应身份状态，不冒充已核验患者。

另给原恶性估计0–.25/.25–.5/.5–.75/.75–1分层中的MP10描述，帮助判断总体梯度是否随肿瘤比例变化。但这些是共享表达推断的粗分层，**不是肿瘤细胞内状态的纯化或因果校正**。本轮不开展CellChat/COMMOT和机制推断。

## 运行与检查

主流程两张切片并发，每进程1个BLAS线程；按片保存，`complete.json` 存在则跳过，重算加 `--force`。无新增大型表达/CNA缓存。当前项目不是可用Git仓库，因此项目脚本以本目录实际文件为准，不虚构commit；第三方版本和运行session单独保存。

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/run_cohort.py
# 单片重算示例
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/run_v3.py Inhousedata_T4413900 --force
# Cottrazm诊断（逐片运行）
/usr/local/bin/Rscript --vanilla results/mp10_pathology_regions_v3/scripts/cottrazm_diagnostic.R Inhousedata_T4413900
/usr/local/bin/Rscript --vanilla results/mp10_pathology_regions_v3/scripts/cottrazm_diagnostic.R GSE319536_BC01
/usr/local/bin/Rscript --vanilla results/mp10_pathology_regions_v3/scripts/cottrazm_diagnostic.R GSE285715_T3
# 完成队列后生成总表/图、sf对象和必要几何检查
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/summarize_v3.py
/usr/local/bin/Rscript --vanilla results/mp10_pathology_regions_v3/scripts/export_sf.R
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/check_geometry.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/write_report.py
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 results/mp10_pathology_regions_v3/scripts/finalize_review.py
```

R库路径沿用v2独立Rlib、用户R4.5/R4.4和/refdir/Rlib；Python优先使用项目现有 `spatial_ecology_top20_v2/runtime`。版本详见 `logs/python_versions.json`、`logs/sf_sessionInfo.txt`、每个Cottrazm目录sessionInfo；随机种子20261029。主几何算法确定性。几何检查包括已知1000μm正方形内外符号、100μm中心/50μm轮廓关系、ROI有效性、无重复barcode、距离状态与数值一致、原MP10完整匹配。

## 核对的官方资源

- [SpaCET 官方](https://github.com/data2intelligence/SpaCET)，[interface说明](https://data2intelligence.github.io/SpaCET/reference/SpaCET.identify.interface.html)。
- [BANKSY 官方](https://github.com/prabhakarlab/Banksy)。
- [SPATA2 创建注释](https://themilolab.github.io/SPATA2/articles/creating-spatial-annotations.html)，[使用空间注释](https://themilolab.github.io/SPATA2/articles/using-spatial-annotations.html)。
- [Gut文章](https://doi.org/10.1136/gutjnl-2025-335084)，[补充方法](https://gut.bmj.com/content/gutjnl/early/2025/07/30/gutjnl-2025-335084/DC1/embed/inline-supplementary-material-1.pdf?download=true)。
- [Morph](https://github.com/ding-lab/morph)，[POSST](https://github.com/zaozaozaonan/POSST)。
- [Cottrazm官方](https://github.com/Yelab2020/Cottrazm)，[论文](https://doi.org/10.1038/s41467-023-36560-7)。
- [10x Visium平台点径/间距说明](https://cdn.10xgenomics.com/image/upload/v1758134784/support-documents/CG000686_VisiumHDSpatialApplicationsPerformance_TechnicalNote_RevB.pdf)，[Visium CytAssist 11mm说明](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000495_VisiumCytAssist_GeneExpressionUserGuide_RevA.pdf)。
