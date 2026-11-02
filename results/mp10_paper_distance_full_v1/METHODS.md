# 全队列空间距离分析 (full v1)

## 研究对象、位置和分数

登记对象来自 mp10_pathology_regions_v2 的冻结 section_manifest，使用其中的输入路径。正式区域直接来自 v3.4 final_spot_labels，保留 lesion_location、region_role、tissue_component_id、signed_distance_um、精确距离状态和 distance_band_200um。图像展示采用配准的低分辨率 H&E、原 x/y 及 v3.4 像素边段；目标距离单独以 x_um/y_um 计算。所有逐点连接使用 section_id 和完整 barcode，保留 position_id、spot_uid 与后缀。患者、分析单位、来源和可用身份政策取自原 section_identity_crosswalk；[输入表](tables/input_manifest.tsv) 提供路径和来源依据。

主查询直接使用原 analysis_eligible_legacy、可用 EpiMP10_Top20_AUCell 和未标记 quality_excluded 的位置。该值是上皮支持位置的 **mixed-spot** AUCell。保留有效零值、原始模型 universe 分数、epi_score_valid、common_anchor、score_eligible 和 E0。旧七张先导属于同一 manifest，各进入一次。

## 距离与条件

采用原 pilot 对 Budden 等 Oral/cSCC 最近目标 spot 中心 min(Euclidean) 的适配。每片在未质量排除、有 center_qc、所有 30 项 posterior mean 有限非负且总和大于零的组织 spot 上计算组成。CAF 汇总 Matrix_CAF、iCAF、rCAF、vCAF；TAM 汇总 C1QC TAM 和 SPP1 TAM；Treg 为单因子；分母是该位置全部 30 项和。FOLR2 没有单独空间因子。按片取线性分位数，纳入比例 ≥ 阈值且 >0 的所有点，保留并列。无变异或没有有效组成时不可定义。target 内估计丰度中位数 <1 仅作为探索性标记。原始作者脚本保存在相邻 pilot 的 authors/。

每片每个组织块使用精确 cKDTree，重叠为零，缺目标为 NA 并存原因。病灶精确有符号距离 (内负外正) 与目标距离独立。P80 全组织含零为主条件；分别比较 P80 全组织正距离、P90 全组织含零、P80 与原 inner_margin/outer_margin 相交的目标含零。margin 沿用全组织阈值。每组计算 Pearson r、Spearman ρ；至少三个有效配对且两变量可变才计算。正值表示离目标越远，评分越高。分箱依次为 0、(0,200]、(200,500]、(500,1000]、>1000 μm。端点四舍五入到六位小数仅用于分箱。分层关系沿用至少 20 有效配对的显示规则；common_anchor/E0 子集单独保留，同样至少需 20 个有效配对才展示 r。所有分析均以切片为相关计算单位。

**固定调整**仅用于 P80 全组织含零。用同一完整集合，对 mixed-spot MP10 和距离分别以截距、该 query 的 Epithelial/all30、log1p(nCount_RNA)、v3.4 lesion_location 哑变量作最小二乘拟合。深度沿用冻结 v8 nodes 的 nCount。两套残差以 Pearson 相关，另报同集合未调整 r、实际协变量及缺失原因。跳过常量列，残余自由度不足或残差近乎恒定时保留 NA。模型层对照只用原 query 内距离、spot_MP10_model_universe、epi_MP10_model_universe 均可用且 epi_score_valid 为真的交集，两个模型层用完全相同位置。调整不处理全部空间依赖，结果是描述性敏感性。

## 分区、niche 和程序

分区组成对各有效 spot 比例等权平均；Treg/T_NK 在每个 spot 先按正式 major_map 计算，再取分区均值。R12 (SAMPLE_R12_LEIDEN 0–8，对应旧展示 1–9) 的分母是该区有 R12 标签的合格 Epi anchors；dataset-local 分母是该区有该版本标签的有效组织点。分别列出 eligible、labelled、coverage；分母为零时比例为 NA。局部编号只在所属 dataset 内比较。R1/R2 参考环 profile 沿用原 R12 q05 汇总，与本轮中心 spot posterior mean 分开。各 niche 的区域分布另以同片同 niche 总数作分母。

两个版本主 MP10_median 一律以该区该标签中冻结 Epi queries 的有效 MP10 为统计集合；单列 n_MP10_valid，至少 5 点才给中位数。缺标签、无外缘、未定位与远区状态分别保留。分区程序分别针对有效组织 spots 和 Epi queries，给出 n_valid、中位数、Q1、Q3。区域的通用 MP10 摘要也取 Epi queries。

直接读取 pilot 冻结的 gene_sets.tsv.gz，human MSigDB 2026.1.Hs，50 Hallmark、3 Reactome、2 C7。逐片读取完整 gene × spot RNA counts，重复 symbol 合并计数，成员去重并与完整表达轴取交集。复用原七片 UCell 缓存，其余片用 UCell 2.14.0、maxRank 1500、平均并列秩、chunk.size 50、ncores 1。记录原基因数、匹配数及比例、状态和 universe 大小。首两张原 pilot 的三项 AddModuleScore 仅保留已有比较和逐点分数。C7 DN 表示 activated 相对高表达，UP 表示 resting 相对高表达。程序是混合 spot 转录量，未按细胞丰度作除法。胆固醇相关与 MP10 为共享表达及可能共享成员的程序一致性；Treg C7/CTLA4 与其组成比较在全组织进行，程序–目标距离在 Epi queries 计算。

## 分层汇总与图形

所有效应逐片计算，肿瘤研究来源与邻近正常始终分开。对有既有可用单位的同单位切片，取 Fisher-z 后等权平均并反变换。单位等权与 verified patient 子集分别输出；后者要求既有患者政策资格和可核验 patient_uid。GSE238145_B5/B18 使用原 P1；未解析的 GSE285715 T/N 不配对或造新单位。逐 dataset 报告切片中位数、四分位数、正负方向数；另计算单位等权效应。跨 dataset 描述对有有效结果的 dataset 等权，先对各 dataset 描述值做 Fisher-z 平均，再反变换，同时列出 dataset 中位数与覆盖。Fisher 变换在精确 ±1 处以 1−1e−12 作数值截断。无 spot-level P 值，无按 spot 数精度加权或跨片 spot 池化。

每片复用 pilot 空间面板，给出三条路线和全部 55 个功能的紧凑效应图；14 项重点程序制图。跨片同一特征采用全队列统一展示范围，常用图从 0 到该特征 q99，比例、热图等固定范围明确标注。适配 jins-routine-plotting 的全局白底、Arial、矢量矩形色标、拆分导出规范；区域组成独立图直接调用其 Stacked_composition 模板和 jin_export。专用空间图册按用户要求沿用 pilot 面板。只有原 H&E 作为位图嵌入 PDF。总览中的空间示例直接使用原七片；各 dataset 图册收录每张切片的全部图页。

计算使用 Python 3.12、NumPy 1.26.4、SciPy 1.11.4、pandas 2.1.4、PyArrow 21.0.0、Matplotlib 3.6.3，R 4.5.3、UCell 2.14.0。逐片读取；UCell 采用 6 个按片分工的 R 进程、每进程 1 个 UCell 核心和最多 2 个 data.table 线程。绘图用少量按片进程，底层 BLAS/OMP 设为 1。
