# NPJ_sample6

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample6；patient: nan。

组织spots 3497，有效Epi 2414，内部Epi 1808，冻结P95 91。精确前沿Epi 957，前沿内部Epi 159，R12标注Epi 2072，未标注 342。

## 候选距离原始Pearson

- 核心面板较低r：GC B / HALLMARK_E2F_TARGETS: r=-0.539, n=91, 0–954 μm, spatial CI [-0.693, -0.361], P=3.49e-08, across-section q=1.8e-06。
- 核心面板较高r：Plasma / HALLMARK_UNFOLDED_PROTEIN_RESPONSE: r=0.541, n=87, 200–2498 μm, spatial CI [0.322, 0.686], P=6.41e-08, across-section q=3.08e-06。

核心面板最大绝对共同范围差异：Plasma / HALLMARK_OXIDATIVE_PHOSPHORYLATION，r_P95=-0.152 (n=87)，r_nonP95=0.333 (n=1496)，delta=-0.485，空间CI [-0.739, -0.151]，空间P=0.00188，跨片空间q=0.0452，共同范围200–2498 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.192；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.107。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.217；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.117。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.196；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.143。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.122；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.141。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 FLORES_CRISPR_PRO_INVASION_10, β=0.058；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.128。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.062, CI [-0.174, 0.299], n=689, 新85上下文族q=0.837。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.154, CI [-0.400, 0.092], n=689, 新85上下文族q=0.481。
- C2_P95_enrichment / P95_probability: β=0.068, CI [-0.037, 0.174], n=689, 新85上下文族q=0.433。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.073, CI [-0.470, 0.615], n=159, 新85上下文族q=0.794。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.016, CI [-0.190, 0.159], n=957, 新85上下文族q=0.433。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：35.16个百分点，CI [nan, nan]，标注P95 n=12、非P95 n=137，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
