# NPJ_sample4

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample4；patient: nan。

组织spots 2639，有效Epi 1885，内部Epi 1569，冻结P95 79。精确前沿Epi 701，前沿内部Epi 136，R12标注Epi 1655，未标注 230。

## 候选距离原始Pearson

- 核心面板较低r：CD8 GZMK/Effector-memory-like / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.751, n=79, 0–1825 μm, spatial CI [-0.874, -0.486], P=1.53e-15, across-section q=5.74e-13。
- 核心面板较高r：cDC2 / BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18: r=0.627, n=79, 0–1136 μm, spatial CI [0.439, 0.726], P=6.56e-10, across-section q=5.85e-08。

核心面板最大绝对共同范围差异：Classical Mono / HALLMARK_GLYCOLYSIS，r_P95=0.437 (n=79)，r_nonP95=-0.141 (n=1465)，delta=0.579，空间CI [0.186, 0.742]，空间P=5.53e-06，跨片空间q=0.00138，共同范围0–1082 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.315；较低 HALLMARK_E2F_TARGETS, β=-0.030。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 HALLMARK_GLYCOLYSIS, β=0.194；较低 HALLMARK_E2F_TARGETS, β=-0.031。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.192；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.050。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.206；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.075。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.251；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.107。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.231, CI [-0.487, 0.025], n=587, 新85上下文族q=0.25。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.065, CI [-0.111, 0.241], n=587, 新85上下文族q=0.683。
- C2_P95_enrichment / P95_probability: β=-0.015, CI [-0.035, 0.006], n=587, 新85上下文族q=0.228。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.510, CI [-1.070, 0.050], n=136, 新85上下文族q=0.153。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.125, CI [-0.068, 0.317], n=701, 新85上下文族q=0.153。

## 限制与入口

29候选中 23 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
