# NPJ_sample3

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample3；patient: nan。

组织spots 4881，有效Epi 3182，内部Epi 2955，冻结P95 148。精确前沿Epi 1207，前沿内部Epi 184，R12标注Epi 2767，未标注 415。

## 候选距离原始Pearson

- 核心面板较低r：GC B / HALLMARK_E2F_TARGETS: r=-0.538, n=148, 0–872 μm, spatial CI [-0.698, -0.295], P=1.75e-12, across-section q=2.89e-10。
- 核心面板较高r：cDC2 / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.721, n=148, 100–2598 μm, spatial CI [0.446, 0.833], P=4.66e-25, across-section q=2.12e-21。

核心面板最大绝对共同范围差异：Lymphatic EC / HALLMARK_GLYCOLYSIS，r_P95=-0.323 (n=148)，r_nonP95=0.538 (n=2767)，delta=-0.861，空间CI [-1.073, -0.534]，空间P=2.37e-09，跨片空间q=8.53e-06，共同范围200–3158 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.761；较低 BUDDEN_S1P_RESPONSE_18, β=-0.193。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.734；较低 BUDDEN_S1P_RESPONSE_18, β=-0.198。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.662；较低 BUDDEN_S1P_RESPONSE_18, β=-0.169。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.551；较低 BUDDEN_S1P_RESPONSE_18, β=-0.092。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.434；较低 BUDDEN_S1P_RESPONSE_18, β=-0.015。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.185, CI [-0.371, -0.000], n=1101, 新85上下文族q=0.179。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.002, CI [-0.163, 0.159], n=1101, 新85上下文族q=0.991。
- C2_P95_enrichment / P95_probability: β=-0.003, CI [-0.007, 0.001], n=1101, 新85上下文族q=0.316。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.221, CI [-0.211, 0.652], n=184, 新85上下文族q=0.349。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.044, CI [-0.090, 0.178], n=1207, 新85上下文族q=0.16。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
