# Inhousedata_T4413900

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4413900；patient: nan。

组织spots 4969，有效Epi 3842，内部Epi 3818，冻结P95 191。精确前沿Epi 2026，前沿内部Epi 226，R12标注Epi 3559，未标注 283。

## 候选距离原始Pearson

- 核心面板较低r：Naive/Memory B / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.527, n=191, 100–5977 μm, spatial CI [-0.658, -0.356], P=4.88e-15, across-section q=1.55e-12。
- 核心面板较高r：iCAF / BUDDEN_TNF_EMT_RESPONSE_15: r=0.429, n=191, 0–1670 μm, spatial CI [0.249, 0.568], P=6.03e-10, across-section q=5.17e-08。

核心面板最大绝对共同范围差异：Naive/Memory B / BUDDEN_S1P_RESPONSE_18，r_P95=-0.247 (n=191)，r_nonP95=0.188 (n=3625)，delta=-0.435，空间CI [-0.590, -0.262]，空间P=1.33e-06，跨片空间q=0.000555，共同范围100–5977 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.375；较低 PURAM_PEMT_TOP30, β=-0.160。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.265；较低 PURAM_PEMT_TOP30, β=-0.156。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.188；较低 PURAM_PEMT_TOP30, β=-0.115。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.120；较低 PURAM_PEMT_TOP30, β=-0.068。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.104；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.005。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.171, CI [-0.312, -0.030], n=2004, 新85上下文族q=0.0844。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.012, CI [-0.353, 0.329], n=2004, 新85上下文族q=0.972。
- C2_P95_enrichment / P95_probability: β=0.044, CI [-0.005, 0.093], n=2004, 新85上下文族q=0.0979。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.791, CI [-1.109, -0.474], n=226, 新85上下文族q=6.28e-05。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.162, CI [0.091, 0.233], n=2026, 新85上下文族q=6.28e-05。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：11.26个百分点，CI [-7.95, 29.78]，标注P95 n=15、非P95 n=183，跨片q=0.531。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 27 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
