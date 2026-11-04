# CanDis_2026_Sample2

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample2；patient: nan。

组织spots 10062，有效Epi 4928，内部Epi 4363，冻结P95 219。精确前沿Epi 2460，前沿内部Epi 335，R12标注Epi 4580，未标注 348。

## 候选距离原始Pearson

- 核心面板较低r：Matrix_CAF / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.688, n=205, 1908–9629 μm, spatial CI [-0.770, -0.564], P=4.15e-30, across-section q=2.08e-25。
- 核心面板较高r：GC B / BUDDEN_TNF_EMT_RESPONSE_15: r=0.733, n=205, 0–2343 μm, spatial CI [0.605, 0.818], P=8.09e-36, across-section q=8.11e-31。

核心面板最大绝对共同范围差异：Pericyte / BUDDEN_S1P_RESPONSE_18，r_P95=0.543 (n=205)，r_nonP95=0.046 (n=4085)，delta=0.497，空间CI [0.353, 0.624]，空间P=6.73e-09，跨片空间q=1.64e-05，共同范围100–2960 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.418；较低 PURAM_PEMT_TOP30, β=-0.187。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.311；较低 PURAM_PEMT_TOP30, β=-0.211。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.306；较低 PURAM_PEMT_TOP30, β=-0.206。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.269；较低 PURAM_PEMT_TOP30, β=-0.178。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=0.158；较低 PURAM_PEMT_TOP30, β=-0.127。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.119, CI [-0.248, 0.010], n=2226, 新85上下文族q=0.177。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.072, CI [-0.410, 0.266], n=2226, 新85上下文族q=0.808。
- C2_P95_enrichment / P95_probability: β=0.007, CI [-0.011, 0.025], n=2226, 新85上下文族q=0.757。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.083, CI [-0.338, 0.172], n=335, 新85上下文族q=0.757。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.048, CI [-0.128, 0.032], n=2460, 新85上下文族q=0.000398。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：22.36个百分点，CI [nan, nan]，标注P95 n=5、非P95 n=313，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 26 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
