# CanDis_2026_Sample5

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample5；patient: nan。

组织spots 3125，有效Epi 2455，内部Epi 1983，冻结P95 100。精确前沿Epi 1186，前沿内部Epi 205，R12标注Epi 2120，未标注 335。

## 候选距离原始Pearson

- 核心面板较低r：Matrix_CAF / BUDDEN_S1P_RESPONSE_18: r=-0.661, n=95, 0–2563 μm, spatial CI [-0.751, -0.517], P=3.26e-13, across-section q=6.56e-11。
- 核心面板较高r：LAMP3 DC / BUDDEN_CHOLESTEROL_SYNTHESIS_27: r=0.655, n=100, 0–964 μm, spatial CI [0.314, 0.804], P=1.51e-13, across-section q=3.28e-11。

核心面板最大绝对共同范围差异：SMC / HALLMARK_E2F_TARGETS，r_P95=-0.232 (n=95)，r_nonP95=0.247 (n=1847)，delta=-0.480，空间CI [-0.646, -0.234]，空间P=3.86e-06，跨片空间q=0.0011，共同范围0–2587 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.255；较低 BUDDEN_S1P_RESPONSE_18, β=-0.066。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.245；较低 BUDDEN_S1P_RESPONSE_18, β=-0.091。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.282；较低 BUDDEN_S1P_RESPONSE_18, β=-0.099。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.278；较低 BUDDEN_S1P_RESPONSE_18, β=-0.119。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.200；较低 BUDDEN_S1P_RESPONSE_18, β=-0.113。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.001, CI [-0.171, 0.170], n=1028, 新85上下文族q=0.997。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.132, CI [-0.092, 0.355], n=1028, 新85上下文族q=0.647。
- C2_P95_enrichment / P95_probability: β=0.018, CI [-0.016, 0.052], n=1028, 新85上下文族q=0.403。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.299, CI [-0.561, -0.036], n=205, 新85上下文族q=0.121。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.155, CI [0.043, 0.267], n=1186, 新85上下文族q=0.116。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 3：-22.97个百分点，CI [-46.17, -0.51]，标注P95 n=10、非P95 n=185，跨片q=0.263。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 18 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
