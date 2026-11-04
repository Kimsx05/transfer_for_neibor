# NPJ_sample7

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample7；patient: nan。

组织spots 3386，有效Epi 2382，内部Epi 1957，冻结P95 99。精确前沿Epi 1178，前沿内部Epi 237，R12标注Epi 2129，未标注 253。

## 候选距离原始Pearson

- 核心面板较低r：Treg / PURAM_PEMT_TOP30: r=-0.369, n=99, 100–1127 μm, spatial CI [-0.548, -0.137], P=0.00017, across-section q=0.00231。
- 核心面板较高r：SMC / HALLMARK_GLYCOLYSIS: r=0.520, n=99, 1311–4115 μm, spatial CI [0.284, 0.693], P=3.52e-08, across-section q=1.82e-06。

核心面板最大绝对共同范围差异：iCAF / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=-0.319 (n=99)，r_nonP95=0.272 (n=1390)，delta=-0.591，空间CI [-0.810, -0.312]，空间P=2.85e-05，跨片空间q=0.00368，共同范围2600–5226 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.238；较低 HALLMARK_E2F_TARGETS, β=-0.117。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.266；较低 HALLMARK_G2M_CHECKPOINT, β=-0.141。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.271；较低 HALLMARK_G2M_CHECKPOINT, β=-0.120。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.281；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.176。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.289；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.155。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.023, CI [-0.178, 0.224], n=913, 新85上下文族q=0.906。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.058, CI [-0.252, 0.136], n=913, 新85上下文族q=0.763。
- C2_P95_enrichment / P95_probability: β=-0.018, CI [-0.051, 0.014], n=913, 新85上下文族q=0.412。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.219, CI [-0.581, 0.142], n=237, 新85上下文族q=0.412。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.071, CI [-0.261, 0.120], n=1178, 新85上下文族q=0.779。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 7：85.71个百分点，CI [nan, nan]，标注P95 n=4、非P95 n=224，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 25 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
