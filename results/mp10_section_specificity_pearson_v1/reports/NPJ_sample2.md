# NPJ_sample2

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample2；patient: nan。

组织spots 4506，有效Epi 3217，内部Epi 3006，冻结P95 151。精确前沿Epi 1026，前沿内部Epi 172，R12标注Epi 2832，未标注 385。

## 候选距离原始Pearson

- 核心面板较低r：SPP1 TAM / HALLMARK_GLYCOLYSIS: r=-0.573, n=151, 0–1249 μm, spatial CI [-0.669, -0.414], P=1.46e-14, across-section q=4.19e-12。
- 核心面板较高r：Arterial EC / HALLMARK_GLYCOLYSIS: r=0.585, n=151, 100–3081 μm, spatial CI [0.425, 0.683], P=3.22e-15, across-section q=1.08e-12。

核心面板最大绝对共同范围差异：rCAF / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=0.434 (n=151)，r_nonP95=-0.204 (n=2828)，delta=0.638，空间CI [0.411, 0.822]，空间P=4.66e-07，跨片空间q=0.000285，共同范围100–3372 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.407；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.178。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.351；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.175。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.288；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.184。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.195；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.157。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 FLORES_CRISPR_PRO_INVASION_10, β=0.137；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.097。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.249, CI [-0.403, -0.096], n=970, 新85上下文族q=0.00805。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.389, CI [0.185, 0.594], n=970, 新85上下文族q=0.00234。
- C2_P95_enrichment / P95_probability: β=-0.009, CI [-0.025, 0.006], n=970, 新85上下文族q=0.286。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.020, CI [-0.326, 0.287], n=172, 新85上下文族q=0.896。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.450, CI [-0.567, -0.333], n=1026, 新85上下文族q=1.42e-09。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 2：35.00个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=140，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 25 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
