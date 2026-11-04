# CanDis_2026_Sample3

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample3；patient: nan。

组织spots 4818，有效Epi 4553，内部Epi 3805，冻结P95 191。精确前沿Epi 1607，前沿内部Epi 493，R12标注Epi 4239，未标注 314。

## 候选距离原始Pearson

- 核心面板较低r：Neutrophil / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.374, n=191, 0–721 μm, spatial CI [-0.493, -0.239], P=9.95e-08, across-section q=4.5e-06。
- 核心面板较高r：Capillary EC / BUDDEN_TNF_EMT_RESPONSE_15: r=0.418, n=191, 0–1637 μm, spatial CI [0.274, 0.524], P=1.81e-09, across-section q=1.37e-07。

核心面板最大绝对共同范围差异：pDC / HALLMARK_GLYCOLYSIS，r_P95=0.238 (n=191)，r_nonP95=-0.128 (n=3579)，delta=0.366，空间CI [0.198, 0.521]，空间P=2.84e-05，跨片空间q=0.00368，共同范围0–872 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.341；较低 PURAM_PEMT_TOP30, β=-0.194。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.287；较低 PURAM_PEMT_TOP30, β=-0.197。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.239；较低 PURAM_PEMT_TOP30, β=-0.169。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.202；较低 BUDDEN_TG_CATABOLISM_24, β=-0.142。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.168；较低 PURAM_PEMT_TOP30, β=-0.100。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.114, CI [-0.276, 0.048], n=1067, 新85上下文族q=0.345。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.038, CI [-0.135, 0.212], n=1067, 新85上下文族q=0.761。
- C2_P95_enrichment / P95_probability: β=-0.030, CI [-0.067, 0.007], n=1067, 新85上下文族q=0.164。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.509, CI [-0.767, -0.252], n=493, 新85上下文族q=0.00168。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.084, CI [-0.291, 0.124], n=1607, 新85上下文族q=0.265。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：6.71个百分点，CI [-4.85, 15.51]，标注P95 n=18、非P95 n=473，跨片q=0.531。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 20 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
