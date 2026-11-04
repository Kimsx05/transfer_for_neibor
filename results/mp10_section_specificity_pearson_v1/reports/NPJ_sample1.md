# NPJ_sample1

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample1；patient: nan。

组织spots 2076，有效Epi 996，内部Epi 589，冻结P95 30。精确前沿Epi 508，前沿内部Epi 161，R12标注Epi 766，未标注 230。

## 候选距离原始Pearson

- 核心面板较低r：CD4 Helper/Conventional / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.715, n=30, 520–1735 μm, spatial CI [-0.891, -0.373], P=9.08e-06, across-section q=0.000206。
- 核心面板较高r：CD8 GZMK/Effector-memory-like / BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18: r=0.855, n=30, 458–1682 μm, spatial CI [0.526, 0.920], P=1.86e-09, across-section q=1.43e-07。

核心面板最大绝对共同范围差异：Treg / HALLMARK_MTORC1_SIGNALING__NO_MP10，r_P95=-0.478 (n=30)，r_nonP95=0.425 (n=524)，delta=-0.903，空间CI [-1.232, -0.066]，空间P=0.00103，跨片空间q=0.0288，共同范围100–1044 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.578；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.146。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.411；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.176。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.336；较低 BUDDEN_S1P_RESPONSE_18, β=-0.171。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.181；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.099。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.158；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.071。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.030, CI [-0.301, 0.241], n=280, 新85上下文族q=0.936。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.113, CI [-0.134, 0.360], n=280, 新85上下文族q=0.651。
- C2_P95_enrichment / P95_probability: β=-0.004, CI [-0.039, 0.031], n=280, 新85上下文族q=0.945。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.020, CI [-0.605, 0.565], n=161, 新85上下文族q=0.945。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.194, CI [-0.302, -0.085], n=508, 新85上下文族q=0.00768。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：22.79个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=136，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 22 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
