# CanDis_2026_Sample10

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample10；patient: nan。

组织spots 2360，有效Epi 2111，内部Epi 1726，冻结P95 87。精确前沿Epi 321，前沿内部Epi 105，R12标注Epi 1500，未标注 611。

## 候选距离原始Pearson

- 核心面板较低r：Classical Mono / BUDDEN_S1P_RESPONSE_18: r=-0.448, n=86, 173–4331 μm, spatial CI [-0.610, -0.212], P=1.5e-05, across-section q=0.000315。
- 核心面板较高r：Naive/Memory B / HALLMARK_MTORC1_SIGNALING__NO_MP10: r=0.622, n=86, 0–1646 μm, spatial CI [0.452, 0.755], P=1.63e-10, across-section q=1.77e-08。

核心面板最大绝对共同范围差异：Classical Mono / HALLMARK_E2F_TARGETS，r_P95=0.174 (n=86)，r_nonP95=-0.296 (n=1580)，delta=0.470，空间CI [0.195, 0.684]，空间P=0.000181，跨片空间q=0.0113，共同范围173–4331 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.507；较低 HALLMARK_E2F_TARGETS, β=-0.091。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.338；较低 HALLMARK_E2F_TARGETS, β=-0.135。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.291；较低 HALLMARK_E2F_TARGETS, β=-0.126。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.233；较低 HALLMARK_E2F_TARGETS, β=-0.170。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.115；较低 BUDDEN_TG_CATABOLISM_24, β=-0.085。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.054, CI [-0.360, 0.251], n=265, 新85上下文族q=0.834。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.336, CI [-0.001, 0.673], n=265, 新85上下文族q=0.169。
- C2_P95_enrichment / P95_probability: β=0.055, CI [0.008, 0.101], n=265, 新85上下文族q=0.102。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.187, CI [-0.567, 0.194], n=105, 新85上下文族q=0.386。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.038, CI [-0.240, 0.164], n=321, 新85上下文族q=0.254。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：35.91个百分点，CI [nan, nan]，标注P95 n=6、非P95 n=97，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 25 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
