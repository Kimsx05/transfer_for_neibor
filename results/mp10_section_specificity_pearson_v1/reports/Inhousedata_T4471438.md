# Inhousedata_T4471438

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4471438；patient: nan。

组织spots 4373，有效Epi 2256，内部Epi 2024，冻结P95 102。精确前沿Epi 801，前沿内部Epi 147，R12标注Epi 2093，未标注 163。

## 候选距离原始Pearson

- 核心面板较低r：cDC1 / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.399, n=102, 0–1300 μm, spatial CI [-0.622, -0.147], P=3.19e-05, across-section q=0.000587。
- 核心面板较高r：CD8 Exhausted / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.596, n=102, 0–1992 μm, spatial CI [0.422, 0.753], P=4.01e-11, across-section q=4.77e-09。

核心面板最大绝对共同范围差异：Pericyte / HALLMARK_E2F_TARGETS，r_P95=-0.295 (n=102)，r_nonP95=0.211 (n=1807)，delta=-0.506，空间CI [-0.784, -0.225]，空间P=0.00168，跨片空间q=0.0422，共同范围100–1212 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.480；较低 PURAM_PEMT_TOP30, β=-0.221。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.402；较低 BUDDEN_TG_CATABOLISM_24, β=-0.096。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.317；较低 HALLMARK_E2F_TARGETS, β=-0.117。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.317；较低 BUDDEN_TG_CATABOLISM_24, β=-0.216。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.368；较低 BUDDEN_TG_CATABOLISM_24, β=-0.280。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.222, CI [-0.395, -0.048], n=749, 新85上下文族q=0.0767。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.253, CI [-0.417, -0.089], n=749, 新85上下文族q=0.0349。
- C2_P95_enrichment / P95_probability: β=0.034, CI [-0.022, 0.091], n=749, 新85上下文族q=0.373。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.058, CI [-0.306, 0.189], n=147, 新85上下文族q=0.629。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.057, CI [-0.055, 0.169], n=801, 新85上下文族q=0.19。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：12.10个百分点，CI [nan, nan]，标注P95 n=12、非P95 n=124，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 26 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
