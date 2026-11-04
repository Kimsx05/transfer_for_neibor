# CanDis_2026_Sample11

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample11；patient: nan。

组织spots 577，有效Epi 247，内部Epi 93，冻结P95 5。精确前沿Epi 77，前沿内部Epi 9，R12标注Epi 200，未标注 47。

## 候选距离原始Pearson

此片没有达到n≥10的核心面板相关展示；全部可计算n≥2结果仍保留在完整长表，不因覆盖不足删除。

核心面板最大绝对共同范围差异：Tfh/CXCL13-like / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=-0.969 (n=4)，r_nonP95=0.561 (n=27)，delta=-1.530，空间CI [nan, nan]，空间P=nan，跨片空间q=nan，共同范围500–721 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.444；较低 BUDDEN_S1P_RESPONSE_18, β=-0.141。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=0.206；较低 BUDDEN_S1P_RESPONSE_18, β=-0.071。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=0.168；较低 PURAM_PEMT_TOP30, β=-0.088。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=0.127；较低 PURAM_PEMT_TOP30, β=-0.101。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_S1P_RESPONSE_18, β=0.136；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.032。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=nan, CI [nan, nan], n=9, 新85上下文族q=nan。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=nan, CI [nan, nan], n=9, 新85上下文族q=nan。
- C2_P95_enrichment / P95_probability: β=nan, CI [nan, nan], n=9, 新85上下文族q=nan。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=nan, CI [nan, nan], n=9, 新85上下文族q=nan。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.252, CI [-0.074, 0.578], n=77, 新85上下文族q=0.213。

## 限制与入口

29候选中 23 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 本片P95少于30，相关取值仍报告，但其稳定性有限，不能据较大|r|认定样本特异表型。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
