# Inhousedata_T4384654

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4384654；patient: nan。

组织spots 4755，有效Epi 110，内部Epi 73，冻结P95 4。精确前沿Epi 26，前沿内部Epi 22，R12标注Epi 98，未标注 12。

## 候选距离原始Pearson

此片没有达到n≥10的核心面板相关展示；全部可计算n≥2结果仍保留在完整长表，不因覆盖不足删除。

核心面板最大绝对共同范围差异：CD4 Naive/Memory-like / BUDDEN_TG_CATABOLISM_24，r_P95=-0.999 (n=3)，r_nonP95=0.337 (n=55)，delta=-1.337，空间CI [nan, nan]，空间P=nan，跨片空间q=nan，共同范围100–265 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_S1P_RESPONSE_18, β=0.281；较低 HALLMARK_G2M_CHECKPOINT, β=-0.054。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 HALLMARK_G2M_CHECKPOINT, β=0.287；较低 BUDDEN_TG_CATABOLISM_24, β=-0.155。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=0.162；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=-0.141。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.139；较低 HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, β=-0.123。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 FLORES_CRISPR_PRO_INVASION_10, β=0.175；较低 PURAM_PEMT_TOP30, β=-0.111。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=nan, CI [nan, nan], n=22, 新85上下文族q=nan。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=nan, CI [nan, nan], n=22, 新85上下文族q=nan。
- C2_P95_enrichment / P95_probability: β=nan, CI [nan, nan], n=22, 新85上下文族q=nan。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.758, CI [-1.338, -0.178], n=22, 新85上下文族q=0.0232。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.220, CI [-1.021, 0.581], n=26, 新85上下文族q=0.00629。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：0.00个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=20，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 16 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 本片P95少于30，相关取值仍报告，但其稳定性有限，不能据较大|r|认定样本特异表型。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
