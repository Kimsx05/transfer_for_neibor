# NPJ_sample5

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: NPJ__UNVERIFIED::anonymous_unit::NPJ_sample5；patient: nan。

组织spots 2815，有效Epi 2186，内部Epi 2178，冻结P95 109。精确前沿Epi 656，前沿内部Epi 122，R12标注Epi 2006，未标注 180。

## 候选距离原始Pearson

- 核心面板较低r：GC B / HALLMARK_E2F_TARGETS: r=-0.508, n=108, 0–819 μm, spatial CI [-0.670, -0.325], P=1.97e-08, across-section q=1.09e-06。
- 核心面板较高r：iCAF / HALLMARK_GLYCOLYSIS: r=0.450, n=108, 346–2946 μm, spatial CI [0.191, 0.642], P=1.02e-06, across-section q=3.3e-05。

核心面板最大绝对共同范围差异：CD8 GZMK/Effector-memory-like / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=-0.050 (n=108)，r_nonP95=-0.446 (n=1998)，delta=0.396，空间CI [0.162, 0.667]，空间P=0.00131，跨片空间q=0.0364，共同范围100–2946 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.302；较低 PURAM_PEMT_TOP30, β=-0.015。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 HALLMARK_MTORC1_SIGNALING__NO_MP10, β=0.282；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.080。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 HALLMARK_MTORC1_SIGNALING__NO_MP10, β=0.236；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.072。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.253；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.085。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.177；较低 HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY, β=-0.063。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.100, CI [-0.525, 0.324], n=652, 新85上下文族q=0.805。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.321, CI [-0.799, 0.156], n=652, 新85上下文族q=0.347。
- C2_P95_enrichment / P95_probability: β=0.011, CI [-0.045, 0.067], n=652, 新85上下文族q=0.842。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.131, CI [-0.356, 0.618], n=122, 新85上下文族q=0.842。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.283, CI [-0.449, -0.117], n=656, 新85上下文族q=0.0079。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 5：84.54个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=97，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 25 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
