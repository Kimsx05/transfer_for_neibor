# CanDis_2026_Sample1

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample1；patient: nan。

组织spots 4823，有效Epi 4486，内部Epi 3636，冻结P95 183。精确前沿Epi 2759，前沿内部Epi 375，R12标注Epi 4227，未标注 259。

## 候选距离原始Pearson

- 核心面板较低r：Treg / PURAM_PEMT_TOP30: r=-0.681, n=183, 0–2615 μm, spatial CI [-0.755, -0.582], P=2.52e-26, across-section q=2.3e-22。
- 核心面板较高r：SMC / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.631, n=183, 0–3857 μm, spatial CI [0.501, 0.729], P=1.08e-21, across-section q=2.47e-18。

核心面板最大绝对共同范围差异：cDC1 / HALLMARK_OXIDATIVE_PHOSPHORYLATION，r_P95=0.423 (n=183)，r_nonP95=0.124 (n=3444)，delta=0.299，空间CI [0.146, 0.410]，空间P=4.65e-05，跨片空间q=0.00485，共同范围0–1353 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.166；较低 HALLMARK_G2M_CHECKPOINT, β=-0.087。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.090；较低 HALLMARK_E2F_TARGETS, β=-0.099。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.099；较低 HALLMARK_E2F_TARGETS, β=-0.109。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.092；较低 HALLMARK_G2M_CHECKPOINT, β=-0.113。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.100；较低 HALLMARK_G2M_CHECKPOINT, β=-0.052。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.051, CI [-0.159, 0.056], n=2225, 新85上下文族q=0.564。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.142, CI [0.002, 0.282], n=2225, 新85上下文族q=0.215。
- C2_P95_enrichment / P95_probability: β=0.054, CI [0.012, 0.095], n=2225, 新85上下文族q=0.0167。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.067, CI [-0.216, 0.349], n=375, 新85上下文族q=0.639。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.141, CI [0.071, 0.212], n=2759, 新85上下文族q=0.000525。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：8.77个百分点，CI [-9.14, 26.57]，标注P95 n=30、非P95 n=333，跨片q=0.632。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
