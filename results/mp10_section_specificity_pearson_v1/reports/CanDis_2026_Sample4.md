# CanDis_2026_Sample4

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample4；patient: nan。

组织spots 3118，有效Epi 2997，内部Epi 2286，冻结P95 115。精确前沿Epi 1114，前沿内部Epi 315，R12标注Epi 2417，未标注 580。

## 候选距离原始Pearson

- 核心面板较低r：Neutrophil / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.561, n=113, 0–1114 μm, spatial CI [-0.705, -0.351], P=1.06e-10, across-section q=1.12e-08。
- 核心面板较高r：cDC1 / PURAM_PEMT_TOP30: r=0.548, n=113, 0–1386 μm, spatial CI [0.225, 0.730], P=3.32e-10, across-section q=3.09e-08。

核心面板最大绝对共同范围差异：SPP1 TAM / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=0.096 (n=114)，r_nonP95=-0.476 (n=1991)，delta=0.572，空间CI [0.342, 0.796]，空间P=8.56e-06，跨片空间q=0.00177，共同范围0–781 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.326；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.006。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_S1P_RESPONSE_18, β=0.260；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.014。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_S1P_RESPONSE_18, β=0.230；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.044。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_S1P_RESPONSE_18, β=0.189；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.084。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.115；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.084。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.015, CI [-0.162, 0.193], n=767, 新85上下文族q=0.921。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.094, CI [-0.246, 0.058], n=767, 新85上下文族q=0.477。
- C2_P95_enrichment / P95_probability: β=0.066, CI [0.016, 0.116], n=767, 新85上下文族q=0.0401。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.165, CI [-0.455, 0.125], n=315, 新85上下文族q=0.464。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.183, CI [-0.312, -0.053], n=1114, 新85上下文族q=0.0253。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 2：-28.47个百分点，CI [-42.55, -14.81]，标注P95 n=22、非P95 n=288，跨片q=0.0185。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 26 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
