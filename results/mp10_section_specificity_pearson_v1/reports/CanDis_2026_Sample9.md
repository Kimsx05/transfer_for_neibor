# CanDis_2026_Sample9

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample9；patient: nan。

组织spots 2893，有效Epi 2310，内部Epi 2204，冻结P95 111。精确前沿Epi 1588，前沿内部Epi 425，R12标注Epi 2130，未标注 180。

## 候选距离原始Pearson

- 核心面板较低r：Tfh/CXCL13-like / PURAM_PEMT_TOP30: r=-0.433, n=111, 0–1100 μm, spatial CI [-0.561, -0.273], P=2.1e-06, across-section q=6.06e-05。
- 核心面板较高r：Neutrophil / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.398, n=111, 0–1249 μm, spatial CI [0.188, 0.580], P=1.5e-05, across-section q=0.000314。

核心面板最大绝对共同范围差异：Capillary EC / BUDDEN_S1P_RESPONSE_18，r_P95=0.222 (n=110)，r_nonP95=-0.144 (n=2089)，delta=0.366，空间CI [0.197, 0.551]，空间P=0.000127，跨片空间q=0.00915，共同范围0–1473 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.324；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.088。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.159；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.059。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.124；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.064。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18, β=0.059；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.060。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.039；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.090。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.005, CI [-0.103, 0.114], n=1512, 新85上下文族q=0.96。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.193, CI [0.068, 0.319], n=1512, 新85上下文族q=0.0191。
- C2_P95_enrichment / P95_probability: β=0.015, CI [-0.010, 0.039], n=1512, 新85上下文族q=0.297。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.163, CI [-0.484, 0.158], n=425, 新85上下文族q=0.352。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.057, CI [-0.034, 0.149], n=1588, 新85上下文族q=0.284。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：9.27个百分点，CI [-14.68, 32.23]，标注P95 n=20、非P95 n=383，跨片q=0.68。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 21 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
