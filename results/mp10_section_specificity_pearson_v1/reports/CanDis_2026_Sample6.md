# CanDis_2026_Sample6

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample6；patient: nan。

组织spots 1985，有效Epi 1753，内部Epi 1138，冻结P95 57。精确前沿Epi 729，前沿内部Epi 131，R12标注Epi 1407，未标注 346。

## 候选距离原始Pearson

- 核心面板较低r：vCAF / PURAM_PEMT_TOP30: r=-0.675, n=31, 0–1044 μm, spatial CI [-0.826, -0.450], P=3.14e-05, across-section q=0.000579。
- 核心面板较高r：Lymphatic EC / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.681, n=28, 0–721 μm, spatial CI [0.418, 0.858], P=6.6e-05, across-section q=0.00107。

核心面板最大绝对共同范围差异：C1QC TAM / PURAM_PEMT_TOP30，r_P95=0.148 (n=41)，r_nonP95=-0.573 (n=1017)，delta=0.721，空间CI [0.458, 0.928]，空间P=1.09e-09，跨片空间q=6.09e-06，共同范围0–1852 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.318；较低 BUDDEN_S1P_RESPONSE_18, β=-0.087。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.213；较低 PURAM_PEMT_TOP30, β=-0.099。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.211；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.102。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.207；较低 BUDDEN_TNF_EMT_RESPONSE_15, β=-0.145。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.052；较低 BUDDEN_S1P_RESPONSE_18, β=-0.029。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.041, CI [-0.347, 0.265], n=518, 新85上下文族q=0.86。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.406, CI [-0.574, -0.238], n=518, 新85上下文族q=0.00307。
- C2_P95_enrichment / P95_probability: β=0.025, CI [-0.016, 0.066], n=518, 新85上下文族q=0.692。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.075, CI [-0.586, 0.435], n=131, 新85上下文族q=0.859。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.103, CI [-0.280, 0.074], n=729, 新85上下文族q=0.0228。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：39.02个百分点，CI [nan, nan]，标注P95 n=3、非P95 n=123，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 27 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
