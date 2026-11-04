# CanDis_2026_Sample7

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample7；patient: nan。

组织spots 2295，有效Epi 1858，内部Epi 1366，冻结P95 69。精确前沿Epi 1033，前沿内部Epi 325，R12标注Epi 1576，未标注 282。

## 候选距离原始Pearson

- 核心面板较低r：SPP1 TAM / BUDDEN_S1P_RESPONSE_18: r=-0.687, n=68, 0–1836 μm, spatial CI [-0.802, -0.376], P=9.5e-11, across-section q=1.02e-08。
- 核心面板较高r：CD8 Exhausted / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.614, n=50, 173–2081 μm, spatial CI [0.261, 0.766], P=2.13e-06, across-section q=6.12e-05。

核心面板最大绝对共同范围差异：Classical Mono / HALLMARK_GLYCOLYSIS，r_P95=0.204 (n=68)，r_nonP95=-0.167 (n=1273)，delta=0.371，空间CI [0.182, 0.547]，空间P=0.000165，跨片空间q=0.0108，共同范围0–1905 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.268；较低 HALLMARK_G2M_CHECKPOINT, β=-0.092。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.082；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.068。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.086；较低 HALLMARK_G2M_CHECKPOINT, β=-0.059。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.072；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.013。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 PURAM_PEMT_TOP30, β=0.082；较低 HALLMARK_E2F_TARGETS, β=-0.025。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.101, CI [-0.301, 0.098], n=644, 新85上下文族q=0.452。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.004, CI [-0.194, 0.187], n=644, 新85上下文族q=0.979。
- C2_P95_enrichment / P95_probability: β=0.018, CI [-0.008, 0.044], n=644, 新85上下文族q=0.645。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.152, CI [-0.129, 0.432], n=325, 新85上下文族q=0.645。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.019, CI [-0.074, 0.113], n=1033, 新85上下文族q=0.775。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 5：-7.97个百分点，CI [-16.44, -1.11]，标注P95 n=15、非P95 n=301，跨片q=0.296。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 26 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
