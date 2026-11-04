# CanDis_2026_Sample13

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample13；patient: nan。

组织spots 4835，有效Epi 3430，内部Epi 3151，冻结P95 158。精确前沿Epi 603，前沿内部Epi 141，R12标注Epi 3182，未标注 248。

## 候选距离原始Pearson

- 核心面板较低r：GC B / HALLMARK_E2F_TARGETS: r=-0.569, n=158, 0–624 μm, spatial CI [-0.687, -0.404], P=6.2e-15, across-section q=1.89e-12。
- 核心面板较高r：Neutrophil / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.557, n=158, 0–900 μm, spatial CI [0.362, 0.716], P=2.81e-14, across-section q=7.42e-12。

核心面板最大绝对共同范围差异：Naive/Memory B / HALLMARK_E2F_TARGETS，r_P95=0.120 (n=158)，r_nonP95=-0.331 (n=2805)，delta=0.451，空间CI [0.236, 0.701]，空间P=0.000108，跨片空间q=0.0082，共同范围0–1905 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.321；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.097。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 HALLMARK_MTORC1_SIGNALING__NO_MP10, β=0.376；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.166。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.403；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.181。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 HALLMARK_E2F_TARGETS, β=0.412；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.290。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.438；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.382。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.044, CI [-0.226, 0.314], n=500, 新85上下文族q=0.901。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.016, CI [-0.183, 0.214], n=500, 新85上下文族q=0.98。
- C2_P95_enrichment / P95_probability: β=-0.002, CI [-0.028, 0.023], n=500, 新85上下文族q=0.96。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.264, CI [-0.546, 0.017], n=141, 新85上下文族q=0.192。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.136, CI [-0.115, 0.388], n=603, 新85上下文族q=0.275。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：24.81个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=133，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 19 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
