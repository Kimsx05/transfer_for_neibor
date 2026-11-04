# Inhousedata_T4465935

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4465935；patient: nan。

组织spots 4067，有效Epi 581，内部Epi 276，冻结P95 14。精确前沿Epi 409，前沿内部Epi 175，R12标注Epi 490，未标注 91。

## 候选距离原始Pearson

- 核心面板较低r：cDC1 / HALLMARK_UNFOLDED_PROTEIN_RESPONSE: r=-0.629, n=14, 0–200 μm, spatial CI [-0.916, 0.206], P=0.016, across-section q=0.077。
- 核心面板较高r：GC B / BUDDEN_TNF_EMT_RESPONSE_15: r=0.713, n=14, 0–346 μm, spatial CI [0.177, 0.952], P=0.0042, across-section q=0.0283。

核心面板最大绝对共同范围差异：cDC1 / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=0.685 (n=14)，r_nonP95=-0.052 (n=229)，delta=0.737，空间CI [0.484, 1.048]，空间P=1.18e-05，跨片空间q=0.00206，共同范围0–200 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.227；较低 HALLMARK_E2F_TARGETS, β=-0.031。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.095；较低 PURAM_PEMT_TOP30, β=-0.026。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.090；较低 BUDDEN_TG_CATABOLISM_24, β=-0.005。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 HALLMARK_GLYCOLYSIS, β=0.094；较低 PURAM_PEMT_TOP30, β=0.037。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TG_CATABOLISM_24, β=0.071；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.001。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.043, CI [-0.325, 0.412], n=257, 新85上下文族q=0.899。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.482, CI [0.171, 0.793], n=257, 新85上下文族q=0.0305。
- C2_P95_enrichment / P95_probability: β=0.004, CI [-0.060, 0.068], n=257, 新85上下文族q=0.99。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.092, CI [-0.348, 0.532], n=175, 新85上下文族q=0.99。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.092, CI [-0.247, 0.063], n=409, 新85上下文族q=0.99。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 5：24.79个百分点，CI [nan, nan]，标注P95 n=9、非P95 n=156，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 19 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 本片P95少于30，相关取值仍报告，但其稳定性有限，不能据较大|r|认定样本特异表型。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
