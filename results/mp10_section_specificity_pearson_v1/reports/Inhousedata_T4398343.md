# Inhousedata_T4398343

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4398343；patient: nan。

组织spots 4691，有效Epi 1517，内部Epi 1483，冻结P95 75。精确前沿Epi 695，前沿内部Epi 224，R12标注Epi 1411，未标注 106。

## 候选距离原始Pearson

- 核心面板较低r：cDC2 / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.516, n=75, 100–1015 μm, spatial CI [-0.674, -0.361], P=2.2e-06, across-section q=6.32e-05。
- 核心面板较高r：pDC / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.599, n=75, 100–1015 μm, spatial CI [0.437, 0.726], P=1.36e-08, across-section q=7.97e-07。

核心面板最大绝对共同范围差异：C1QC TAM / PURAM_PEMT_TOP30，r_P95=-0.402 (n=75)，r_nonP95=0.087 (n=1177)，delta=-0.489，空间CI [-0.676, -0.257]，空间P=5.58e-05，跨片空间q=0.00546，共同范围100–900 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.462；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.098。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.355；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.110。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.295；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.167。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.211；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.196。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 PURAM_PEMT_TOP30, β=0.151；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.219。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.211, CI [-0.427, 0.004], n=675, 新85上下文族q=0.139。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.196, CI [-0.094, 0.486], n=675, 新85上下文族q=0.316。
- C2_P95_enrichment / P95_probability: β=0.002, CI [-0.015, 0.018], n=675, 新85上下文族q=0.968。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.016, CI [-0.274, 0.242], n=224, 新85上下文族q=0.968。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.104, CI [-0.050, 0.259], n=695, 新85上下文族q=0.0741。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 0：8.41个百分点，CI [nan, nan]，标注P95 n=4、非P95 n=211，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 19 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
