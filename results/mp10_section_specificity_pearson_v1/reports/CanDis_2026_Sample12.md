# CanDis_2026_Sample12

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample12；patient: nan。

组织spots 3508，有效Epi 405，内部Epi 390，冻结P95 20。精确前沿Epi 189，前沿内部Epi 66，R12标注Epi 370，未标注 35。

## 候选距离原始Pearson

- 核心面板较低r：pDC / BUDDEN_TNF_EMT_RESPONSE_15: r=-0.808, n=20, 0–854 μm, spatial CI [-0.950, -0.459], P=1.62e-05, across-section q=0.000336。
- 核心面板较高r：Naive/Memory B / BUDDEN_TNF_EMT_RESPONSE_15: r=0.587, n=20, 100–1277 μm, spatial CI [0.053, 0.856], P=0.0065, across-section q=0.0395。

核心面板最大绝对共同范围差异：Matrix_CAF / BUDDEN_TNF_EMT_RESPONSE_15，r_P95=-0.547 (n=20)，r_nonP95=0.402 (n=291)，delta=-0.949，空间CI [-1.292, -0.547]，空间P=2.55e-05，跨片空间q=0.00345，共同范围173–1082 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.934；较低 HALLMARK_G2M_CHECKPOINT, β=-0.376。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.863；较低 HALLMARK_G2M_CHECKPOINT, β=-0.431。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.708；较低 HALLMARK_G2M_CHECKPOINT, β=-0.390。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.544；较低 HALLMARK_G2M_CHECKPOINT, β=-0.451。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.504；较低 HALLMARK_G2M_CHECKPOINT, β=-0.458。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=0.072, CI [-0.237, 0.380], n=185, 新85上下文族q=0.734。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=0.080, CI [-0.310, 0.469], n=185, 新85上下文族q=0.769。
- C2_P95_enrichment / P95_probability: β=-0.011, CI [-0.035, 0.013], n=185, 新85上下文族q=0.395。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.301, CI [-0.688, 0.085], n=66, 新85上下文族q=0.207。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.585, CI [-0.896, -0.274], n=189, 新85上下文族q=0.000141。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 2：57.38个百分点，CI [nan, nan]，标注P95 n=1、非P95 n=61，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 本片P95少于30，相关取值仍报告，但其稳定性有限，不能据较大|r|认定样本特异表型。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
