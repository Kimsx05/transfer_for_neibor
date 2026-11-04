# Inhousedata_T4469479

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4469479；patient: nan。

组织spots 3917，有效Epi 1846，内部Epi 1583，冻结P95 80。精确前沿Epi 1040，前沿内部Epi 199，R12标注Epi 1619，未标注 227。

## 候选距离原始Pearson

- 核心面板较低r：iCAF / HALLMARK_MTORC1_SIGNALING__NO_MP10: r=-0.648, n=80, 0–819 μm, spatial CI [-0.770, -0.504], P=8.27e-11, across-section q=9.59e-09。
- 核心面板较高r：Naive/Memory B / HALLMARK_MTORC1_SIGNALING__NO_MP10: r=0.677, n=80, 0–781 μm, spatial CI [0.579, 0.764], P=5.48e-12, across-section q=9.61e-10。

核心面板最大绝对共同范围差异：cDC1 / HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY，r_P95=0.273 (n=80)，r_nonP95=-0.242 (n=1502)，delta=0.515，空间CI [0.224, 0.713]，空间P=3.3e-05，跨片空间q=0.00399，共同范围0–1418 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.610；较低 HALLMARK_OXIDATIVE_PHOSPHORYLATION, β=-0.098。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.535；较低 HALLMARK_E2F_TARGETS, β=-0.090。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.501；较低 HALLMARK_E2F_TARGETS, β=-0.065。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.436；较低 HALLMARK_E2F_TARGETS, β=-0.053。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.292；较低 HALLMARK_E2F_TARGETS, β=-0.029。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.129, CI [-0.302, 0.045], n=975, 新85上下文族q=0.357。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.071, CI [-0.311, 0.169], n=975, 新85上下文族q=0.766。
- C2_P95_enrichment / P95_probability: β=-0.012, CI [-0.046, 0.022], n=975, 新85上下文族q=0.859。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.195, CI [-0.694, 0.305], n=199, 新85上下文族q=0.859。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.030, CI [-0.145, 0.084], n=1040, 新85上下文族q=0.859。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 2：21.20个百分点，CI [nan, nan]，标注P95 n=4、非P95 n=184，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
