# CanDis_2026_Sample8

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: CanDis__UNVERIFIED::anonymous_unit::CanDis_2026_Sample8；patient: nan。

组织spots 11958，有效Epi 10711，内部Epi 10037，冻结P95 502。精确前沿Epi 2755，前沿内部Epi 320，R12标注Epi 10188，未标注 523。

## 候选距离原始Pearson

- 核心面板较低r：SPP1 TAM / PURAM_PEMT_TOP30: r=-0.331, n=502, 0–2007 μm, spatial CI [-0.426, -0.232], P=2.56e-14, across-section q=6.79e-12。
- 核心面板较高r：CD4 Naive/Memory-like / HALLMARK_OXIDATIVE_PHOSPHORYLATION: r=0.421, n=502, 0–1664 μm, spatial CI [0.327, 0.501], P=5.69e-23, across-section q=1.73e-19。

核心面板最大绝对共同范围差异：pDC / HALLMARK_E2F_TARGETS，r_P95=0.286 (n=502)，r_nonP95=0.078 (n=9523)，delta=0.207，空间CI [0.111, 0.292]，空间P=1.67e-05，跨片空间q=0.00261，共同范围0–1552 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.286；较低 PURAM_PEMT_TOP30, β=-0.176。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.241；较低 PURAM_PEMT_TOP30, β=-0.181。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.249；较低 PURAM_PEMT_TOP30, β=-0.170。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.265；较低 PURAM_PEMT_TOP30, β=-0.154。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.238；较低 FLORES_CRISPR_PRO_INVASION_10, β=-0.154。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.065, CI [-0.201, 0.071], n=2507, 新85上下文族q=0.454。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.188, CI [-0.424, 0.048], n=2507, 新85上下文族q=0.207。
- C2_P95_enrichment / P95_probability: β=0.019, CI [-0.015, 0.053], n=2507, 新85上下文族q=0.379。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=-0.265, CI [-0.580, 0.050], n=320, 新85上下文族q=0.29。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=0.069, CI [-0.003, 0.140], n=2755, 新85上下文族q=0.00081。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 1：7.94个百分点，CI [-16.67, 33.11]，标注P95 n=20、非P95 n=293，跨片q=0.78。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 25 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
