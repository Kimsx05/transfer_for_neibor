# Inhousedata_T4463209

来源: tumor_study_specimen_not_region_diagnosis；analysis unit: Inhousedata__UNVERIFIED::anonymous_unit::Inhousedata_T4463209；patient: nan。

组织spots 4441，有效Epi 1829，内部Epi 1496，冻结P95 75。精确前沿Epi 891，前沿内部Epi 174，R12标注Epi 1637，未标注 192。

## 候选距离原始Pearson

- 核心面板较低r：Mast / PURAM_PEMT_TOP30: r=-0.452, n=75, 0–1308 μm, spatial CI [-0.640, -0.219], P=4.73e-05, across-section q=0.000814。
- 核心面板较高r：Venous EC / HALLMARK_GLYCOLYSIS: r=0.616, n=75, 520–3148 μm, spatial CI [0.425, 0.766], P=4.11e-09, across-section q=2.78e-07。

核心面板最大绝对共同范围差异：Mast / PURAM_PEMT_TOP30，r_P95=-0.452 (n=75)，r_nonP95=-0.036 (n=1410)，delta=-0.416，空间CI [-0.610, -0.202]，空间P=0.00061，跨片空间q=0.0232，共同范围0–1308 μm。效应选择是探索性展示，不按显著性筛片。

## MP10中心及邻域

- R0: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.429；较低 HALLMARK_E2F_TARGETS, β=-0.160。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R1: 核心面板较高 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=0.188；较低 HALLMARK_E2F_TARGETS, β=-0.129。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R2: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.157；较低 HALLMARK_E2F_TARGETS, β=-0.094。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R3_R5: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.147；较低 HALLMARK_E2F_TARGETS, β=-0.056。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。
- R6_R10: 核心面板较高 BUDDEN_TNF_EMT_RESPONSE_15, β=0.194；较低 BUDDEN_CHOLESTEROL_SYNTHESIS_27, β=-0.020。这是UCell环平均程序的调整标准化系数，不是候选距离Pearson。

## 前沿及niche

- A_front_vs_internal / FLORES_CRISPR_PRO_INVASION_10: β=-0.093, CI [-0.372, 0.185], n=672, 新85上下文族q=0.684。
- A_front_vs_internal / PURAM_PEMT_TOP30: β=-0.175, CI [-0.346, -0.005], n=672, 新85上下文族q=0.148。
- C2_P95_enrichment / P95_probability: β=-0.011, CI [-0.038, 0.016], n=672, 新85上下文族q=0.53。
- C_active_front_MP10 / EpiMP10_Top20_AUCell: β=0.431, CI [-0.153, 1.015], n=174, 新85上下文族q=0.366。
- C1_all_Epi_continuous / EpiMP10_Top20_AUCell: β=-0.079, CI [-0.240, 0.082], n=891, 新85上下文族q=0.529。

仅前沿带[-200,0] μm，最大绝对R12比例差为niche 2：-3.55个百分点，CI [nan, nan]，标注P95 n=3、非P95 n=169，跨片q=nan。所有九类和未标注覆盖同时保留；这里没有重新定义Top5%。

## 限制与入口

29候选中 24 类的目标核心posterior mean丰度中位数<1；这是定位支持描述，不是剔除标准。P95只是片内相对排名。 精确前沿NA不替代为距离带中点。普通Pearson P/q不处理空间自相关；500 μm块区间单独列。缺区间不是缺效应；niche未标注不是该niche生物学为零。

三页主报告见同名PDF，完整范围散点见figures/sections，全部数值见该片pearson_q/delta_q和紧凑spot源数据。
