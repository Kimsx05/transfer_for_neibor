# MP10逐切片Pearson与样本差异：实际结果

本轮已登记并完成全部 63 张切片，265,002 个组织spot，85个程序(61主、24敏感性)和29个非Epi候选。所有切片均有三页报告与模块状态；没有依据效应方向、显著性或前沿覆盖筛片。MP10仍为冻结Top20 AUCell，新AddModuleScore不改变P95或pEMT-high front。

| source_context | sections | tissue_spots | Epi | internal_Epi | P95 | exact_front_Epi |
| --- | --- | --- | --- | --- | --- | --- |
| adjacent_normal_source | 3 | 5628 | 2309 | 0 | 0 | 0 |
| tumor_study_specimen_not_region_diagnosis | 60 | 259374 | 141152 | 112547 | 5664 | 58556 |

## 1. MP10高的Epi周围有哪些功能变化？

以下为与上一轮兼容的UCell中心/邻域模型，不是本轮候选距离Pearson。旧80项点级结果复用，新增5项同法补齐，并在85程序面板内重新进行族校正。效应先在既有analysis unit内等权，再在dataset间等权；95%区间为既有层级bootstrap，P来自dataset均值检验，因此二者不是同一检验的对偶。

| program | effect | low | high | q_context85 | n_units | n_datasets | positive_units |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HALLMARK_UV_RESPONSE_UP | 0.04894 | 0.01895 | 0.07787 | 0.01356 | 54 | 6 | 41 |
| HALLMARK_MTORC1_SIGNALING | 0.1585 | 0.1005 | 0.226 | 0.01472 | 54 | 6 | 53 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | 0.08151 | 0.02824 | 0.1433 | 0.09324 | 54 | 6 | 47 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | 0.04021 | 0.01154 | 0.07107 | 0.06958 | 54 | 6 | 38 |
| HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.03886 | 0.005935 | 0.07252 | 0.07637 | 54 | 6 | 37 |
| HALLMARK_E2F_TARGETS | 0.02552 | -0.01941 | 0.06792 | 0.2426 | 54 | 6 | 29 |
| HALLMARK_G2M_CHECKPOINT | 0.01852 | -0.02378 | 0.05808 | 0.3293 | 54 | 6 | 28 |
| PURAM_PEMT_TOP30 | -0.05084 | -0.08385 | -0.0251 | 0.01032 | 54 | 6 | 11 |
| FLORES_CRISPR_PRO_INVASION_10 | -0.009466 | -0.03095 | 0.01175 | 0.3178 | 54 | 6 | 24 |
| BUDDEN_CHOLESTEROL_SYNTHESIS_27 | 0.3564 | 0.2661 | 0.4402 | 0.005073 | 54 | 6 | 54 |
| BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | 0.1296 | 0.06695 | 0.2092 | 0.02982 | 54 | 6 | 49 |
| BUDDEN_S1P_RESPONSE_18 | -0.004411 | -0.04217 | 0.03837 | 0.777 | 54 | 6 | 22 |
| BUDDEN_TG_CATABOLISM_24 | 0.006189 | -0.03311 | 0.0467 | 0.7902 | 54 | 6 | 28 |
| BUDDEN_TNF_EMT_RESPONSE_15 | 0.04817 | -0.003575 | 0.09943 | 0.0643 | 54 | 6 | 35 |

非脂代谢中心信号中，UV_RESPONSE_UP的调整标准化β为0.04894 (95%区间0.01895–0.07787，85面板族q=0.01356；54单位、6个dataset中41单位为正)，效应较小但具有重复性。UPR、ROS和去MP10重叠mTORC1也偏正，但该层级q分别约0.070、0.076和0.093，暂作为支持不足的后续线索。Puram pEMT中心关系偏负(β=−0.05084，区间−0.08385至−0.02510，q=0.01032)，不能把MP10高直接等同于pEMT高。新增S1P和TG程序没有一致的总体中心升高。


原始脂代谢集合与MP10的共享基因会增加一致性；Budden胆固醇27与Top20重叠9个，去重叠18成员单列。mTORC1去重叠版本与原版也分别保留；OXPHOS原集合没有Top20重叠，因此不制造额外“去重叠”集合。程序名称不是代谢通量、细胞活性或暴露本身的直接测量。

| program | analysis | effect | low | high | q_context85 | n_units |
| --- | --- | --- | --- | --- | --- | --- |
| BUDDEN_S1P_RESPONSE_18 | B1_R0 | -0.004411 | -0.04217 | 0.03837 | 0.777 | 54 |
| BUDDEN_TNF_EMT_RESPONSE_15 | B1_R0 | 0.04817 | -0.003575 | 0.09943 | 0.0643 | 54 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | B1_R0 | 0.04021 | 0.01154 | 0.07107 | 0.06958 | 54 |
| PURAM_PEMT_TOP30 | B1_R0 | -0.05084 | -0.08385 | -0.0251 | 0.01032 | 54 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | B1_R0 | 0.08151 | 0.02824 | 0.1433 | 0.09324 | 54 |
| BUDDEN_S1P_RESPONSE_18 | B1_R1 | -0.0002347 | -0.03847 | 0.03954 | 0.9946 | 54 |
| BUDDEN_TNF_EMT_RESPONSE_15 | B1_R1 | 0.05176 | -0.002889 | 0.1069 | 0.1914 | 54 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | B1_R1 | 0.03962 | 0.005597 | 0.07388 | 0.2766 | 54 |
| PURAM_PEMT_TOP30 | B1_R1 | -0.03821 | -0.0717 | -0.008747 | 0.1321 | 54 |
| BUDDEN_S1P_RESPONSE_18 | B1_R2 | 0.002446 | -0.03281 | 0.03481 | 0.9027 | 54 |
| BUDDEN_TNF_EMT_RESPONSE_15 | B1_R2 | 0.04051 | -0.00792 | 0.08637 | 0.134 | 54 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | B1_R2 | 0.02701 | -0.003679 | 0.06293 | 0.4187 | 54 |
| PURAM_PEMT_TOP30 | B1_R2 | -0.02843 | -0.06608 | 0.001246 | 0.2888 | 54 |
| BUDDEN_S1P_RESPONSE_18 | B1_R3_R5 | 0.01149 | -0.01754 | 0.04148 | 0.3158 | 54 |
| BUDDEN_TNF_EMT_RESPONSE_15 | B1_R3_R5 | 0.04265 | -0.001603 | 0.0878 | 0.1137 | 54 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | B1_R3_R5 | 0.02033 | -0.009928 | 0.05566 | 0.561 | 54 |
| PURAM_PEMT_TOP30 | B1_R3_R5 | -0.01456 | -0.05756 | 0.02262 | 0.6832 | 54 |
| BUDDEN_S1P_RESPONSE_18 | B1_R6_R10 | 0.0158 | -0.02506 | 0.05983 | 0.6251 | 54 |
| BUDDEN_TNF_EMT_RESPONSE_15 | B1_R6_R10 | 0.04219 | -0.006248 | 0.09726 | 0.2978 | 54 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | B1_R6_R10 | 0.01643 | -0.01503 | 0.05093 | 0.5872 | 54 |
| PURAM_PEMT_TOP30 | B1_R6_R10 | -0.00154 | -0.05155 | 0.04812 | 0.9724 | 54 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | B1_R1 | 0.08439 | 0.0272 | 0.1508 | 0.2287 | 54 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | B1_R2 | 0.06608 | 0.008033 | 0.1321 | 0.3204 | 54 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | B1_R3_R5 | 0.05303 | 0.004255 | 0.112 | 0.3268 | 54 |
| HALLMARK_MTORC1_SIGNALING__NO_MP10 | B1_R6_R10 | 0.03852 | -0.0009287 | 0.08918 | 0.3243 | 54 |

这些是中心MP10与对应环平均程序的联系。R0不混入邻域平均；缺环不填0。旧重叠邻域仍采用原2000 μm空间块处理，未硬套新点级500 μm方案。P95/匹配对照、d_MP10曲线与前沿模型也在context85完整表中；二次项不解读成全程单调斜率。

## 2. P95 Epi自身程序随候选距离如何变化？

AddModuleScore主结果共有 139,655 个可计算的P95全范围相关(含敏感性程序)。在主61程序且n≥30、贡献空间块≥8的组合中，r<−0.1有 24,632 个，r>0.1有 26,264 个，|r|≤0.1有 30,844 个。这是重叠候选/程序扫描的描述性计数，不能当作独立生物学发现数。

下列覆盖充分的例子分别保留负向、正向和接近零，按效应和覆盖展示，并非只列显著结果。距离为候选相对丰度P95核心，而非到MP10-P95集合。

| example_type | section_id | candidate_type | program | pearson_r | n_spots | distance_min_um | distance_max_um | spatial_ci_low | spatial_ci_high | p_pearson | q_pearson_across_sections |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| negative | NPJ_sample4 | CD8 GZMK/Effector-memory-like | BUDDEN_TNF_EMT_RESPONSE_15 | -0.7512 | 79 | 0 | 1825 | -0.8743 | -0.4857 | 1.535e-15 | 5.739e-13 |
| negative | GSE319536_BC14 | GC B | HALLMARK_E2F_TARGETS | -0.7474 | 41 | 0 | 1114 | -0.8759 | -0.4618 | 1.983e-08 | 1.095e-06 |
| negative | CanDis_2026_Sample2 | Matrix_CAF | BUDDEN_TNF_EMT_RESPONSE_15 | -0.6882 | 205 | 1908 | 9629 | -0.7702 | -0.5644 | 4.147e-30 | 2.078e-25 |
| positive | NPJ_sample1 | CD8 GZMK/Effector-memory-like | BUDDEN_CHOLESTEROL_SYNTHESIS_27 | 0.8064 | 30 | 458.3 | 1682 | 0.3391 | 0.8905 | 7.459e-08 | 3.513e-06 |
| positive | CanDis_2026_Sample2 | GC B | BUDDEN_TNF_EMT_RESPONSE_15 | 0.733 | 205 | 0 | 2343 | 0.6047 | 0.8178 | 8.093e-36 | 8.111e-31 |
| positive | GSE319536_BC14 | Neutrophil | HALLMARK_E2F_TARGETS | 0.699 | 41 | 300 | 2835 | 0.3545 | 0.8489 | 3.708e-07 | 1.389e-05 |
| near_zero | NPJ_sample3 | Matrix_CAF | HALLMARK_E2F_TARGETS | -1.482e-05 | 148 | 100 | 755 | -0.279 | 0.1567 | 0.9999 | 0.9999 |
| near_zero | GSE285715_T1 | Venous EC | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | -1.645e-05 | 177 | 0 | 1136 | -0.1437 | 0.1659 | 0.9998 | 0.9999 |
| near_zero | CanDis_2026_Sample3 | SPP1 TAM | FLORES_CRISPR_PRO_INVASION_10 | 2.192e-05 | 191 | 0 | 818.5 | -0.1333 | 0.136 | 0.9998 | 0.9999 |

负r表示总体负线性关系，不证明每个距离段都下降；近零也可能含弯曲或平台。报告中的散点为原始spot自身AddModuleScore，负分保留，完整范围不因结果裁切。普通P/q基于独立spot假设；空间区间单独列，不能以q替代空间相关处理。

## 3. P95与非P95是否不同？

共同距离范围内、两组n≥30且各≥8空间块的组合中，18,124 个delta_r的联合块bootstrap区间不跨0；其中空间差异P的跨切片同族BH q<0.05为 5,565 个(主61程序 3,834 个)。计数含29候选及相关程序的重叠扫描，区间与Fisher-z空间检验不是完全相同的检验，因此支持数不必相同。

| section_id | candidate_type | program | program_family | n_P95 | n_nonP95 | common_min_um | common_max_um | r_P95 | r_nonP95 | delta_r | spatial_ci_low | spatial_ci_high | p_spatial_fisher_z | q_spatial_across_sections |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NPJ_sample1 | Treg | HALLMARK_MTORC1_SIGNALING__NO_MP10 | sensitivity | 30 | 524 | 100 | 1044 | -0.4778 | 0.4251 | -0.9029 | -1.232 | -0.06567 | 0.001027 | 0.02877 |
| GSE319536_BC02 | Lymphatic EC | HALLMARK_OXIDATIVE_PHOSPHORYLATION | main | 65 | 1058 | 100 | 1992 | -0.3311 | 0.4544 | -0.7855 | -1.004 | -0.4395 | 5.504e-08 | 7.273e-05 |
| CanDis_2026_Sample6 | C1QC TAM | PURAM_PEMT_TOP30 | main | 41 | 1017 | 0 | 1852 | 0.1478 | -0.5729 | 0.7208 | 0.4576 | 0.9279 | 1.094e-09 | 6.093e-06 |
| GSE285715_T3 | NK | BUDDEN_CHOLESTEROL_SYNTHESIS_27 | main | 97 | 1735 | 0 | 854.4 | -0.5423 | 0.06555 | -0.6078 | -0.7555 | -0.2556 | 1.525e-05 | 0.002443 |
| GSE171351_B1 | Neutrophil | FLORES_CRISPR_PRO_INVASION_10 | main | 33 | 527 | 100 | 556.8 | -0.5913 | -0.02037 | -0.571 | -0.7301 | -0.3485 | 1.144e-05 | 0.002048 |
| GSE238145_B15 | Mast | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | main | 73 | 1371 | 100 | 2052 | -0.3731 | 0.156 | -0.5291 | -0.6975 | -0.3325 | 1.207e-07 | 0.0001192 |
| Inhousedata_T4469479 | cDC1 | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | main | 80 | 1502 | 0 | 1418 | 0.2728 | -0.2418 | 0.5146 | 0.2243 | 0.7134 | 3.296e-05 | 0.003991 |

同一空间块抽样同时用于两组，重复抽中块保留重复权重。delta_r是相关强度差，不是单位距离下降速度；没有以“一组显著、另一组不显著”替代比较。full_range delta仅作描述，正式结果使用common_support，跨样本另提供common_support∩[0,1000] μm。两组实际距离方差、断档、零距离比例及n均保留。

## 4. 同dataset内的切片/样本异质性

| dataset_id | candidate_type | program | n_sections | n_units | r_min | r_max | n_negative | n_positive | n_near_zero | p_heterogeneity | q_heterogeneity |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NPJ__UNVERIFIED | SMC | BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | 3 | 3 | -1 | 1 | 1 | 2 | 0 | NA | NA |
| Inhousedata__UNVERIFIED | SMC | HALLMARK_GLYCOLYSIS | 7 | 7 | -1 | 0.9956 | 3 | 4 | 0 | NA | NA |
| GSE238145 | CD4 Naive/Memory-like | HALLMARK_G2M_CHECKPOINT | 9 | 8 | -0.8923 | 1 | 3 | 2 | 4 | 0.1797 | 0.3099 |
| CanDis__UNVERIFIED | C1QC TAM | BUDDEN_TG_CATABOLISM_24 | 13 | 13 | -0.9044 | 0.9255 | 5 | 3 | 5 | 0.0708 | 0.1249 |
| GSE171351 | Venous EC | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 3 | 3 | -0.7328 | 1 | 2 | 1 | 0 | NA | NA |
| GSE319536 | CD4 Helper/Conventional | HALLMARK_GLYCOLYSIS | 14 | 14 | -0.7461 | 0.9841 | 2 | 6 | 6 | 0.8473 | 0.8709 |

上表使用统一0–1000 μm内各片两组共同支持范围。异质性Q使用空间bootstrap Fisher-z标准误，未使用n_spot−3精度权重；描述性单位均值仍按切片等权Fisher-z合并。身份缺失者保留切片描述，不补造patient。表中出现双向效应说明值得逐片核查，只有直接区间/异质性检验有支持时才讨论统计差异。

| dataset_id | section_id | candidate_type | program | pearson_r | n_spots | distance_min_um | distance_max_um | spatial_ci_low | spatial_ci_high | q_pearson_across_sections |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NPJ__UNVERIFIED | NPJ_sample4 | SMC | BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | -1 | 2 | 888.8 | 984.9 | NA | NA | 1 |
| NPJ__UNVERIFIED | NPJ_sample6 | SMC | BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | 1 | 2 | 953.9 | 984.9 | NA | NA | 1 |
| Inhousedata__UNVERIFIED | Inhousedata_T4398343 | SMC | HALLMARK_GLYCOLYSIS | -1 | 2 | 916.5 | 984.9 | NA | NA | 1 |
| Inhousedata__UNVERIFIED | Inhousedata_T4471438 | SMC | HALLMARK_GLYCOLYSIS | 0.9956 | 4 | 556.8 | 1000 | NA | NA | 0.03964 |
| GSE238145 | GSE238145_B5 | CD4 Naive/Memory-like | HALLMARK_G2M_CHECKPOINT | -0.8923 | 6 | 200 | 655.7 | NA | NA | 0.1015 |
| GSE238145 | GSE238145_B18 | CD4 Naive/Memory-like | HALLMARK_G2M_CHECKPOINT | 1 | 2 | 100 | 173.2 | NA | NA | 1 |
| CanDis__UNVERIFIED | CanDis_2026_Sample11 | C1QC TAM | BUDDEN_TG_CATABOLISM_24 | -0.9044 | 5 | 100 | 346.4 | NA | NA | 0.164 |
| CanDis__UNVERIFIED | CanDis_2026_Sample2 | C1QC TAM | BUDDEN_TG_CATABOLISM_24 | 0.9255 | 3 | 608.3 | 916.5 | NA | NA | 0.5247 |
| GSE171351 | GSE171351_C1 | Venous EC | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | -0.7328 | 23 | 0 | 800 | -0.8868 | -0.491 | 0.001768 |
| GSE171351 | GSE171351_A1 | Venous EC | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 1 | 2 | 173.2 | 200 | NA | NA | 1 |
| GSE319536 | GSE319536_BC05 | CD4 Helper/Conventional | HALLMARK_GLYCOLYSIS | -0.7461 | 5 | 264.6 | 781 | NA | NA | 0.3935 |
| GSE319536 | GSE319536_BC12 | CD4 Helper/Conventional | HALLMARK_GLYCOLYSIS | 0.9841 | 3 | 692.8 | 866 | NA | NA | 0.3396 |

样本单位与已核实patient分别汇总，B5/B18按已有P1关系合并。重复切片有限，以上优先称“切片/样本间异质性”，不能仅凭一张切片命名稳定患者表型。GSE285715缺analysis_unit映射的切片保留全部效应及图册，不混进患者推断。

### 可用于机制对照的具体例子：GSE319536 / vCAF / Budden TNF–EMT response

| section_id | range_type | pearson_r | n_spots | distance_min_um | distance_max_um | spatial_ci_low | spatial_ci_high | p_pearson | q_pearson_across_sections |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE319536_BC11 | full_range | -0.04594 | 222 | 0 | 1418 | -0.2059 | 0.1334 | 0.4959 | 0.7101 |
| GSE319536_BC11 | common_support_0_1000 | -0.04269 | 214 | 0 | 953.9 | -0.2208 | 0.1544 | 0.5345 | 0.7709 |
| GSE319536_BC12 | full_range | -0.5949 | 39 | 0 | 916.5 | -0.7576 | -0.395 | 6.501e-05 | 0.001056 |
| GSE319536_BC12 | common_support_0_1000 | -0.5949 | 39 | 0 | 916.5 | -0.7576 | -0.395 | 6.501e-05 | 0.001669 |
| GSE319536_BC13 | full_range | -0.4693 | 45 | 100 | 900 | -0.7128 | -0.1286 | 0.001144 | 0.01046 |
| GSE319536_BC13 | common_support_0_1000 | -0.4693 | 45 | 100 | 900 | -0.7128 | -0.1286 | 0.001144 | 0.01482 |
| GSE319536_BC18 | full_range | 0.1871 | 111 | 0 | 2052 | -0.061 | 0.4627 | 0.04925 | 0.17 |
| GSE319536_BC18 | common_support_0_1000 | 0.2326 | 90 | 0 | 984.9 | 0.04047 | 0.4218 | 0.02736 | 0.1401 |

BC13完整范围r=−0.469 (n=45，100–900 μm)，非P95共同范围r=+0.286 (n=776)，delta_r=−0.755，联合空间区间[−1.043,−0.380]，空间P=3.5×10⁻⁵、跨片同族q=0.00416。BC12也负相关，但P95与非P95差异区间跨0，不能称为P95特有。BC18完整范围为弱正向且空间区间跨0；限定预设共同0–1000 μm窗口时r=+0.233 (n=90)，区间[0.040,0.422]，提示相反方向但对距离范围敏感。BC11接近零。这四片共同说明局部关系存在异质性，不能泛化为所有MP10高Epi都呈同样梯度。

同dataset统一窗口的vCAF–TNF程序P95相关异质性Q=42.116，空间标准误检验P=6.3×10⁻⁵、族q=0.000627 (15单位，14个有可用空间区间)。BC12与BC18的直接Fisher-z差异为−0.922，区间[−1.275,−0.570]，固定核心面板配对校正q=0.000289。这个对照优先于仅比较两片各自是否显著。BC13候选核心posterior mean丰度中位数1.76，但仅6个P95点有精确前沿距离、其中3个在前沿带，因此不能判断其关系是否独立于前沿位置。

## 5. 前沿和niche能否解释？

| analysis | program | effect | low | high | q_context85 | n_units | n_datasets |
| --- | --- | --- | --- | --- | --- | --- | --- |
| A_front_vs_internal | BUDDEN_CHOLESTEROL_SYNTHESIS_27 | -0.01589 | -0.1237 | 0.1037 | 0.9503 | 51 | 6 |
| A_front_vs_internal | BUDDEN_S1P_RESPONSE_18 | -0.00275 | -0.112 | 0.1369 | 0.9658 | 51 | 6 |
| A_front_vs_internal | BUDDEN_TG_CATABOLISM_24 | -0.1088 | -0.2029 | -0.02036 | 0.1421 | 51 | 6 |
| A_front_vs_internal | BUDDEN_TNF_EMT_RESPONSE_15 | 0.05461 | -0.07914 | 0.2215 | 0.6804 | 51 | 6 |
| A_front_vs_internal | FLORES_CRISPR_PRO_INVASION_10 | -0.04421 | -0.1127 | 0.02616 | 0.5571 | 51 | 6 |
| A_front_vs_internal | PURAM_PEMT_TOP30 | -0.01837 | -0.144 | 0.1082 | 0.9282 | 51 | 6 |
| A_front_vs_internal | BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | -0.03935 | -0.1442 | 0.05589 | 0.789 | 51 | 6 |
| C3_P95_state | BUDDEN_CHOLESTEROL_SYNTHESIS_27 | 0.04347 | -0.08209 | 0.1803 | 0.8709 | 23 | 5 |
| C3_P95_state | BUDDEN_S1P_RESPONSE_18 | -0.07728 | -0.2268 | 0.05761 | 0.7713 | 23 | 5 |
| C3_P95_state | BUDDEN_TG_CATABOLISM_24 | -0.1841 | -0.3851 | -0.01726 | 0.7713 | 23 | 5 |
| C3_P95_state | BUDDEN_TNF_EMT_RESPONSE_15 | 0.07251 | -0.1257 | 0.2453 | 0.8709 | 23 | 5 |
| C3_P95_state | FLORES_CRISPR_PRO_INVASION_10 | -0.1056 | -0.2101 | 0.01523 | 0.7713 | 23 | 5 |
| C3_P95_state | PURAM_PEMT_TOP30 | -0.0573 | -0.2643 | 0.1433 | 0.8709 | 23 | 5 |
| C1_all_Epi_continuous | EpiMP10_Top20_AUCell | -0.01331 | -0.08434 | 0.06496 | 0.8841 | 53 | 6 |
| C2_P95_enrichment | P95_probability | 0.008717 | -0.00666 | 0.02259 | 0.665 | 51 | 6 |
| C_active_front_MP10 | EpiMP10_Top20_AUCell | -0.08522 | -0.231 | 0.08011 | 0.665 | 50 | 6 |
| C3_P95_state | BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | 0.06453 | -0.06386 | 0.1773 | 0.9239 | 23 | 5 |

跨单位/数据集没有支持所有前沿都具有更高pEMT或Flores分数：前沿对内部的调整β分别−0.0184 (q=0.928)和−0.0442 (q=0.557)。MP10连续前沿距离、P95前沿富集及分子活跃前沿的MP10总体比较亦无一致支持(q分别0.884、0.665、0.665)。这不排除局部边段或单片现象，但不能把全队列边界统一命名为已验证侵袭前沿。


原有精确v3.4前沿、pEMT-high定义与R12标签保持冻结。缺精确前沿距离只限制位置分析，不限制候选距离相关。补充partial r先控制Epi比例与深度，再在同一精确前沿子集比较加入d_front前后；不能把换了spot集合导致的变化叫作前沿调整效应。


仅前沿带内仍存在片别的niche差异，例如CanDis_2026_Sample4的R12标签2在P95中低28.47个百分点 (标注P95/非P95 n=22/288，区间−42.55至−14.81个百分点，空间跨片q=0.0185)；GSE238145_B14的标签7低22.04个百分点 (n=12/313，区间−35.49至−9.84，q=0.0378)。这是同一前沿带内的比例差，不代表已匹配带内每一点的几何位置，也不是所有切片的共同niche偏好。

本轮额外完成严格[-200,0] μm带内的原P95 vs非P95九类niche比较，分母为各组已有正式标签的Epi，不在前沿内重新取Top5%。以下展示实际比例差，缺少空间区间的行仅为描述：

| section_id | niche | n_P95 | n_nonP95 | fraction_P95 | fraction_nonP95 | difference_percentage_points | ci_low_percentage_points | ci_high_percentage_points | q_spatial_across_sections |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CanDis_2026_Sample4 | 2 | 22 | 288 | 0 | 0.2847 | -28.47 | -42.55 | -14.81 | 0.01847 |
| GSE238145_B17 | 4 | 13 | 152 | 0.3846 | 0.1447 | 23.99 | -0.5053 | 51.05 | 0.3042 |
| Inhousedata_T4413900 | 0 | 15 | 183 | 0.2 | 0.08743 | 11.26 | -7.953 | 29.78 | 0.531 |
| GSE285715_T2 | 1 | 35 | 563 | 0.9143 | 0.8028 | 11.14 | -2.255 | 22.02 | 0.3042 |
| GSE319536_BC17 | 0 | 24 | 269 | 0.875 | 0.829 | 4.6 | -6.412 | 14.84 | 0.6534 |

新增文献程序的非配对单细胞来源支持(样本层UCell中位数的中位数)：

| program | CAF | Epithelial | Myeloid | SMC |
| --- | --- | --- | --- | --- |
| BUDDEN_CHOLESTEROL_SYNTHESIS_27 | 0.01145 | 0.03115 | 0.01845 | 0.01641 |
| BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18 | 0.0007152 | 0.01489 | 0.006295 | 0 |
| BUDDEN_S1P_RESPONSE_18 | 0.2534 | 0.1128 | 0.1243 | 0.2965 |
| BUDDEN_TG_CATABOLISM_24 | 0.03978 | 0.05869 | 0.04058 | 0.05858 |
| BUDDEN_TNF_EMT_RESPONSE_15 | 0.3233 | 0.2206 | 0.2688 | 0.4007 |

S1P响应和TNF–EMT响应集合可在CAF/SMC等来源显著表达，空间相关不能直接称为Epi专属信号。胆固醇去重叠版本单列用于区别共享基因一致性；这些UCell来源分数与空间SCT AddModuleScore不能直接比较绝对数值。

上文共同范围示例中的组成/深度及前沿敏感性(前沿加入前后严格同点)：

| section_id | candidate_type | program | pearson_r | partial_r_epi_depth | n_partial | partial_r_front_subset_before | partial_r_front_subset_after | n_front_subset |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NPJ_sample1 | Treg | HALLMARK_MTORC1_SIGNALING__NO_MP10 | -0.4778 | -0.2895 | 30 | NA | NA | 2 |
| CanDis_2026_Sample6 | C1QC TAM | PURAM_PEMT_TOP30 | 0.1478 | 0.07963 | 41 | -0.2879 | -0.2433 | 17 |
| GSE171351_B1 | Neutrophil | FLORES_CRISPR_PRO_INVASION_10 | -0.5913 | -0.5441 | 33 | -0.543 | -0.5759 | 16 |
| GSE238145_B15 | Mast | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | -0.3731 | -0.2934 | 73 | -0.2517 | -0.2612 | 62 |
| GSE285715_T3 | NK | BUDDEN_CHOLESTEROL_SYNTHESIS_27 | -0.5423 | -0.5278 | 97 | -0.04177 | -0.1126 | 22 |
| GSE319536_BC02 | Lymphatic EC | HALLMARK_OXIDATIVE_PHOSPHORYLATION | -0.3311 | -0.4477 | 65 | -0.4966 | -0.1586 | 19 |
| Inhousedata_T4469479 | cDC1 | HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY | 0.2728 | 0.04277 | 80 | -0.3158 | -0.3159 | 32 |

协变量完整点数小于10的partial r单独标记为稀疏描述；特别是5–6个点在加入多个协变量后可因剩余自由度极少接近±1，不能用来证明调整稳健性。原始Pearson仍按可计算点数保留，不因无法稳定调整而删除。


全部九类及未标注覆盖保留在niche_counts_all；原前沿位置匹配和一般前沿分布保留原名称，未改名为仅前沿带比较。R12解释使用原R1/R2 profile，不用中心spot组成冒充邻域组成。

## 6. 优先审阅与证据限制

| dataset_id | n_sections | n_P95_ge30 | n_units | n_verified_patients | n_exact_front_sections | n_core_pairs_with_opposite_section_directions |
| --- | --- | --- | --- | --- | --- | --- |
| GSE319536 | 16 | 14 | 16 | 16 | 15 | 393 |
| CanDis__UNVERIFIED | 13 | 11 | 13 | 0 | 13 | 386 |
| NPJ__UNVERIFIED | 7 | 7 | 7 | 0 | 7 | 328 |
| GSE238145 | 10 | 6 | 9 | 9 | 9 | 349 |
| Inhousedata__UNVERIFIED | 7 | 5 | 7 | 0 | 7 | 336 |
| GSE285715 | 6 | 3 | 0 | 0 | 3 | 135 |
| GSE171351 | 4 | 1 | 4 | 4 | 3 | 197 |

优先顺序按P95数量、共同距离覆盖、同dataset重复单位及可用于位置解释的前沿覆盖安排；具体候选–程序和切片见section_relationship_priorities。方向相反的切片是机制对照线索，未因“难看”删除。罕见类型的P95只是相对排名，候选核心绝对丰度和实际目标比例见candidate_targets，不能仅凭相对P95称大量存在。

主要覆盖限制：1 张无冻结Epi；2 张仅1个内部P95 Epi；47 张达到P95≥30。GSE238145_B13保留为无Epi模块状态；GSE171351_D1仅1个内部Epi；GSE319536_BC15有4个内部Epi、仅1个P95，两片均不能计算P95相关。三张邻近正常只作合格Epi参照，没有伪造肿瘤P95或前沿。

所有85项AddModuleScore以完整可用SCT data基因池评分，SCT拟合要求至少5个细胞检出，故少数低检测成员可能不进入该切片SCT池；各片实际成员/缺失和RNA检测数完整记录。新5项UCell按完整RNA排名宇宙补齐。单细胞来源支持复用上一轮19,053个分层抽样细胞，按真实样本和原注释汇总；与空间不配对者仅支持可能来源。混合spot信号不能直接写成Epi特异表达、脂质转移、S1P浓度或代谢通量。空间500 μm块近似、共同范围缺口和少量高端P95仍限制精度。

63片评分、四项候选比较、旧邻域/前沿/niche整理、新前沿带内niche比较、单细胞来源及逐片/数据集图册均已实际运行。没有重跑cell2location、niche、CNV或通讯。旧结果保持不变。
