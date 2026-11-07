# MP10 数据集组合、距离关系和样本异质性：本轮实际结果

本轮已完成。新增的是**小表层面的全组合搜索、真实r/partial r的统一单位汇总、跨输入异质性比较以及niche生态轴对距离结果的解释模型**。没有重跑MP10/UCell/AddModuleScore、cell2location、点级距离或niche聚类，没有恢复点级parquet。既有空间niche标签仍为R12 0–8；本轮Ward分组是分析单元层的候选分组，不是新的spot niche。

登记63片(60肿瘤来源、3正常参照)。全Epi中心端点实际可用58张肿瘤切片，对应54个已有analysis units加3个缺单位映射的切片记录，共57条合并记录；不能把这57条全部称为已核实患者或独立分析单元。单位推断主参照因此为6个dataset；七dataset描述保留GSE285715的三张T切片。不同端点有不同覆盖，完整记录见input_coverage和section_manifest。

主名称61个，按完整成员合并Budden/Reactome胆固醇别名后为60个canonical主程序：22个有现成noMP10版本，38个原集合本就无共享成员。相同成员的别名不在同一矩阵中重复加权；以Budden实际评分为统一代表，没有把两个AddModuleScore数值强行当成完全相同。Top100仅替换Top30作敏感性。

共登记19,419个统计端点配置，保存3,582,729条数据集组合汇总。七dataset所有127个非空组合均在支持它们的端点中遍历；仅6个dataset有单位身份的矩阵有63个非空组合。第一层异质性覆盖87类输入配置、24,090条组合/K或状态记录。包括单dataset、无需Inhouse的组合；没有在dataset内部按效应挑片。

## 1. 新增和直接复用的部分

A使用旧80与新增5各自真实保存的R0/R1/R2 r、partial r，按origin去重；B使用B3_Epi_distance/B3_tissue_distance的原逐spot r/partial r，0–1000 μm和零距离保持；C1/C2复用已算出的AMS/UCell候选距离相关、共同支持及联合空间delta；C3/C4复用距离摘要和全Epi MP10相关。新结果的条件P/q来自真实单位层，不借用旧二次模型、β或原始Pearson的P/CI。原q列仍在input_index所指源文件中，未改名冒充本轮校正。

每个统计量先按冻结单位合并重复切片，r在Fisher-z上等权；配对delta直接平均各切片的配对delta，不由两组汇总r相减。对于缺单位映射的记录，只保留section-only描述；正常来源不进入肿瘤异质性建轴。未找到以dataset_reanalysis/mp10_dataset_reanalysis/mp10_inhouse_pca_sweep命名的当前本地缓存，已由原始效应小表重建，不以目录名称缺失阻塞。

## 2. 两类距离轴及全数据参照

下表均为全部六个有身份dataset的等权参照。B3是到**MP10-P95 Epi参考集合**的距离；C1是P95 Epi自身程序到**候选细胞相对富集核心**的距离。二者不是同一个问题。conditional区间/P以所列dataset固定、单位为重复；between检验以dataset均值为重复。所有q_across_search覆盖同端点族的全部组合扫描，不能将条件均值支持等同于可推广的跨队列证据。

| module | canonical_program | candidate_type | estimate | conditional_low | conditional_high | n_analysis_units | unit_direction_fraction | q_conditional_across_search | q_between_across_search |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B1_R0 | HALLMARK_ANDROGEN_RESPONSE | none | 0.2245 | 0.1713 | 0.2764 | 54 | 0.8889 | 8.221e-07 | 0.04632 |
| B3_Epi_distance | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | none | -0.2044 | -0.2487 | -0.1592 | 53 | 0.9057 | 1.141e-07 | 0.04563 |
| B3_Epi_distance | HALLMARK_MTORC1_SIGNALING | none | -0.2758 | -0.3337 | -0.2159 | 53 | 0.9623 | 3.713e-06 | 0.04563 |
| B3_tissue_distance | HALLMARK_MYC_TARGETS_V1 | none | -0.3381 | -0.3858 | -0.2886 | 54 | 0.9444 | 5.793e-10 | 0.005177 |
| C1 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | Classical Mono | -0.3082 | -0.375 | -0.2382 | 47 | 0.9149 | 6.785e-05 | 0.17 |
| C1 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | Neutrophil | 0.309 | 0.1751 | 0.4316 | 45 | 0.9111 | 0.03386 | 0.1714 |
| C2 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | Plasma | 0.09287 | 0.03884 | 0.1469 | 46 | 0.6304 | 0.3968 | 0.9462 |
| C3 | distance | Treg | 63.49 | 10.15 | 116.8 | 54 | 0.5741 | 0.3515 | 0.772 |

B3非脂代谢方向值得关注：Epi中的UPR随离MP10高分参考点更远而总体降低(r=−0.2044，53单位、6dataset，48/53方向为负)；mTORC1 r=−0.2758，51/53为负。全组织MYC_TARGETS_V1 r=−0.3381，54单位，51/54为负。它们是混合spot程序的总体线性关联，不能据此推断所有距离段单调或细胞内代谢通量。

C1中Classical Mono–TNFA/NFKB为负(r=−0.3082，47单位，43/47负向)，Neutrophil–OXPHOS为正(r=+0.3090，45单位，41/45正向)。前者是MP10-P95 Epi内部的炎症程序梯度，后者意味着远离该候选核心时OXPHOS程序总体较高；不能把常见细胞组成关系自动称为MP10特异机制。它们的条件q分别约6.8×10⁻⁵和0.0339，但dataset-replicate跨扫描q约0.170/0.171，跨队列外推证据弱于条件单位层证据。

## 3. P95内部梯度与P95/nonP95差异必须分开

相同Classical Mono–TNFA关系的六dataset配对delta约−0.0107，条件区间[−0.0796,0.0582]，跨搜索q=0.979；Neutrophil–OXPHOS的delta约+0.0996，区间[−0.0155,0.2148]，q=0.833。故其P95内部梯度不能升级为已证实的P95特有关系。

C2存在选定dataset组合的探索性差异：CanDis、GSE171351、GSE238145(D007)的Plasma–OXPHOS配对delta=+0.1993，21单位，条件区间[0.1629,0.2356]，跨该扫描族q=0.0009256。但以3个dataset均值为重复的跨搜索q=0.9462，全部六dataset也没有C2主面板q<0.05(最低约0.1149)。这是**指定组合中的条件差异候选**，不是独立验证或全队列P95差异。所有单片旧联合空间区间及其字段在源表保留；本轮条件区间是单位层变异，不冒充原空间区间。

C3并未给出普遍“P95更靠近Treg”的证据。全六dataset均值差为+63.49 μm(P95更远)，但跨搜索q≈0.352，不能宣称排斥。均值差、中位数、零比例差和匹配摘要分别保留；匹配源结果只有描述性时，没有借用未匹配距离的空间CI/P。

## 4. 哪些完整dataset组合更支持？

| module | canonical_program | candidate_type | combination | datasets | estimate | n_analysis_units | q_conditional_across_search | q_between_across_search |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B1_R0 | HALLMARK_ANDROGEN_RESPONSE | none | D117 | CanDis__UNVERIFIED / GSE238145 / GSE319536 / Inhousedata__UNVERIFIED / NPJ__UNVERIFIED | 0.2406 | 51 | 2.681e-08 | 0.05818 |
| B3_Epi_distance | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | none | D115 | CanDis__UNVERIFIED / GSE171351 / GSE319536 / Inhousedata__UNVERIFIED / NPJ__UNVERIFIED | -0.2323 | 45 | 1.15e-07 | 0.04456 |
| B3_Epi_distance | HALLMARK_MTORC1_SIGNALING | none | D113 | CanDis__UNVERIFIED / GSE319536 / Inhousedata__UNVERIFIED / NPJ__UNVERIFIED | -0.3233 | 42 | 1.547e-08 | 0.04757 |
| B3_tissue_distance | HALLMARK_MYC_TARGETS_V1 | none | D119 | CanDis__UNVERIFIED / GSE171351 / GSE238145 / GSE319536 / Inhousedata__UNVERIFIED / NPJ__UNVERIFIED | -0.3381 | 54 | 5.793e-10 | 0.005177 |
| C1 | HALLMARK_TNFA_SIGNALING_VIA_NFKB | Classical Mono | D085 | CanDis__UNVERIFIED / GSE238145 / GSE319536 / NPJ__UNVERIFIED | -0.3578 | 38 | 0.000295 | 0.1976 |
| C1 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | Neutrophil | D083 | CanDis__UNVERIFIED / GSE171351 / GSE319536 / NPJ__UNVERIFIED | 0.3853 | 34 | 0.07005 | 0.1599 |
| C2 | HALLMARK_OXIDATIVE_PHOSPHORYLATION | Plasma | D007 | CanDis__UNVERIFIED / GSE171351 / GSE238145 | 0.1993 | 21 | 0.0009256 | 0.9462 |
| C3 | distance | Treg | D081 | CanDis__UNVERIFIED / GSE319536 / NPJ__UNVERIFIED | 104.5 | 36 | 0.3491 | 0.772 |

上表按效应、单位/数据集方向一致性及覆盖综合排序，至少两个dataset和15单位用于优先讨论；这只是展示规则，完整搜索仍包括所有单dataset和小组合。Classical Mono–TNFA的较强组合D085为CanDis、GSE238145、GSE319536、NPJ，**不包含Inhouse**；其r=−0.3578，但dataset层q=0.198。UPR较强组合D115只排除GSE238145；mTORC1较强组合D113排除GSE171351/GSE238145。这些排除并未改变全数据参照，也不能被解释成某dataset“错误”。七dataset含section-only记录的描述与六dataset单位推断分列，不以缺失填0。

## 5. 哪种输入刻画异质性更清楚？

| input_family | n_units | n_features | PC1_variance | PC2_variance | PC1_dataset_eta2 | PC1_log_coverage_r | K | silhouette | cluster_sizes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B1_R0_partial_r_full | 54 | 60 | 0.245 | 0.184 | 0.02334 | -0.3343 | 2 | 0.5673 | 53 / 1 |
| B3_Epi_distance_partial_r_full | 53 | 60 | 0.3237 | 0.09781 | 0.03556 | 0.2473 | 2 | 0.1622 | 36 / 17 |
| C1_CAF_full | 49 | 26 | 0.2229 | 0.1259 | 0.123 | 0.1588 | 2 | 0.1194 | 32 / 17 |
| niche_all_Epi_all | 54 | 9 | 0.3985 | 0.2417 | 0.07603 | 0.2041 | 4 | 0.3492 | 20 / 14 / 4 / 16 |
| cell_composition_whole_tissue | 56 | 9 | 0.5383 | 0.2065 | 0.3436 | 0.6367 | 2 | 0.3167 | 42 / 14 |

在深入检查的五类代表输入中，**全Epi九类niche组成最稳定且最容易解释**。54单位，PC1/PC2解释39.85%/24.17%，PC1的dataset解释比例约7.6%。PC1主要是标签0正向(载荷+0.736)与标签1负向(−0.656)的生态构成对比。原R1参考中，标签0的Epithelial比例约82.1%；标签1约59.1%，同时B比例约18.6%。这是原邻域profile背景，不是本轮中心spot纯度，更不能叫全组织细胞组成。

Ward的描述性最优K为4，组大小20/14/4/16，silhouette≈0.349；80%单位留出平均ARI≈0.833、特征扰动≈0.705、K-means敏感性≈0.765。这支持连续生态轴和若干内部稳定候选组，尚不支持正式生物学亚型。组内niche组成图使用建组变量本身，只是解释，不是独立验证。

R0 partial谱虽有较高silhouette≈0.567，却是53+1；单独单位是CanDis Sample12(n_Epi=405)，K-means一致性低，不能包装为稳定二分类。B3 Epi partial谱和CAF候选梯度更接近连续变化，候选聚类稳定性有限。全组织九大谱系的PC1同时与dataset(η²≈0.344)和Epi点数(r≈0.637)相关，解释时需警惕组成/覆盖驱动。其余核心及扩展60程序分块、P95 R0状态、匹配状态差异及internal/front/near/deep niche均完整进入第一层表，没有因结构弱而删除。

| input_family | check | count | mean | min | max |
| --- | --- | --- | --- | --- | --- |
| B1_R0_partial_r_full | feature_subset_80percent | 100 | 0.7155 | -0.02242 | 1 |
| B1_R0_partial_r_full | kmeans_same_K | 10 | 0.2221 | -0.0274 | 1 |
| B1_R0_partial_r_full | unit_subset_80percent | 100 | 0.77 | 0 | 1 |
| B3_Epi_distance_partial_r_full | feature_subset_80percent | 100 | 0.4946 | -0.0741 | 1 |
| B3_Epi_distance_partial_r_full | kmeans_same_K | 10 | 0.2816 | -0.0435 | 0.8501 |
| B3_Epi_distance_partial_r_full | unit_subset_80percent | 100 | 0.4647 | -0.1022 | 1 |
| C1_CAF_full | feature_subset_80percent | 100 | 0.4939 | -0.02432 | 0.9177 |
| C1_CAF_full | kmeans_same_K | 10 | 0.3829 | -0.05984 | 0.6255 |
| C1_CAF_full | unit_subset_80percent | 100 | 0.5517 | 0.08065 | 1 |
| cell_composition_whole_tissue | feature_subset_80percent | 100 | 0.3391 | 0.03509 | 0.9236 |
| cell_composition_whole_tissue | kmeans_same_K | 10 | 0.3045 | 0.2769 | 0.3164 |
| cell_composition_whole_tissue | unit_subset_80percent | 100 | 0.4215 | 0.006463 | 1 |
| niche_all_Epi_all | feature_subset_80percent | 100 | 0.7055 | 0.2897 | 0.9149 |
| niche_all_Epi_all | kmeans_same_K | 10 | 0.7647 | 0.6051 | 0.8804 |
| niche_all_Epi_all | unit_subset_80percent | 100 | 0.833 | 0.5816 | 1 |

## 6. niche生态背景能解释MP10距离关系吗？

已实际拟合29,664个生态轴–端点关系，生态PC1/PC2只用all_Epi niche比例构建，未使用距离r/delta筛轴。模型使用真实分析单位，控制dataset项，另检查技术覆盖、全组织Epithelial比例、可用空间SE权重和dataset留出。**没有关系通过全部生态搜索BH q<0.05**。

| module | candidate_type | canonical_program | statistic | n_units | beta | low | high | p_HC3 | q_within_endpoint_family | q_across_ecology_search | between_dataset_r | composition_adjusted_beta | technical_adjusted_beta |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C1 | vCAF | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | partial_r | 47 | -0.0959 | -0.1365 | -0.05534 | 2.482e-05 | 0.08636 | 0.5487 | -0.3038 | -0.08613 | -0.09501 |
| C2 | NK | HALLMARK_DNA_REPAIR | paired_delta_r | 48 | 0.1139 | 0.05856 | 0.1693 | 0.0001642 | 0.5715 | 0.9743 | 0.03777 | 0.1218 | 0.1121 |

两个可保留的预验证线索：vCAF–UPR的**partial r**随生态PC1增加更负(每1 SD PC1的Fisher-z系数−0.0959，47单位，区间[−0.1365,−0.0553])；NK–DNA_REPAIR配对delta随PC1增加更正(β=+0.1139，48单位，区间[0.0586,0.1693])。六次dataset留出分别保持负向/正向，但全搜索q为0.549/0.974，不能声称生态背景已经解释了这些差异。前者在换回原始r、UCell或其他距离范围时减弱；partial r没有可借用的原始r空间SE。后者可用空间SE加权在同点子集约+0.101，仍属于内部敏感性。

两例的between-dataset相关分别约−0.304和+0.038，不能将内部生态轴关系直接等同于dataset均值差异；逐dataset单位关联和留出结果已单独保存。已有两PC和dataset模型的增量R²是样本内拟合解释量，不是外部预测性能。

## 7. signature、调整、评分、距离和单dataset敏感性

| canonical_program | signature | statistic | estimate | conditional_low | conditional_high | q_conditional_across_search | q_between_across_search |
| --- | --- | --- | --- | --- | --- | --- | --- |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | full | r | -0.2044 | -0.2487 | -0.1592 | 1.141e-07 | 0.04563 |
| HALLMARK_UNFOLDED_PROTEIN_RESPONSE | full | partial_r | -0.04752 | -0.07833 | -0.01662 | 0.03226 | 0.2629 |
| HALLMARK_MTORC1_SIGNALING | full | r | -0.2758 | -0.3337 | -0.2159 | 3.713e-06 | 0.04563 |
| HALLMARK_MTORC1_SIGNALING | noMP10 | r | -0.2497 | -0.3071 | -0.1904 | 1.589e-05 | 0.03277 |
| HALLMARK_MTORC1_SIGNALING | full | partial_r | -0.146 | -0.1834 | -0.1082 | 2.378e-05 | 0.2231 |
| HALLMARK_MTORC1_SIGNALING | noMP10 | partial_r | -0.1101 | -0.1437 | -0.07636 | 0.0001003 | 0.1625 |

mTORC1去共享后原始r仍为−0.2497；再控制组成/深度后的partial r减至−0.1101，条件q=0.000100，但dataset层q=0.163。UPR无共享成员，调整后由−0.2044减至−0.0475，dataset层q=0.263。故部分梯度与组成/深度共同变化，不能仅凭noMP10结果宣称细胞内独立机制。

| input_family | variant | n_paired_units | n_common_features | K | adjusted_rand | absolute_PC1_correlation |
| --- | --- | --- | --- | --- | --- | --- |
| C1_CAF_full | UCell_RNA / common_support_0_1000 / r / n10 | 49 | 26 | 2 | -0.01551 | 0.8125 |
| C1_CAF_full | AddModuleScore_SCT_data / common_support / r / n10 | 49 | 26 | 2 | -0.01901 | 0.6685 |
| C1_CAF_full | AddModuleScore_SCT_data / full_range / r / n10 | 49 | 26 | 2 | -0.01901 | 0.6608 |
| C1_CAF_full | AddModuleScore_SCT_data / common_support_0_1000 / partial_r / n10 | 49 | 26 | 2 | -0.01188 | 0.4351 |
| C1_CAF_full | AddModuleScore_SCT_data / common_support_0_1000 / r / n30 | 40 | 13 | 2 | 0.06749 | 0.8503 |
| B1_R0_partial_r_full | Puram_Top100_replaces_Top30 | 54 | 60 | 2 | 1 | 0.9998 |
| B3_Epi_distance_partial_r_full | Puram_Top100_replaces_Top30 | 53 | 60 | 2 | 1 | 0.9999 |

CAF候选梯度的连续PC1在AMS/UCell之间仍有一定对应(|r|≈0.81)，但离散分组ARI接近0；范围和partial处理同样改变分组。n≥30检查在保留下来的单位上比较原完整拟合标签与新拟合，不能将两个完全相同的子矩阵ARI=1当稳定性。Top100替换单一Top30主题对整体轴影响很小，这只是内部冗余敏感性。

C4另有明确的稀疏影响：原全部Epi汇总中，GSE171351_D1仅3个查询点，近±1的r在Fisher-z上对dataset均值影响很大。没有删除原主结果；单列总Epi n≥20敏感性后，Venous EC r由约0.394降至0.116，Matrix_CAF由约0.368降至0.083；进一步partial后约0.043/0.040，跨该敏感性搜索q均>0.1。这是全Epi覆盖问题，不是P95门槛，也不应据原大r给样本命名特有表型。

## 8. 下一轮优先验证的五个方向

1. **MP10-P95周围的mTORC1/UPR梯度。** 多dataset方向支持，mTORC1去共享仍存在；需要新队列预设验证和更明确的上皮来源表达/组成控制，以区分混合spot背景。
2. **P95 Epi中的Classical Mono–TNFA/NFKB关系。** 目前支持P95内部梯度，不支持P95相对nonP95的特殊性；下一轮需同时测量两组、控制空间位置并验证候选核心定位。
3. **Neutrophil–OXPHOS的正向距离关系。** 可作为与炎症梯度不同方向的对照；需核查稀有细胞定位、评分方法和低氧/组织结构，不能推断代谢转移。
4. **Plasma–OXPHOS的指定dataset组合delta。** 条件单位层支持，但跨dataset推断和全数据不支持泛化；需在未参与搜索的队列中固定组合对应问题，不能重复挑阳性dataset。
5. **R12连续生态轴及候选四组。** 先在独立空间队列验证生态组成轴/分组，再预先限定检验vCAF–UPR、NK–DNA_REPAIR等少数端点；当前FDR未通过，应明确保持探索性。

## 9. 可选缺项和最小提取需求

| input | status | minimum_needed |
| --- | --- | --- |
| 14_Epi_MP_sample_activity | not_found_as_compatible_small_table | section_id, frozen Epi population, MP1..MP14 mean/median/q25/q75/n; no full expression required |
| whole_tissue_major_lineage | available_exact_weighted_disjoint_region_means | none; 9 major lineages, not 30-subtype composition |
| formal_R12_spatial_adjacency | no_compatible_whole_cohort_small_summary | section_id, R12 pair adjacency counts, eligible edges, tissue-component-aware denominator |
| molecular_subtype_clinical | no_verified_crosswalk_linked_summary_used | section/sample to verified subtype/pathology/clinical label and provenance |
| existing_dataset_combination_cache | named reanalysis/sweep directories not found; rebuilt from original effect summaries | none |

已纳入的是真实内侧Epi UCell中位数，以及原P95参考R0均值、P95匹配R0程序差；不是14个Epi MPs的活性替代品。全组织组成来自分区内9大谱系均值按有效spot数精确加权，不用candidate core丰度或niche比例冒充全组织组成。若后续需要全30类型组成，最小提取为section_id、有效组织n、各类型mean fraction/abundance及覆盖；若需要空间组织，需额外输出正式R12邻接计数及边分母。没有验证的临床/分子分型或scRNA–空间对应关系不用于本轮解释。

本轮没有声称重新评估所有样本的曲线形状，也没有生成或虚构逐spot散点。图中的点是已有分析单位效应，缺身份的section-only记录明确区别；曲线摘要不完整不妨碍已有真实r的组合遍历。

## 10. 文件入口和体积

- `tables/input_coverage.tsv`、`canonical_programs.tsv`、`section_manifest.tsv`: 输入、成员映射及全部切片状态。
- `tables/unit_endpoints.tsv.gz`、`unit_background.tsv.gz`: 唯一统一单位/背景输入；含无法形成单位的section-only记录。
- `tables/endpoints.tsv.gz`、`combinations.tsv`: 统计端点和完整dataset组合字典。
- `tables/combination_associations.tsv.gz`: 全部组合效应及两层P/q；字段明确条件推断与dataset-replicate推断。
- `tables/heterogeneity_configurations.tsv.gz`: 核心/扩展面板全组合PC1/PC2、带符号主载荷、K/规模及技术关系。
- `tables/ecology_associations.tsv.gz`、`ecology_leave_dataset_out.tsv`、`ecology_sensitivity.tsv.gz`: 完整生态模型及代表敏感性。
- `tables/discussion_relations.tsv`、`discussion_selected_datasets.tsv`、`discussion_sensitivity.tsv.gz`: 精简讨论入口。
- `tables/figure_index.tsv`、`figures/`: 九组、16张独立代表图，Arial白底、矢量PDF及300dpi PNG；有图例时三版PDF。
- `scripts/run_all.sh`、`config.json`、`METHODS.md`: 可重启运行入口和实际定义。

未创建ZIP或重复上传目录。全部源结果只读。本轮不恢复点级缓存、不复制表达/SCT/反卷积对象，只保留必要矩阵和汇总。新增文件实际合计约 **94.15 MB (89.78 MiB)**，低于100 MB；未删除任何旧项目结果。
