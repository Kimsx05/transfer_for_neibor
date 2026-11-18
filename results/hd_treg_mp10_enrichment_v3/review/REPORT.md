# hd_treg_mp10_enrichment_v3

完成日期：2026-11-18 15:30 CST。这是查看过v1/v2后的探索性区域空间分析。五张切片全部纳入；T32是表达聚合区域，不是Treg细胞。没有运行通讯或上游重分析。

## 主要结论

简化标签增加了可描述的FOXP3支持T区域，但没有形成跨切片一致的MP10-high偏好。H05主FE计划10项，8项可计算限制置换P，BH后均未检出明确富集(最小q=0.8916)。50μm有严重的测量机会不足；MPBC1的阳性可交换支持不足，MPBC2无H05曝露，不能将NA解释为没有生物邻近。

200μm连续Treg程序在D1呈小幅负向调整效应(β≈−0.036)，400/800μm区块区间均偏负；其他四片区间跨零。它与FOXP3标签的FE接近1同时成立，不能把两个读数当作同一身份指标或互相独立的验证。本轮没有新增连续检验P；v1正式连续P/q另表原样保留。

## 定义与测量

T_eligible核验与历史一致：nFeature≥325；CD3D/CD3E/CD3G/TRAC/CD2中至少2个检出且原始UMI合计≥3。连续11基因UCell使用全部T_eligible，包括零值；完整背景、average ties、maxRank1300完全复用。新标签仅T_eligible且raw FOXP3>0，没有任何附加门槛。旧strict仍为原3基因/core/旧cutoff联合规则，旧cutoff=0.19955389201641047。

| section | T_eligible | FOXP3_supported | strict | intersection | FOXP3_only | strict_only |
| --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 280 | 39 | 2 | 1 | 38 | 1 |
| D1_Johnson | 7103 | 1254 | 51 | 27 | 1227 | 24 |
| UC_SCC | 3579 | 285 | 13 | 3 | 282 | 10 |
| MPBC1 | 547 | 35 | 2 | 1 | 34 | 1 |
| MPBC2 | 2270 | 179 | 2 | 2 | 177 | 0 |

差异不是重新评分或更多真实Treg被证实。旧strict不要求FOXP3，因此两标签不是简单包含关系。FOXP3-only大量未超过旧UCell门槛，一部分不足3个panel基因；strict-only则可以凭其他core/程序信号满足旧规则。支持marker、检测量和程序分布见`label_expression_distributions.tsv`及skill小提琴图。

| section | n | zero_program | median_program | median_nFeature | less3_signature_genes | below_old_cutoff |
| --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 39 | 12 | 0.06112 | 2309 | 14 | 38 |
| D1_Johnson | 1254 | 594 | 0.007921 | 3456 | 384 | 1227 |
| UC_SCC | 285 | 72 | 0.04001 | 1903 | 165 | 282 |
| MPBC1 | 35 | 3 | 0.04682 | 1742 | 27 | 34 |
| MPBC2 | 179 | 39 | 0.05364 | 1655 | 114 | 177 |

上述FOXP3阳性零程序分别为12/39、594/1254、72/285、3/35、39/179。实际秩例子(Cervilla_T32_000033)FOXP3=1但平均秩1356.5，超出1300窗口，UCell=0。`historical_zero_UCell_rank_examples.tsv`保留原v2例子。FOXP3单个原始计数、程序窗口及宽度要求衡量不同事实。活化/Tfh可贡献共享程序是竞争解释，不能仅靠这些分布证明真实Treg身份；本轮没有增加新硬排除条件。

## 校准分母修复

复用v2的21600个冻结混合物及原训练/样本留出，未重生成或重划分。逐条件先应用完整候选规则，以**所有T_eligible阴性混合物**作为分母。最小可行非负`score > cutoff`为0.27019811718420683；此前v2有错误条件分母的候选0.319938没有被复用。本轮候选单独保存，既不替换旧strict，也不进入FOXP3定义。

| partition | rule | metric | numerator | T_eligible_denominator | rate |
| --- | --- | --- | --- | --- | --- |
| sample_holdout | denominator_fixed | FPR | 1 | 1249 | 0.0008006 |
| sample_holdout | denominator_fixed | detection | 160 | 3681 | 0.04347 |
| sample_holdout | historical | FPR | 6 | 1249 | 0.004804 |
| sample_holdout | historical | detection | 396 | 3681 | 0.1076 |
| train | denominator_fixed | FPR | 1 | 1050 | 0.0009524 |
| train | denominator_fixed | detection | 248 | 3187 | 0.07782 |
| train | historical | FPR | 13 | 1050 | 0.01238 |
| train | historical | detection | 522 | 3187 | 0.1638 |

表中池化率仅是摘要，阈值选择约束是每个非空训练条件。训练唯一FP在Epi90/Tfh/5000条件，为1/30；留出唯一FP在Epi40/Tfh/5000，为1/59。小条件区间很宽；无T_eligible阴性的条件率为NA。完整逐条件分子/分母、Wilson描述区间与最小性验证见calibration表。留出已在v2看过，不是新的独立验证；模拟重复也不是患者重复。修复后仍低检出，不能据此宣称区域身份问题已解决。

## MP10和有效空间总体

历史Top20、18085基因背景、平均并列秩、K=905、Epi64测量资格和分数均未变。H05/H10分别为每片合格Epi64的严格大于q95/q90(linear/type7)，没有拆分并列值。

| section | definition | threshold | n_eligible | n_high | actual_fraction |
| --- | --- | --- | --- | --- | --- |
| Cervilla | H05 | 0.1094 | 1707 | 86 | 0.05038 |
| Cervilla | H10 | 0.09263 | 1707 | 171 | 0.1002 |
| D1_Johnson | H05 | 0.1158 | 2884 | 145 | 0.05028 |
| D1_Johnson | H10 | 0.1053 | 2884 | 289 | 0.1002 |
| UC_SCC | H05 | 0.1046 | 4141 | 207 | 0.04999 |
| UC_SCC | H10 | 0.08924 | 4141 | 414 | 0.09998 |
| MPBC1 | H05 | 0.1181 | 1635 | 82 | 0.05015 |
| MPBC1 | H10 | 0.1034 | 1635 | 164 | 0.1003 |
| MPBC2 | H05 | 0.09879 | 236 | 12 | 0.05085 |
| MPBC2 | H10 | 0.08477 | 236 | 24 | 0.1017 |

仅同片同输入组件的真实μm质心距离；50/200均累计圆。共享任何原始bin的整对T32–Epi64剔除，原分数不作分摊。U_r需要T_eligible、原supplied circle coverage≥0.5、至少1个独立可评分Epi64。所有观察/null/曲线沿用同一来源排除。

| section | radius_um | T_eligible | U_r | FOXP3_in_U | strict_in_U |
| --- | --- | --- | --- | --- | --- |
| Cervilla | 50 | 280 | 43 | 9 | 0 |
| Cervilla | 200 | 280 | 261 | 36 | 2 |
| D1_Johnson | 50 | 7103 | 2640 | 462 | 0 |
| D1_Johnson | 200 | 7103 | 5156 | 919 | 17 |
| UC_SCC | 50 | 3579 | 418 | 31 | 0 |
| UC_SCC | 200 | 3579 | 2214 | 165 | 2 |
| MPBC1 | 50 | 547 | 37 | 1 | 0 |
| MPBC1 | 200 | 547 | 349 | 24 | 0 |
| MPBC2 | 50 | 2270 | 46 | 3 | 0 |
| MPBC2 | 200 | 2270 | 396 | 32 | 0 |

50μm首拒绝原因：

| section | reason | n | denominator |
| --- | --- | --- | --- |
| Cervilla | no_qualified_Epi_within_radius | 204 | 280 |
| Cervilla | eligible | 43 | 280 |
| Cervilla | all_pairs_shared_source | 33 | 280 |
| D1_Johnson | no_qualified_Epi_within_radius | 3108 | 7103 |
| D1_Johnson | eligible | 2640 | 7103 |
| D1_Johnson | all_pairs_shared_source | 1355 | 7103 |
| UC_SCC | no_qualified_Epi_within_radius | 2727 | 3579 |
| UC_SCC | all_pairs_shared_source | 433 | 3579 |
| UC_SCC | eligible | 418 | 3579 |
| UC_SCC | coverage_below0.5 | 1 | 3579 |
| MPBC1 | no_qualified_Epi_within_radius | 503 | 547 |
| MPBC1 | eligible | 37 | 547 |
| MPBC1 | all_pairs_shared_source | 7 | 547 |
| MPBC2 | no_qualified_Epi_within_radius | 2206 | 2270 |
| MPBC2 | eligible | 46 | 2270 |
| MPBC2 | all_pairs_shared_source | 18 | 2270 |

粒度与来源排除确实限制近距离分析：Cervilla有33个中心的全部候选对因共享来源而失去资格，D1为1355，UC_SCC433，MPBC1为7，MPBC2为18。另有大量位置本来就没有达到Epi64测量资格的邻居。不能将这些位置放入非high对照。

## 主FE及历史标签比较

FE分母按唯一T32中心计数，PR比较high与有测量但非high组的标签阳性率；OR不加0.5。完整a/b/c/d、0/Inf/NA及其分母在`contingency_tables.tsv`和`enrichment_effects.tsv`。

| section | radius_um | N | label_positive | a | near_high_denominator | FE | FE_low | FE_high | P | q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 50 | 43 | 9 | 0 | 1 | 0 | NA | NA | 1 | 1 |
| Cervilla | 200 | 261 | 36 | 10 | 56 | 1.295 | 0.6371 | 1.898 | 0.814 | 1 |
| D1_Johnson | 50 | 2640 | 462 | 18 | 137 | 0.7508 | 0.351 | 1.213 | 0.3969 | 1 |
| D1_Johnson | 200 | 5156 | 919 | 427 | 2324 | 1.031 | 0.943 | 1.112 | 0.2229 | 0.8916 |
| UC_SCC | 50 | 418 | 31 | 1 | 12 | 1.124 | 0 | 4.355 | 1 | 1 |
| UC_SCC | 200 | 2214 | 165 | 38 | 458 | 1.113 | 0.7867 | 1.445 | 0.2082 | 0.8916 |
| MPBC1 | 50 | 37 | 1 | 0 | 4 | 0 | NA | NA | NA | NA |
| MPBC1 | 200 | 349 | 24 | 8 | 109 | 1.067 | 0.4185 | 1.88 | 0.8099 | 1 |
| MPBC2 | 50 | 46 | 3 | 0 | 0 | NA | NA | NA | NA | NA |
| MPBC2 | 200 | 396 | 32 | 2 | 28 | 0.8839 | 0 | 2.785 | 1 | 1 |

以上区间为400μm空间块重采样，NA包括样本支持不足、过多退化抽样。800μm敏感性中，全部五片200μm FOXP3 FE区间仍跨1。主FE的双侧P基于9999次同组件/检测量/T信号/上皮机会/coverage分层内标签置换。分层不使用MP10高低或结果选择；每条标签的可交换数和不可交换比例单独保存。该零模型不是完整空间自相关保持模型，因此区块/热点检查不可省略。

严格标签在Cervilla/UC_SCC的200μm有效总体仅各2个，MPBC1/2为0，推断不足。D1有17个，H05 FE≈0.522；这也没有经辅助家族校正的明确证据。H10及strict共30个辅助FE计划行完整保留，没有依据方向替换H05主定义。

## 邻接边、距离和混合

| section | radius_um | O | null_mean | null_sd | Z | P | q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 50 | 0 | 0.1486 | 0.3557 | -0.4178 | 1 | 1 |
| Cervilla | 200 | 18 | 20.56 | 5.947 | -0.4311 | 0.7079 | 1 |
| D1_Johnson | 50 | 19 | 23.26 | 4.503 | -0.9454 | 0.3773 | 1 |
| D1_Johnson | 200 | 958 | 942 | 46.88 | 0.3407 | 0.7431 | 1 |
| UC_SCC | 50 | 1 | 0.8705 | 0.874 | 0.1482 | 1 | 1 |
| UC_SCC | 200 | 57 | 49.11 | 10.53 | 0.7498 | 0.474 | 1 |
| MPBC1 | 50 | 0 | 0 | 0 | NA | NA | NA |
| MPBC1 | 200 | 9 | 8.959 | 2.862 | 0.01443 | 1 | 1 |
| MPBC2 | 50 | 0 | 0 | 0 | NA | NA | NA |
| MPBC2 | 200 | 2 | 2.743 | 1.77 | -0.4198 | 0.7981 | 1 |

边数对有更多上皮邻居的中心赋予更多配对机会，不能混称FE。全部40项辅助边检验共同BH；没有明确校正后信号。零方差或可交换支持不足时Z/P=NA，观察数保留。

来源排除后的累计交叉邻近曲线固定U200，0–200μm每10μm，报告每阳性中心的高上皮邻居数及同mask限制置换基线；不在曲线峰值追加新检验。200μm观察/随机均值边比为：

| section | N_fixed_U200 | n_positive | observed_edges | null_mean_edges | observed_over_null |
| --- | --- | --- | --- | --- | --- |
| Cervilla | 261 | 36 | 18 | 20.56 | 0.8753 |
| D1_Johnson | 5156 | 919 | 958 | 942 | 1.017 |
| UC_SCC | 2214 | 165 | 57 | 49.11 | 1.161 |
| MPBC1 | 349 | 24 | 9 | 8.959 | 1.005 |
| MPBC2 | 396 | 32 | 2 | 2.743 | 0.7291 |

曲线包络是逐距离点的条件随机包络，不是同时包络，不能据个别越界点宣称新的特异距离。主50μm总体与曲线固定U200不同，数字不应直接替换。

标准形式K/L仅对前三片已知mask、supplied in-tissue观察窗做同源全局不重叠子集的8μm像素border估计。D1因此从2884个Epi64减少至1001个，子集损失很大。保留各组件及不足点数的NA；这是自实现的描述性homogeneous border估计，未解决组织非均匀强度，不是未经修改的原生Kcross显著性检验。MPBC1/2缺真实组织窗，标准K/L未执行并有明确NA；没有用包围矩形替代。适配曲线五片均完成。

双向MS/NMS全部80个组合已计算，包含观察范围内无边点。A-A为无序对、A-B计一次；零A-A分母NA，合法零跨边才为0。MS/NMS有不同的总体(coverage合格T及其组件内全部合格Epi)，不冒充FE的U_r；999次同规则随机基线仅供描述，没有额外P。NMS=1不是非均匀组织中自动的无关联基准。见`mixing_scores.tsv`与双向图。

## 连续程序及稳健性

| section | radius_um | N | spearman | adjusted_rank_beta | low | high |
| --- | --- | --- | --- | --- | --- | --- |
| Cervilla | 50 | 43 | 0.1639 | 0.1421 | NA | NA |
| Cervilla | 200 | 261 | 0.04151 | -0.005822 | -0.09757 | 0.0746 |
| D1_Johnson | 50 | 2640 | -0.09099 | -0.02353 | -0.06558 | 0.01914 |
| D1_Johnson | 200 | 5156 | -0.1358 | -0.03601 | -0.06805 | -0.0009054 |
| UC_SCC | 50 | 418 | 0.09252 | 0.04195 | -0.06359 | 0.1299 |
| UC_SCC | 200 | 2214 | 0.02065 | 0.0399 | -0.01951 | 0.09358 |
| MPBC1 | 50 | 37 | -0.1074 | -0.09824 | NA | NA |
| MPBC1 | 200 | 349 | 0.009137 | 0.02752 | -0.08728 | 0.1469 |
| MPBC2 | 50 | 46 | -0.03087 | 0.08062 | NA | NA |
| MPBC2 | 200 | 396 | 0.01484 | 0.0084 | -0.09313 | 0.1165 |

原始Spearman与固定秩标准化调整均保留。调整候选顺序为T-lineage、T侧UMI、Epi侧深度、独立Epi区域数、coverage和原T区域背景；顺序残差方差≤0.1的共线列按统一规则不纳入。机会数是本轮明确加入的诊断调整，因此不把系数冒充v1原系数。实际设计矩阵保存。D1的200μm负向效应在800μm区间约[−0.0697,−0.0033]，去除原固定候选ROI后β≈−0.034；只有片内、探索性支持，未跨片重复。

200μm逐个阳性热点排除后的FE范围：

| section | FE_min | FE_max | max_removed | hotspots |
| --- | --- | --- | --- | --- |
| Cervilla | 0.7229 | 1.36 | 8 | 28 |
| D1_Johnson | 1.018 | 1.041 | 14 | 501 |
| UC_SCC | 1.028 | 1.145 | 5 | 127 |
| MPBC1 | 0.8845 | 1.11 | 2 | 22 |
| MPBC2 | 0 | 1 | 4 | 23 |

热点是预定64μm距离连接的同组件阳性区域组；仅删除组内中心，剩余表达/邻域不变。另有原固定候选ROI合并范围内全部T中心排除，以及每个400/800μm空间块leave-one-out。Cervilla FE可从0.72至1.36，MPBC2从0至1，说明局部正例少、影响较大；D1/UC_SCC范围更窄但总体仍接近1。没有按影响大小删除对象后选择最终模型。

## 固定ROI与局部验证建议

原24个ROI ID/中心/候选对照身份全部保留，以512μm视野显示50/200圆。新叠图同时展示有测量的Epi64及H05，统计表明确哪些配对被来源排除。原v2图册继续提供H&E/原始多基因/成员支持图，逐ROI链接和解释在`ROI_interpretation.tsv`，没有新挑选最漂亮的重叠区域。

Cervilla两个相邻候选不是独立重复：candidate_1有FOXP3=1、IL2RA=1且strict；candidate_2 FOXP3=0、IKZF2=1而strict。两者50μm均无独立可评分Epi64；200μm各有相同的2个非H05 Epi64。这支持核查局部身份与测量来源，不支持本轮H05邻近富集。

有限后续优先级(本轮未启动)：D1 candidate_1及control_2可用于区分“FOXP3+联合严格程序”与“FOXP3+但不strict”，有原高分辨率图与较多全片测量支持；Cervilla相邻候选与control_1适合保留旧cell来源差异的局部核查；UC_SCC candidate_3核心多基因可检测，但200μm没有独立合格Epi64，主要用于身份/测量核验。MPBC候选缺H&E和原组织mask，只支持表达区域描述，不宜先进行细胞位置归因。优先级依据身份/覆盖与材料可用性，不依据富集P。

## 解释边界与完成状态

本轮在可评价T区域作背景，并控制部分T信号、检测深度和上皮机会，因此检验的是标签/程序相对于一般T区域的偏好；仍不能完整区分活化T、稀疏FOXP3和真实Treg细胞。连续分数不是比例，FOXP3标签不是细胞鉴定，空间混合不是通讯。

已完成：五片定义审核、冻结标签与MP10核验、校准分母修复、40组FE/PR/OR与40组边统计、双向80组MS/NMS、连续与空间稳健性、五片适配距离曲线、已知窗的描述性K/L尝试、24固定ROI新叠图、轻量加载及重读。标准K/L在MPBC无窗及小类型/内部点不足组合不可估计；strict多数分支推断不足；50μm的部分主效应/P/区间不可估计，均保留计划行和原因。没有因缺输入停止独立模块。

没有重分割、重聚合counts、重跑BANKSY或改变历史MP10/Treg分数；不需要新增计数守恒试验，复用原已验证聚合及来源映射。原始文件mtime/尺寸、逐ID数值、来源排除、置换分母、无序配对因子2、退化抽样和输出重读均核验。没有复制全表达矩阵或上传发布。

切片患者关系仍未知/未建立独立性，MPBC tissue_original=NA；不作五患者总体P。结论为：**未检出跨片一致的正向富集；D1有小幅负向连续线索；多个近距离/严格标签分支信息不足。** 有理由进行有限局部身份核验，没有证据据此开展全片细胞重建或宣称通讯增强。

文件入口：`code/load_result.py`、`code/load_result.R`；复跑`rerun.sh`，仅从已存表重绘`replot.sh`。主要表格位于tables，分片紧凑区域表位于各section目录。原counts、T32/Epi64成员及配对依赖v1路径在manifest/loader中明确。资源和磁盘占用见`tables/resource_and_delivery.tsv`。
