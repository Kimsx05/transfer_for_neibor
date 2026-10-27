# MP10 空间细胞通讯试点报告

完成日期：2026-10-27（北京时间）。本轮在同一批两张完整Visium切片上完成SpatialCellChat与COMMOT，另对首片完成一次较宽作用距离敏感性。旧scRNA与空间输入均只读。主要输出为全部纳入LR的总强度，NECTIN2–TIGIT未被拿来替代总网络。

## 主要结果（探索性）

Epi→Treg：在4个方法×切片组合中，主强度High−Low有0项为正、4项为负。排名及谱系份额逐片变化，不能把总强度或份额称为免疫抑制强度，也不能由其推断功能和因果。

Treg→Epi：在4个方法×切片组合中，主强度High−Low有0项为正、4项为负。排名及谱系份额逐片变化，不能把总强度或份额称为免疫抑制强度，也不能由其推断功能和因果。

Epi→EC/基质：两片与两算法方向均增加的可辨标签为iCAF；均下降的为Arterial EC, Lymphatic EC, Venous EC；存在片间或方法间不一致的为Capillary EC, Matrix_CAF, Pericyte, SMC, rCAF, vCAF。这四项比较共享输入/数据库因素，不是四次独立实验。

EC/基质→Epi：两片与两算法方向均增加的可辨标签为iCAF；均下降的为Arterial EC, Lymphatic EC, Venous EC；存在片间或方法间不一致的为Capillary EC, Matrix_CAF, Pericyte, SMC, rCAF, vCAF。这四项比较共享输入/数据库因素，不是四次独立实验。

主强度与流量总和可能方向相反，尤其应看Treg的逐片F与分母。只有一位来源核验患者加一份未核验Inhouse记录，结果是有限切片探索；没有患者级显著性结论。

## 试点与可辨识范围

| section_id | dataset | n_spots | n_epi_high | n_epi_low | patient_id_status |
| --- | --- | --- | --- | --- | --- |
| GSE171351_A1 | GSE171351 | 768 | 244 | 244 | VERIFIED_WITHIN_STUDY |
| Inhousedata_T4469479 | Inhousedata__UNVERIFIED | 3917 | 923 | 923 | USER_ASSUMED_SAMPLE_UNIT_PATIENT_ID_UNVERIFIED |


两片分别为768和3917个混合spot，均20/20 MP10基因可测。GSE171351有来源核验患者标识；Inhouse仅为未核实身份的记录单位，因此不能称为两个已核验患者，也不作患者级显著性检验。选片在通讯之前依据完整性和每数据集最小计算规模；不按Treg结果选择。

每片保留完整背景。使用accepted cell2location 30类型后验均值除以全部建模类型总量的Q，未做0.1截断或dominant硬标签。三个原参考低可辨识T/NK列合并为T_NK_unresolved，保留其质量。Treg仍单列；Pericyte/SMC遵循空间原映射保留在Mural，未改成CAF。Tip EC和Cycling EC未独立建模，主圈图明确NE。其余细群多数仍属于INTERPRET_WITH_CAUTION，不等于已确认单细胞身份。

H/L只在既有Epi比例≥0.25 anchors内按冻结MP10中位数严格划分，等值/其他spot上皮成分留作EpiOther背景。每片H/L共用一个位置网络和一份TME Q。Epi_H/L准确含义是MP10高/低关联spot中的上皮成分，未分辨同spot内两种上皮细胞。

## 模型、总强度与距离

| method | section_id | n_spots | n_LR_selected | n_LR_cross_supported | elapsed_seconds | peak_single_process_rss_GiB | status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | 768 | 675 | 675 | 1118.4 | 1.0161 | complete |
| COMMOT | Inhousedata_T4469479 | 3917 | 1212 | 1212 | 2729.9 | 1.9589 | complete |
| SpatialCellChat | GSE171351_A1 | 768 | 1037 | 593 | 404.45 | 0.90517 | complete |
| SpatialCellChat | Inhousedata_T4469479 | 3917 | 1829 | 1801 | 1546.5 | 1.9623 | complete |


SpatialCellChat主指标total_mean_strength为各通过原生位置QC的LR跨spot矩阵按官方Visium Avg helper计算F/K后求和；K为该LR有链接且两类型有支持的spot-pair数。没有group P二次筛选。COMMOT在完整片中对全部纳入LR联合求竞争性最优传输，主G=QᵀTcrossQ/[sumQa·sumQb−sum(QaQb)]，分母包含跨spot零边。每种方法均另存total_flow和明确分母；major与fine独立聚合，绝不相加fine均值冒充major。

实测counts仅做每spot总量1万+log1p一次，未用c2l比例或MP10替代LR表达。human库覆盖Secreted、ECM、Contact三类蛋白信号，方法实际库版本和筛选后LR数不同；SCC按原DB交互记录计数，COMMOT精确ligand/receptor对去重计数，展开复合物后的交集见cross_method_LR_scope.tsv。主分泌/ECM距离为2×dNN=200μm，接触为1.05×dNN=105μm；采用可靠Visium阵列索引及标称100μm间距，非误用图像倍率。额外1e−6μm仅处理浮点边界。邻spot接触信号不证明膜接触。

所有模型先保留位置对角项推断，汇总副本中再去除同spot项并单列原贡献。这不等于重新求解禁止同位置通讯的模型，也不能把Visium同spot项称为准确的自分泌。

两算法的绝对数值、分母和内部筛选不同，也与scRNA W不同，不能直接相减、平均或比较大小。共享数据库/去卷积时，方向一致并非独立实验验证。

## Major结果

### SpatialCellChat

| section_id | node | direction | strength_low | strength_high | delta_strength | share_all_TME_low | share_all_TME_high |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE171351_A1 | B | epi_to_TME | 1.6837 | 1.7209 | 0.037183 | 0.10601 | 0.10621 |
| GSE171351_A1 | Plasma | epi_to_TME | 0.265 | 0.29511 | 0.030109 | 0.016684 | 0.018214 |
| GSE171351_A1 | T/NK | epi_to_TME | 2.166 | 2.252 | 0.085973 | 0.13637 | 0.13899 |
| GSE171351_A1 | Myeloid | epi_to_TME | 3.3014 | 3.2875 | -0.013906 | 0.20785 | 0.2029 |
| GSE171351_A1 | Fibroblast | epi_to_TME | 2.8905 | 3.0975 | 0.207 | 0.18199 | 0.19118 |
| GSE171351_A1 | Mural | epi_to_TME | 2.7721 | 2.7675 | -0.0046739 | 0.17453 | 0.1708 |
| GSE171351_A1 | Endothelial | epi_to_TME | 2.6202 | 2.5959 | -0.024267 | 0.16497 | 0.16022 |
| GSE171351_A1 | Mast | epi_to_TME | 0.18422 | 0.18601 | 0.0017939 | 0.011598 | 0.01148 |
| GSE171351_A1 | B | TME_to_epi | 1.6897 | 1.7364 | 0.046688 | 0.10519 | 0.10597 |
| GSE171351_A1 | Plasma | TME_to_epi | 0.26605 | 0.2964 | 0.030346 | 0.016563 | 0.018089 |
| GSE171351_A1 | T/NK | TME_to_epi | 2.1882 | 2.2683 | 0.080045 | 0.13622 | 0.13843 |
| GSE171351_A1 | Myeloid | TME_to_epi | 3.32 | 3.2939 | -0.026094 | 0.20668 | 0.20103 |
| GSE171351_A1 | Fibroblast | TME_to_epi | 2.9469 | 3.1471 | 0.20015 | 0.18346 | 0.19206 |
| GSE171351_A1 | Mural | TME_to_epi | 2.8366 | 2.8245 | -0.012033 | 0.17659 | 0.17238 |
| GSE171351_A1 | Endothelial | TME_to_epi | 2.6285 | 2.6276 | -0.00090455 | 0.16364 | 0.16036 |
| GSE171351_A1 | Mast | TME_to_epi | 0.18734 | 0.19136 | 0.0040219 | 0.011662 | 0.011679 |
| Inhousedata_T4469479 | B | epi_to_TME | 1.7827 | 1.2993 | -0.48333 | 0.21218 | 0.19715 |
| Inhousedata_T4469479 | Plasma | epi_to_TME | 0.29297 | 0.26379 | -0.029178 | 0.034871 | 0.040025 |
| Inhousedata_T4469479 | T/NK | epi_to_TME | 1.4921 | 1.0489 | -0.44321 | 0.1776 | 0.15915 |
| Inhousedata_T4469479 | Myeloid | epi_to_TME | 1.9293 | 1.7175 | -0.21176 | 0.22964 | 0.2606 |
| Inhousedata_T4469479 | Fibroblast | epi_to_TME | 1.678 | 1.2346 | -0.4434 | 0.19972 | 0.18732 |
| Inhousedata_T4469479 | Mural | epi_to_TME | 0.38813 | 0.32376 | -0.064374 | 0.046198 | 0.049124 |
| Inhousedata_T4469479 | Endothelial | epi_to_TME | 0.54199 | 0.42351 | -0.11849 | 0.064511 | 0.064258 |
| Inhousedata_T4469479 | Mast | epi_to_TME | 0.29642 | 0.27927 | -0.017147 | 0.035281 | 0.042373 |
| Inhousedata_T4469479 | B | TME_to_epi | 1.8052 | 1.34 | -0.46527 | 0.20938 | 0.1968 |
| Inhousedata_T4469479 | Plasma | TME_to_epi | 0.29267 | 0.26395 | -0.028726 | 0.033946 | 0.038766 |
| Inhousedata_T4469479 | T/NK | TME_to_epi | 1.5567 | 1.1019 | -0.4548 | 0.18056 | 0.16184 |
| Inhousedata_T4469479 | Myeloid | TME_to_epi | 1.9545 | 1.7501 | -0.20444 | 0.2267 | 0.25703 |
| Inhousedata_T4469479 | Fibroblast | TME_to_epi | 1.754 | 1.3019 | -0.45206 | 0.20344 | 0.19122 |
| Inhousedata_T4469479 | Mural | TME_to_epi | 0.40429 | 0.33505 | -0.06924 | 0.046892 | 0.049209 |
| Inhousedata_T4469479 | Endothelial | TME_to_epi | 0.55692 | 0.4352 | -0.12172 | 0.064595 | 0.063918 |
| Inhousedata_T4469479 | Mast | TME_to_epi | 0.29739 | 0.28059 | -0.016806 | 0.034493 | 0.04121 |


### COMMOT

| section_id | node | direction | strength_low | strength_high | delta_strength | share_all_TME_low | share_all_TME_high |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE171351_A1 | B | epi_to_TME | 0.026109 | 0.02684 | 0.00073166 | 0.13144 | 0.13043 |
| GSE171351_A1 | Plasma | epi_to_TME | 0.022436 | 0.025518 | 0.003082 | 0.11295 | 0.124 |
| GSE171351_A1 | T/NK | epi_to_TME | 0.026062 | 0.02736 | 0.0012978 | 0.1312 | 0.13295 |
| GSE171351_A1 | Myeloid | epi_to_TME | 0.028835 | 0.028692 | -0.00014307 | 0.14517 | 0.13942 |
| GSE171351_A1 | Fibroblast | epi_to_TME | 0.024959 | 0.027142 | 0.0021824 | 0.12565 | 0.13189 |
| GSE171351_A1 | Mural | epi_to_TME | 0.020357 | 0.02028 | -7.7048e-05 | 0.10248 | 0.098546 |
| GSE171351_A1 | Endothelial | epi_to_TME | 0.028065 | 0.027721 | -0.00034405 | 0.14129 | 0.1347 |
| GSE171351_A1 | Mast | epi_to_TME | 0.021811 | 0.022238 | 0.0004262 | 0.10981 | 0.10806 |
| GSE171351_A1 | B | TME_to_epi | 0.026534 | 0.027775 | 0.001241 | 0.13078 | 0.1298 |
| GSE171351_A1 | Plasma | TME_to_epi | 0.022975 | 0.026516 | 0.0035405 | 0.11324 | 0.12392 |
| GSE171351_A1 | T/NK | TME_to_epi | 0.026527 | 0.028163 | 0.0016361 | 0.13075 | 0.13161 |
| GSE171351_A1 | Myeloid | TME_to_epi | 0.029272 | 0.029526 | 0.00025441 | 0.14428 | 0.13799 |
| GSE171351_A1 | Fibroblast | TME_to_epi | 0.025565 | 0.028106 | 0.0025412 | 0.12601 | 0.13135 |
| GSE171351_A1 | Mural | TME_to_epi | 0.021185 | 0.021751 | 0.00056659 | 0.10442 | 0.10165 |
| GSE171351_A1 | Endothelial | TME_to_epi | 0.028348 | 0.028441 | 9.2646e-05 | 0.13973 | 0.13291 |
| GSE171351_A1 | Mast | TME_to_epi | 0.022477 | 0.023703 | 0.0012266 | 0.11079 | 0.11077 |
| Inhousedata_T4469479 | B | epi_to_TME | 0.0071096 | 0.0055436 | -0.001566 | 0.16914 | 0.15624 |
| Inhousedata_T4469479 | Plasma | epi_to_TME | 0.005351 | 0.0050308 | -0.00032017 | 0.12731 | 0.14179 |
| Inhousedata_T4469479 | T/NK | epi_to_TME | 0.0042423 | 0.0030548 | -0.0011875 | 0.10093 | 0.086095 |
| Inhousedata_T4469479 | Myeloid | epi_to_TME | 0.0041916 | 0.0038993 | -0.00029232 | 0.099723 | 0.1099 |
| Inhousedata_T4469479 | Fibroblast | epi_to_TME | 0.0061469 | 0.0046932 | -0.0014537 | 0.14624 | 0.13227 |
| Inhousedata_T4469479 | Mural | epi_to_TME | 0.0045374 | 0.0039249 | -0.00061246 | 0.10795 | 0.11062 |
| Inhousedata_T4469479 | Endothelial | epi_to_TME | 0.0049148 | 0.0039372 | -0.00097755 | 0.11693 | 0.11097 |
| Inhousedata_T4469479 | Mast | epi_to_TME | 0.0055393 | 0.0053977 | -0.00014165 | 0.13179 | 0.15213 |
| Inhousedata_T4469479 | B | TME_to_epi | 0.0075175 | 0.0061166 | -0.0014009 | 0.16611 | 0.15743 |
| Inhousedata_T4469479 | Plasma | TME_to_epi | 0.0056166 | 0.0052793 | -0.0003373 | 0.12411 | 0.13588 |
| Inhousedata_T4469479 | T/NK | TME_to_epi | 0.0047763 | 0.0035149 | -0.0012614 | 0.10554 | 0.09047 |
| Inhousedata_T4469479 | Myeloid | TME_to_epi | 0.0043772 | 0.0042 | -0.00017723 | 0.096723 | 0.1081 |
| Inhousedata_T4469479 | Fibroblast | TME_to_epi | 0.0068374 | 0.0053509 | -0.0014865 | 0.15108 | 0.13772 |
| Inhousedata_T4469479 | Mural | TME_to_epi | 0.0051407 | 0.004436 | -0.00070468 | 0.11359 | 0.11418 |
| Inhousedata_T4469479 | Endothelial | TME_to_epi | 0.005307 | 0.0042989 | -0.0010082 | 0.11727 | 0.11065 |
| Inhousedata_T4469479 | Mast | TME_to_epi | 0.0056824 | 0.0056555 | -2.6831e-05 | 0.12556 | 0.14557 |


## T/NK内的Treg

| method | section_id | direction | strength_low | strength_high | delta_strength | rank_strength_high | rank_delta_strength | share_all_TME_low | share_all_TME_high | share_within_lineage_low | share_within_lineage_high |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | epi_to_TME | 0.030508 | 0.030264 | -0.00024473 | 1 | 5 | 0.041522 | 0.040938 | 0.23493 | 0.22146 |
| COMMOT | GSE171351_A1 | TME_to_epi | 0.030806 | 0.030804 | -2.7378e-06 | 1 | 5 | 0.041238 | 0.040377 | 0.23306 | 0.21883 |
| COMMOT | Inhousedata_T4469479 | epi_to_TME | 0.0041036 | 0.002986 | -0.0011176 | 4 | 3 | 0.031001 | 0.026338 | 0.18642 | 0.18298 |
| COMMOT | Inhousedata_T4469479 | TME_to_epi | 0.0046618 | 0.0034863 | -0.0011756 | 4 | 3 | 0.032593 | 0.027965 | 0.18934 | 0.18672 |
| SpatialCellChat | GSE171351_A1 | epi_to_TME | 0.38329 | 0.37941 | -0.0038889 | 2 | 5 | 0.024132 | 0.023417 | 0.17696 | 0.16848 |
| SpatialCellChat | GSE171351_A1 | TME_to_epi | 0.38676 | 0.38218 | -0.0045785 | 2 | 5 | 0.024077 | 0.023324 | 0.17675 | 0.16849 |
| SpatialCellChat | Inhousedata_T4469479 | epi_to_TME | 0.20451 | 0.14742 | -0.057098 | 5 | 3 | 0.024342 | 0.022367 | 0.13706 | 0.14054 |
| SpatialCellChat | Inhousedata_T4469479 | TME_to_epi | 0.21597 | 0.15674 | -0.05924 | 5 | 3 | 0.02505 | 0.02302 | 0.13874 | 0.14224 |


COMMOT / GSE171351_A1 / epi_to_TME：High强度在T/NK内排名1/5，增量排名5/5；谱系内份额23.49%→22.15%，全TME份额4.15%→4.09%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

COMMOT / GSE171351_A1 / TME_to_epi：High强度在T/NK内排名1/5，增量排名5/5；谱系内份额23.31%→21.88%，全TME份额4.12%→4.04%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

COMMOT / Inhousedata_T4469479 / epi_to_TME：High强度在T/NK内排名4/5，增量排名3/5；谱系内份额18.64%→18.30%，全TME份额3.10%→2.63%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

COMMOT / Inhousedata_T4469479 / TME_to_epi：High强度在T/NK内排名4/5，增量排名3/5；谱系内份额18.93%→18.67%，全TME份额3.26%→2.80%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

SpatialCellChat / GSE171351_A1 / epi_to_TME：High强度在T/NK内排名2/5，增量排名5/5；谱系内份额17.70%→16.85%，全TME份额2.41%→2.34%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

SpatialCellChat / GSE171351_A1 / TME_to_epi：High强度在T/NK内排名2/5，增量排名5/5；谱系内份额17.67%→16.85%，全TME份额2.41%→2.33%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

SpatialCellChat / Inhousedata_T4469479 / epi_to_TME：High强度在T/NK内排名5/5，增量排名3/5；谱系内份额13.71%→14.05%，全TME份额2.43%→2.24%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

SpatialCellChat / Inhousedata_T4469479 / TME_to_epi：High强度在T/NK内排名5/5，增量排名3/5；谱系内份额13.87%→14.22%，全TME份额2.50%→2.30%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。

SpatialCellChat / Inhousedata_T4469479 / epi_to_TME的主强度与flow方向相反：主强度0.204514→0.147416，flow 106.585→120.582。这是归一化分母与流量口径不同，不能把主mean下降写成所有通讯总量减少。

SpatialCellChat / Inhousedata_T4469479 / TME_to_epi的主强度与flow方向相反：主强度0.215975→0.156735，flow 115.411→133.793。这是归一化分母与流量口径不同，不能把主mean下降写成所有通讯总量减少。

SpatialCellChat总mean为sum_l(F_l/K_l)，一般不等于sum_l F_l / sum_l K_l；表中summed_LR_K仅作链接覆盖汇总，不是主mean的统一分母。COMMOT的G则有其单一跨位置权重分母D，两者不能混用。SCC的K计有正Q支持的保留链接数，不直接除以目标类型Q总量；COMMOT的D则按Q加权的全部跨位置配对定义。因此稀有类型在两方法中的绝对排名可以不同，不能把某算法较高名次直接解释成真实细胞通讯更强。

T/NK对照保留CD8 Exhausted、Tfh/CXCL13-like、NK及原低可辨识质量合成的T_NK_unresolved；后者是背景集合，不应与单一生物亚型等同。完整双方向强度、增量、两个分母份额见各片comparison_T_NK面板。Treg即使从未占某spot第一也没有被删掉。

| method | section_id | node | direction | strength_high | delta_strength | rank_strength_high | delta_share_within_lineage |
| --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | CD8 Exhausted | epi_to_TME | 0.02581 | 0.0019076 | 5 | 0.0048136 |
| COMMOT | GSE171351_A1 | NK | epi_to_TME | 0.026381 | 0.00083077 | 4 | -0.0036967 |
| COMMOT | GSE171351_A1 | Tfh/CXCL13-like | epi_to_TME | 0.02681 | 0.003301 | 3 | 0.015161 |
| COMMOT | GSE171351_A1 | Treg | epi_to_TME | 0.030264 | -0.00024473 | 1 | -0.013464 |
| COMMOT | GSE171351_A1 | T_NK_unresolved | epi_to_TME | 0.027389 | 0.00099544 | 2 | -0.0028142 |
| COMMOT | GSE171351_A1 | CD8 Exhausted | TME_to_epi | 0.026704 | 0.0023254 | 5 | 0.0052746 |
| COMMOT | GSE171351_A1 | NK | TME_to_epi | 0.027438 | 0.0014196 | 4 | -0.0019171 |
| COMMOT | GSE171351_A1 | Tfh/CXCL13-like | TME_to_epi | 0.027667 | 0.0035411 | 3 | 0.014027 |
| COMMOT | GSE171351_A1 | Treg | TME_to_epi | 0.030804 | -2.7378e-06 | 1 | -0.01423 |
| COMMOT | GSE171351_A1 | T_NK_unresolved | TME_to_epi | 0.028155 | 0.0012997 | 2 | -0.0031546 |
| COMMOT | Inhousedata_T4469479 | CD8 Exhausted | epi_to_TME | 0.0035713 | -0.0013087 | 2 | -0.0028448 |
| COMMOT | Inhousedata_T4469479 | NK | epi_to_TME | 0.0035979 | -0.00094643 | 1 | 0.014037 |
| COMMOT | Inhousedata_T4469479 | Tfh/CXCL13-like | epi_to_TME | 0.0035059 | -0.0010434 | 3 | 0.0081702 |
| COMMOT | Inhousedata_T4469479 | Treg | epi_to_TME | 0.002986 | -0.0011176 | 4 | -0.0034393 |
| COMMOT | Inhousedata_T4469479 | T_NK_unresolved | epi_to_TME | 0.0026577 | -0.0012779 | 5 | -0.015923 |
| COMMOT | Inhousedata_T4469479 | CD8 Exhausted | TME_to_epi | 0.004089 | -0.0013653 | 1 | -0.0025195 |
| COMMOT | Inhousedata_T4469479 | NK | TME_to_epi | 0.0040209 | -0.00096315 | 2 | 0.01293 |
| COMMOT | Inhousedata_T4469479 | Tfh/CXCL13-like | TME_to_epi | 0.0039803 | -0.0010516 | 3 | 0.0088146 |
| COMMOT | Inhousedata_T4469479 | Treg | TME_to_epi | 0.0034863 | -0.0011756 | 4 | -0.0026168 |
| COMMOT | Inhousedata_T4469479 | T_NK_unresolved | TME_to_epi | 0.0030949 | -0.0013953 | 5 | -0.016609 |
| SpatialCellChat | GSE171351_A1 | CD8 Exhausted | epi_to_TME | 0.31947 | 0.01991 | 4 | 0.0035612 |
| SpatialCellChat | GSE171351_A1 | NK | epi_to_TME | 0.28249 | 0.0099135 | 5 | -0.00040217 |
| SpatialCellChat | GSE171351_A1 | Tfh/CXCL13-like | epi_to_TME | 0.34617 | 0.033604 | 3 | 0.0094127 |
| SpatialCellChat | GSE171351_A1 | Treg | epi_to_TME | 0.37941 | -0.0038889 | 2 | -0.0084825 |
| SpatialCellChat | GSE171351_A1 | T_NK_unresolved | epi_to_TME | 0.92445 | 0.026435 | 1 | -0.0040892 |
| SpatialCellChat | GSE171351_A1 | CD8 Exhausted | TME_to_epi | 0.32095 | 0.018395 | 4 | 0.0032306 |
| SpatialCellChat | GSE171351_A1 | NK | TME_to_epi | 0.28508 | 0.0086757 | 5 | -0.00063273 |
| SpatialCellChat | GSE171351_A1 | Tfh/CXCL13-like | TME_to_epi | 0.34946 | 0.031635 | 3 | 0.0088212 |
| SpatialCellChat | GSE171351_A1 | Treg | TME_to_epi | 0.38218 | -0.0045785 | 2 | -0.0082558 |
| SpatialCellChat | GSE171351_A1 | T_NK_unresolved | TME_to_epi | 0.9306 | 0.025918 | 1 | -0.0031633 |
| SpatialCellChat | Inhousedata_T4469479 | CD8 Exhausted | epi_to_TME | 0.15852 | -0.067682 | 3 | -0.00046921 |
| SpatialCellChat | Inhousedata_T4469479 | NK | epi_to_TME | 0.15231 | -0.046141 | 4 | 0.012207 |
| SpatialCellChat | Inhousedata_T4469479 | Tfh/CXCL13-like | epi_to_TME | 0.16472 | -0.054828 | 2 | 0.0098994 |
| SpatialCellChat | Inhousedata_T4469479 | Treg | epi_to_TME | 0.14742 | -0.057098 | 5 | 0.0034785 |
| SpatialCellChat | Inhousedata_T4469479 | T_NK_unresolved | epi_to_TME | 0.42596 | -0.21746 | 1 | -0.025116 |
| SpatialCellChat | Inhousedata_T4469479 | CD8 Exhausted | TME_to_epi | 0.16596 | -0.069201 | 3 | -0.00045154 |
| SpatialCellChat | Inhousedata_T4469479 | NK | TME_to_epi | 0.15769 | -0.046923 | 4 | 0.011665 |
| SpatialCellChat | Inhousedata_T4469479 | Tfh/CXCL13-like | TME_to_epi | 0.17197 | -0.054807 | 2 | 0.010387 |
| SpatialCellChat | Inhousedata_T4469479 | Treg | TME_to_epi | 0.15674 | -0.05924 | 5 | 0.0035001 |
| SpatialCellChat | Inhousedata_T4469479 | T_NK_unresolved | TME_to_epi | 0.44958 | -0.22463 | 1 | -0.025101 |


## Endothelial全部可辨亚群

| method | section_id | node | direction | strength_low | strength_high | delta_strength | rank_strength_high | delta_share_all_TME | delta_share_within_lineage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | Arterial EC | epi_to_TME | 0.028217 | 0.027833 | -0.00038396 | 2 | -0.0007531 | -0.0029688 |
| COMMOT | GSE171351_A1 | Capillary EC | epi_to_TME | 0.024445 | 0.027166 | 0.0027208 | 4 | 0.003478 | 0.02496 |
| COMMOT | GSE171351_A1 | Lymphatic EC | epi_to_TME | 0.029513 | 0.028089 | -0.0014242 | 1 | -0.002171 | -0.012326 |
| COMMOT | GSE171351_A1 | Venous EC | epi_to_TME | 0.028946 | 0.027818 | -0.0011281 | 3 | -0.0017657 | -0.0096655 |
| COMMOT | GSE171351_A1 | Arterial EC | TME_to_epi | 0.028953 | 0.028941 | -1.1851e-05 | 1 | -0.00082066 | -0.0027962 |
| COMMOT | GSE171351_A1 | Capillary EC | TME_to_epi | 0.024653 | 0.027705 | 0.0030517 | 4 | 0.0033146 | 0.024484 |
| COMMOT | GSE171351_A1 | Lymphatic EC | TME_to_epi | 0.030023 | 0.028782 | -0.0012403 | 2 | -0.0024607 | -0.013675 |
| COMMOT | GSE171351_A1 | Venous EC | TME_to_epi | 0.029148 | 0.028544 | -0.00060435 | 3 | -0.0016028 | -0.008013 |
| COMMOT | Inhousedata_T4469479 | Arterial EC | epi_to_TME | 0.0044398 | 0.0036205 | -0.00081932 | 3 | -0.0016065 | 0.0050005 |
| COMMOT | Inhousedata_T4469479 | Capillary EC | epi_to_TME | 0.0041443 | 0.0030895 | -0.0010548 | 4 | -0.0040575 | -0.014093 |
| COMMOT | Inhousedata_T4469479 | Lymphatic EC | epi_to_TME | 0.0058082 | 0.0049912 | -0.00081699 | 1 | 0.00014626 | 0.023028 |
| COMMOT | Inhousedata_T4469479 | Venous EC | epi_to_TME | 0.0049775 | 0.0037569 | -0.0012206 | 2 | -0.0044657 | -0.013936 |
| COMMOT | Inhousedata_T4469479 | Arterial EC | TME_to_epi | 0.0049126 | 0.0040359 | -0.00087664 | 3 | -0.0019721 | 0.0044847 |
| COMMOT | Inhousedata_T4469479 | Capillary EC | TME_to_epi | 0.0044813 | 0.0033994 | -0.0010819 | 4 | -0.004063 | -0.012549 |
| COMMOT | Inhousedata_T4469479 | Lymphatic EC | TME_to_epi | 0.0060281 | 0.0052775 | -0.0007506 | 1 | 0.00018795 | 0.02467 |
| COMMOT | Inhousedata_T4469479 | Venous EC | TME_to_epi | 0.0056214 | 0.0042496 | -0.0013719 | 2 | -0.0052143 | -0.016606 |
| SpatialCellChat | GSE171351_A1 | Arterial EC | epi_to_TME | 0.23301 | 0.22792 | -0.0050927 | 4 | -0.00060334 | -0.0011305 |
| SpatialCellChat | GSE171351_A1 | Capillary EC | epi_to_TME | 0.45335 | 0.49496 | 0.041615 | 2 | 0.0020061 | 0.017648 |
| SpatialCellChat | GSE171351_A1 | Lymphatic EC | epi_to_TME | 0.29645 | 0.2868 | -0.009652 | 3 | -0.00096342 | -0.0026605 |
| SpatialCellChat | GSE171351_A1 | Venous EC | epi_to_TME | 1.6374 | 1.5863 | -0.051137 | 1 | -0.0051871 | -0.013857 |
| SpatialCellChat | GSE171351_A1 | Arterial EC | TME_to_epi | 0.23584 | 0.23091 | -0.004936 | 4 | -0.00058996 | -0.0018476 |
| SpatialCellChat | GSE171351_A1 | Capillary EC | TME_to_epi | 0.4633 | 0.51796 | 0.054663 | 2 | 0.0027689 | 0.020864 |
| SpatialCellChat | GSE171351_A1 | Lymphatic EC | TME_to_epi | 0.29644 | 0.28555 | -0.010886 | 3 | -0.0010272 | -0.004104 |
| SpatialCellChat | GSE171351_A1 | Venous EC | TME_to_epi | 1.633 | 1.5932 | -0.039746 | 1 | -0.0044247 | -0.014912 |
| SpatialCellChat | Inhousedata_T4469479 | Arterial EC | epi_to_TME | 0.07976 | 0.062929 | -0.016831 | 4 | 5.4718e-05 | 0.0014302 |
| SpatialCellChat | Inhousedata_T4469479 | Capillary EC | epi_to_TME | 0.13292 | 0.097295 | -0.035623 | 2 | -0.0010581 | -0.015502 |
| SpatialCellChat | Inhousedata_T4469479 | Lymphatic EC | epi_to_TME | 0.21072 | 0.17613 | -0.034592 | 1 | 0.0016429 | 0.027096 |
| SpatialCellChat | Inhousedata_T4469479 | Venous EC | epi_to_TME | 0.11859 | 0.087152 | -0.031442 | 3 | -0.00089222 | -0.013025 |
| SpatialCellChat | Inhousedata_T4469479 | Arterial EC | TME_to_epi | 0.083513 | 0.065963 | -0.01755 | 4 | 1.802e-06 | 0.001616 |
| SpatialCellChat | Inhousedata_T4469479 | Capillary EC | TME_to_epi | 0.13879 | 0.10226 | -0.036528 | 2 | -0.0010784 | -0.014233 |
| SpatialCellChat | Inhousedata_T4469479 | Lymphatic EC | TME_to_epi | 0.21092 | 0.17646 | -0.034455 | 1 | 0.0014538 | 0.026755 |
| SpatialCellChat | Inhousedata_T4469479 | Venous EC | TME_to_epi | 0.12371 | 0.090519 | -0.033192 | 3 | -0.001054 | -0.014138 |


SpatialCellChat / GSE171351_A1 / TME_to_epi：High绝对强度最高为Venous EC；最大有符号High−Low为Capillary EC（0.054663）；最大有符号谱系份额变化为Capillary EC（2.086个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / GSE171351_A1 / epi_to_TME：High绝对强度最高为Venous EC；最大有符号High−Low为Capillary EC（0.041615）；最大有符号谱系份额变化为Capillary EC（1.765个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为Lymphatic EC；最大有符号High−Low为Arterial EC（-0.01755）；最大有符号谱系份额变化为Lymphatic EC（2.676个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为Lymphatic EC；最大有符号High−Low为Arterial EC（-0.016831）；最大有符号谱系份额变化为Lymphatic EC（2.710个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / TME_to_epi：High绝对强度最高为Arterial EC；最大有符号High−Low为Capillary EC（0.0030517）；最大有符号谱系份额变化为Capillary EC（2.448个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / epi_to_TME：High绝对强度最高为Lymphatic EC；最大有符号High−Low为Capillary EC（0.0027208）；最大有符号谱系份额变化为Capillary EC（2.496个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为Lymphatic EC；最大有符号High−Low为Lymphatic EC（-0.0007506）；最大有符号谱系份额变化为Lymphatic EC（2.467个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为Lymphatic EC；最大有符号High−Low为Lymphatic EC（-0.00081699）；最大有符号谱系份额变化为Lymphatic EC（2.303个百分点）。若最大值仍为负，则表示下降最少，不是增强。

## Fibroblast与独立Mural

| method | section_id | node | direction | strength_low | strength_high | delta_strength | rank_strength_high | delta_share_all_TME | delta_share_within_lineage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | Matrix_CAF | epi_to_TME | 0.029692 | 0.030043 | 0.00035124 | 1 | 0.0002292 | -0.016879 |
| COMMOT | GSE171351_A1 | iCAF | epi_to_TME | 0.017548 | 0.019196 | 0.0016482 | 4 | 0.0020842 | 0.0043641 |
| COMMOT | GSE171351_A1 | rCAF | epi_to_TME | 0.021069 | 0.022889 | 0.0018204 | 3 | 0.002288 | 0.0036619 |
| COMMOT | GSE171351_A1 | vCAF | epi_to_TME | 0.025676 | 0.028336 | 0.0026596 | 2 | 0.003385 | 0.0088533 |
| COMMOT | GSE171351_A1 | Matrix_CAF | TME_to_epi | 0.030355 | 0.03111 | 0.00075516 | 1 | 0.00014575 | -0.018402 |
| COMMOT | GSE171351_A1 | iCAF | TME_to_epi | 0.0181 | 0.020406 | 0.0023066 | 4 | 0.0025202 | 0.0067399 |
| COMMOT | GSE171351_A1 | rCAF | TME_to_epi | 0.02154 | 0.024072 | 0.0025323 | 3 | 0.0027204 | 0.0059917 |
| COMMOT | GSE171351_A1 | vCAF | TME_to_epi | 0.026296 | 0.029215 | 0.0029192 | 2 | 0.0030952 | 0.0056702 |
| COMMOT | Inhousedata_T4469479 | Matrix_CAF | epi_to_TME | 0.0054159 | 0.0046403 | -0.00077554 | 3 | 1.5256e-05 | -0.00042928 |
| COMMOT | Inhousedata_T4469479 | iCAF | epi_to_TME | 0.0057926 | 0.0062919 | 0.00049934 | 1 | 0.011737 | 0.070395 |
| COMMOT | Inhousedata_T4469479 | rCAF | epi_to_TME | 0.0036387 | 0.002979 | -0.00065967 | 4 | -0.0012124 | -0.0076795 |
| COMMOT | Inhousedata_T4469479 | vCAF | epi_to_TME | 0.0070037 | 0.0048431 | -0.0021606 | 2 | -0.010192 | -0.062286 |
| COMMOT | Inhousedata_T4469479 | Matrix_CAF | TME_to_epi | 0.0061553 | 0.0053249 | -0.00083038 | 3 | -0.00032103 | -0.0016017 |
| COMMOT | Inhousedata_T4469479 | iCAF | TME_to_epi | 0.0059626 | 0.0066534 | 0.00069082 | 1 | 0.011683 | 0.069678 |
| COMMOT | Inhousedata_T4469479 | rCAF | TME_to_epi | 0.0042462 | 0.0034565 | -0.00078975 | 4 | -0.0019614 | -0.011438 |
| COMMOT | Inhousedata_T4469479 | vCAF | TME_to_epi | 0.0077486 | 0.0055567 | -0.0021919 | 2 | -0.0096014 | -0.056638 |
| SpatialCellChat | GSE171351_A1 | Matrix_CAF | epi_to_TME | 0.47218 | 0.47923 | 0.0070464 | 2 | -0.00015079 | -0.0086415 |
| SpatialCellChat | GSE171351_A1 | iCAF | epi_to_TME | 0.20542 | 0.22216 | 0.016735 | 4 | 0.00077809 | 0.00065367 |
| SpatialCellChat | GSE171351_A1 | rCAF | epi_to_TME | 0.21008 | 0.22571 | 0.015626 | 3 | 0.00070386 | 0.00018786 |
| SpatialCellChat | GSE171351_A1 | vCAF | epi_to_TME | 2.0028 | 2.1704 | 0.16759 | 1 | 0.0078591 | 0.0078 |
| SpatialCellChat | GSE171351_A1 | Matrix_CAF | TME_to_epi | 0.48329 | 0.49618 | 0.012897 | 2 | 0.00019545 | -0.0063317 |
| SpatialCellChat | GSE171351_A1 | iCAF | TME_to_epi | 0.2087 | 0.22326 | 0.014569 | 4 | 0.00063364 | 0.00012544 |
| SpatialCellChat | GSE171351_A1 | rCAF | TME_to_epi | 0.21235 | 0.22646 | 0.01411 | 3 | 0.00060117 | -9.9101e-05 |
| SpatialCellChat | GSE171351_A1 | vCAF | TME_to_epi | 2.0426 | 2.2012 | 0.15857 | 1 | 0.007177 | 0.0063054 |
| SpatialCellChat | Inhousedata_T4469479 | Matrix_CAF | epi_to_TME | 0.43612 | 0.36526 | -0.070859 | 2 | 0.0035114 | 0.03595 |
| SpatialCellChat | Inhousedata_T4469479 | iCAF | epi_to_TME | 0.1167 | 0.12084 | 0.0041382 | 3 | 0.0044444 | 0.02833 |
| SpatialCellChat | Inhousedata_T4469479 | rCAF | epi_to_TME | 0.1044 | 0.087522 | -0.01688 | 4 | 0.0008532 | 0.0086734 |
| SpatialCellChat | Inhousedata_T4469479 | vCAF | epi_to_TME | 1.0208 | 0.66096 | -0.3598 | 1 | -0.021209 | -0.072953 |
| SpatialCellChat | Inhousedata_T4469479 | Matrix_CAF | TME_to_epi | 0.46653 | 0.39037 | -0.076162 | 2 | 0.0032232 | 0.033857 |
| SpatialCellChat | Inhousedata_T4469479 | iCAF | TME_to_epi | 0.11752 | 0.1228 | 0.0052861 | 3 | 0.0044059 | 0.027324 |
| SpatialCellChat | Inhousedata_T4469479 | rCAF | TME_to_epi | 0.11036 | 0.091627 | -0.018732 | 4 | 0.00065735 | 0.0074592 |
| SpatialCellChat | Inhousedata_T4469479 | vCAF | TME_to_epi | 1.0596 | 0.69714 | -0.36246 | 1 | -0.020508 | -0.068639 |


SpatialCellChat / GSE171351_A1 / TME_to_epi：High绝对强度最高为vCAF；最大有符号High−Low为vCAF（0.15857）；最大有符号谱系份额变化为vCAF（0.631个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / GSE171351_A1 / epi_to_TME：High绝对强度最高为vCAF；最大有符号High−Low为vCAF（0.16759）；最大有符号谱系份额变化为vCAF（0.780个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为vCAF；最大有符号High−Low为iCAF（0.0052861）；最大有符号谱系份额变化为Matrix_CAF（3.386个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为vCAF；最大有符号High−Low为iCAF（0.0041382）；最大有符号谱系份额变化为Matrix_CAF（3.595个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / TME_to_epi：High绝对强度最高为Matrix_CAF；最大有符号High−Low为vCAF（0.0029192）；最大有符号谱系份额变化为iCAF（0.674个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / epi_to_TME：High绝对强度最高为Matrix_CAF；最大有符号High−Low为vCAF（0.0026596）；最大有符号谱系份额变化为vCAF（0.885个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为iCAF；最大有符号High−Low为iCAF（0.00069082）；最大有符号谱系份额变化为iCAF（6.968个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为iCAF；最大有符号High−Low为iCAF（0.00049934）；最大有符号谱系份额变化为iCAF（7.039个百分点）。若最大值仍为负，则表示下降最少，不是增强。

## Pericyte/SMC原映射

| method | section_id | node | direction | strength_low | strength_high | delta_strength | rank_strength_high | delta_share_all_TME | delta_share_within_lineage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | GSE171351_A1 | Pericyte | epi_to_TME | 0.040533 | 0.034282 | -0.0062509 | 1 | -0.0087914 | -0.049471 |
| COMMOT | GSE171351_A1 | SMC | epi_to_TME | 0.016634 | 0.017695 | 0.0010614 | 2 | 0.001298 | 0.049471 |
| COMMOT | GSE171351_A1 | Pericyte | TME_to_epi | 0.040618 | 0.034638 | -0.00598 | 1 | -0.0089681 | -0.056384 |
| COMMOT | GSE171351_A1 | SMC | TME_to_epi | 0.017598 | 0.019372 | 0.0017739 | 2 | 0.0018358 | 0.056384 |
| COMMOT | Inhousedata_T4469479 | Pericyte | epi_to_TME | 0.0052495 | 0.0053437 | 9.4293e-05 | 1 | 0.0074769 | 0.11186 |
| COMMOT | Inhousedata_T4469479 | SMC | epi_to_TME | 0.0037398 | 0.0023359 | -0.0014039 | 2 | -0.0076489 | -0.11186 |
| COMMOT | Inhousedata_T4469479 | Pericyte | TME_to_epi | 0.005371 | 0.005555 | 0.00018394 | 1 | 0.0070075 | 0.11193 |
| COMMOT | Inhousedata_T4469479 | SMC | TME_to_epi | 0.0048826 | 0.0031828 | -0.0016999 | 2 | -0.0086062 | -0.11193 |
| SpatialCellChat | GSE171351_A1 | Pericyte | epi_to_TME | 0.90591 | 0.79903 | -0.10688 | 2 | -0.0077204 | -0.03807 |
| SpatialCellChat | GSE171351_A1 | SMC | epi_to_TME | 1.8662 | 1.9684 | 0.10221 | 1 | 0.0039935 | 0.03807 |
| SpatialCellChat | GSE171351_A1 | Pericyte | TME_to_epi | 0.9052 | 0.80542 | -0.099775 | 2 | -0.0071974 | -0.033965 |
| SpatialCellChat | GSE171351_A1 | SMC | TME_to_epi | 1.9314 | 2.0191 | 0.087743 | 1 | 0.0029905 | 0.033965 |
| SpatialCellChat | Inhousedata_T4469479 | Pericyte | epi_to_TME | 0.24842 | 0.23918 | -0.0092438 | 1 | 0.0067218 | 0.09871 |
| SpatialCellChat | Inhousedata_T4469479 | SMC | epi_to_TME | 0.13971 | 0.084581 | -0.05513 | 2 | -0.0037957 | -0.09871 |
| SpatialCellChat | Inhousedata_T4469479 | Pericyte | TME_to_epi | 0.25062 | 0.24203 | -0.0085848 | 1 | 0.0064796 | 0.10248 |
| SpatialCellChat | Inhousedata_T4469479 | SMC | TME_to_epi | 0.15367 | 0.093018 | -0.060655 | 2 | -0.0041622 | -0.10248 |


SpatialCellChat / GSE171351_A1 / TME_to_epi：High绝对强度最高为SMC；最大有符号High−Low为SMC（0.087743）；最大有符号谱系份额变化为SMC（3.396个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / GSE171351_A1 / epi_to_TME：High绝对强度最高为SMC；最大有符号High−Low为SMC（0.10221）；最大有符号谱系份额变化为SMC（3.807个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为Pericyte；最大有符号High−Low为Pericyte（-0.0085848）；最大有符号谱系份额变化为Pericyte（10.248个百分点）。若最大值仍为负，则表示下降最少，不是增强。

SpatialCellChat / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为Pericyte；最大有符号High−Low为Pericyte（-0.0092438）；最大有符号谱系份额变化为Pericyte（9.871个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / TME_to_epi：High绝对强度最高为Pericyte；最大有符号High−Low为SMC（0.0017739）；最大有符号谱系份额变化为SMC（5.638个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / GSE171351_A1 / epi_to_TME：High绝对强度最高为Pericyte；最大有符号High−Low为SMC（0.0010614）；最大有符号谱系份额变化为SMC（4.947个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / TME_to_epi：High绝对强度最高为Pericyte；最大有符号High−Low为Pericyte（0.00018394）；最大有符号谱系份额变化为Pericyte（11.193个百分点）。若最大值仍为负，则表示下降最少，不是增强。

COMMOT / Inhousedata_T4469479 / epi_to_TME：High绝对强度最高为Pericyte；最大有符号High−Low为Pericyte（9.4293e-05）；最大有符号谱系份额变化为Pericyte（11.186个百分点）。若最大值仍为负，则表示下降最少，不是增强。

## 跨片重复、算法分歧与scRNA对照

| method | lineage | node | direction | n_high_greater | n_low_greater | n_zero | share_within_lineage_increases | share_within_lineage_decreases |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COMMOT | Endothelial | Arterial EC | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| COMMOT | Endothelial | Arterial EC | epi_to_TME | 0 | 2 | 0 | 1 | 1 |
| COMMOT | Endothelial | Capillary EC | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Endothelial | Capillary EC | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Endothelial | Lymphatic EC | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| COMMOT | Endothelial | Lymphatic EC | epi_to_TME | 0 | 2 | 0 | 1 | 1 |
| COMMOT | Endothelial | Venous EC | TME_to_epi | 0 | 2 | 0 | 0 | 2 |
| COMMOT | Endothelial | Venous EC | epi_to_TME | 0 | 2 | 0 | 0 | 2 |
| COMMOT | Fibroblast | Matrix_CAF | TME_to_epi | 1 | 1 | 0 | 0 | 2 |
| COMMOT | Fibroblast | Matrix_CAF | epi_to_TME | 1 | 1 | 0 | 0 | 2 |
| COMMOT | Fibroblast | iCAF | TME_to_epi | 2 | 0 | 0 | 2 | 0 |
| COMMOT | Fibroblast | iCAF | epi_to_TME | 2 | 0 | 0 | 2 | 0 |
| COMMOT | Fibroblast | rCAF | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Fibroblast | rCAF | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Fibroblast | vCAF | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Fibroblast | vCAF | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Mural | Pericyte | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Mural | Pericyte | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Mural | SMC | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | Mural | SMC | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | T/NK | CD8 Exhausted | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | T/NK | CD8 Exhausted | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | T/NK | NK | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| COMMOT | T/NK | NK | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| COMMOT | T/NK | T_NK_unresolved | TME_to_epi | 1 | 1 | 0 | 0 | 2 |
| COMMOT | T/NK | T_NK_unresolved | epi_to_TME | 1 | 1 | 0 | 0 | 2 |
| COMMOT | T/NK | Tfh/CXCL13-like | TME_to_epi | 1 | 1 | 0 | 2 | 0 |
| COMMOT | T/NK | Tfh/CXCL13-like | epi_to_TME | 1 | 1 | 0 | 2 | 0 |
| COMMOT | T/NK | Treg | TME_to_epi | 0 | 2 | 0 | 0 | 2 |
| COMMOT | T/NK | Treg | epi_to_TME | 0 | 2 | 0 | 0 | 2 |
| SpatialCellChat | Endothelial | Arterial EC | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Arterial EC | epi_to_TME | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Capillary EC | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Capillary EC | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Lymphatic EC | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Lymphatic EC | epi_to_TME | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Endothelial | Venous EC | TME_to_epi | 0 | 2 | 0 | 0 | 2 |
| SpatialCellChat | Endothelial | Venous EC | epi_to_TME | 0 | 2 | 0 | 0 | 2 |
| SpatialCellChat | Fibroblast | Matrix_CAF | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Fibroblast | Matrix_CAF | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Fibroblast | iCAF | TME_to_epi | 2 | 0 | 0 | 2 | 0 |
| SpatialCellChat | Fibroblast | iCAF | epi_to_TME | 2 | 0 | 0 | 2 | 0 |
| SpatialCellChat | Fibroblast | rCAF | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Fibroblast | rCAF | epi_to_TME | 1 | 1 | 0 | 2 | 0 |
| SpatialCellChat | Fibroblast | vCAF | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Fibroblast | vCAF | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Mural | Pericyte | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Mural | Pericyte | epi_to_TME | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | Mural | SMC | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | Mural | SMC | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | CD8 Exhausted | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | CD8 Exhausted | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | NK | TME_to_epi | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | NK | epi_to_TME | 1 | 1 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | T_NK_unresolved | TME_to_epi | 1 | 1 | 0 | 0 | 2 |
| SpatialCellChat | T/NK | T_NK_unresolved | epi_to_TME | 1 | 1 | 0 | 0 | 2 |
| SpatialCellChat | T/NK | Tfh/CXCL13-like | TME_to_epi | 1 | 1 | 0 | 2 | 0 |
| SpatialCellChat | T/NK | Tfh/CXCL13-like | epi_to_TME | 1 | 1 | 0 | 2 | 0 |
| SpatialCellChat | T/NK | Treg | TME_to_epi | 0 | 2 | 0 | 1 | 1 |
| SpatialCellChat | T/NK | Treg | epi_to_TME | 0 | 2 | 0 | 1 | 1 |


只有两片，逐片效应方向比试点合并排名更有解释性。表中次数是切片/记录次数，不是已核验患者重复。相同方向也可能来自Q组成与混合表达共享因素，不能由总强度推出免疫抑制、差预后或因果机制。

GSE171351_A1：重点细群双向比较共30项，COMMOT与SpatialCellChat有30项High−Low方向一致；其余0项保留为算法分歧，详见cross_method_direction_comparison.tsv。

Inhousedata_T4469479：重点细群双向比较共30项，COMMOT与SpatialCellChat有28项High−Low方向一致；其余2项保留为算法分歧，详见cross_method_direction_comparison.tsv。

scRNA重绘结果作为独立尺度对照：ALL Epi→Treg强度和增量均为T/NK第10/11，全TME份额下降而T/NK内部略升；五队列内部份额却仅1增4减。Tip/Cycling EC在scRNA pooled强度突出但可计算队列少；空间参考缺少它们，故不能以空间NE否定scRNA信号。Capillary在scRNA两方向4/4增加；Matrix_CAF→Epi pooled强度和增量领先但三个队列内部份额都下降；Pericyte两方向4/4增加。空间本表需要按自身方向与可测标签判断是否相符，不复制scRNA名次。

## 首片距离敏感性

COMMOT_GSE171351_A1_range_sensitivity：重点双向方向30/30保持一致。

| node | direction | delta_strength_primary | delta_strength_wide | direction_same | delta_share_within_lineage_primary | delta_share_within_lineage_wide |
| --- | --- | --- | --- | --- | --- | --- |
| Treg | epi_to_TME | -0.00024473 | -0.00074247 | True | -0.013464 | -0.014265 |
| Treg | TME_to_epi | -2.7378e-06 | -0.00048273 | True | -0.01423 | -0.014682 |


SpatialCellChat_GSE171351_A1_range_sensitivity：重点双向方向30/30保持一致。

| node | direction | delta_strength_primary | delta_strength_wide | direction_same | delta_share_within_lineage_primary | delta_share_within_lineage_wide |
| --- | --- | --- | --- | --- | --- | --- |
| Treg | epi_to_TME | -0.0038889 | -0.0044563 | True | -0.0084825 | -0.0086661 |
| Treg | TME_to_epi | -0.0045785 | -0.0049387 | True | -0.0082558 | -0.0086778 |


敏感性只有一次较宽范围：COMMOT仅Secreted改到3dNN；SpatialCellChat原生同一个interaction.range同时作用Secreted与ECM，故两者一起扩到3dNN，Contact不变。没有把两次扰动冒称完全相同的参数试验，也未大量扫描阈值。SpatialCellChat的K会随保留链接范围变化，因此mean下降不等同于total_flow下降；两种量分别查表。

## 覆盖、局限与交付

所有Q列达到原验证类型覆盖规则（≥10个spot的mean abundance≥0.01且总abundance>1），逐位置比例不截断。独立组成列存在不等于可确定每个spot中某个LR的细胞来源。原参考只涵盖已建模类型；Q是在其内部归一化。实测spot背景与Q软归属保留，其总强度仅为此模型与采样尺度下的通讯归属。

位置内部筛选、LR数量、spot数量都不当患者H/L检验。未运行额外空间位置置换、tradeSeq、impact、新niche、CNV或预后模型。现成可靠的完整细胞级注释ROI不足，未把HD bin冒充单细胞开展可选ROI。

交付108张High/Low/Both主圈图，全部只保留Epi↔目标边、无群体自环，三图固定位置和共同线宽；另有重点比较、空间分布与敏感性面板。共553项数值/图形检查通过。原推断位置矩阵保留在各方法目录；review.zip只含PDF、报告和重点小表，不含表达、qs/RDS/npz或PNG。

入口：[图索引](FIGURE_INDEX.md)、[方法](methods.md)、[完整逐片表](comparison/epi_state_type_contrasts.tsv)、[跨片方向](comparison/cross_section_summary.tsv)。
