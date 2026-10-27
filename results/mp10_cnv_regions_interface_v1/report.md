# MP10 CNV 区域与界面分析报告



完成63张切片、7个dataset、265,002个质控保留位置；生成1,916个Epi-rich连通区域、57,146条有效界面线段。256,788个位置可给有符号距离；排除支持边影响后217,538个位置进入主分带描述。输出为组成支持的区域，尚非病理确认的肿瘤边界。



## 阅读重点

- 主比较为同片“外侧邻域−内侧边缘”。混合spot MP10的dataset中位差在6/7个dataset为负，GSE171351为正；这是混合组成层面的描述，不能据此认定上皮细胞内状态降低。
- Inhouse的7/7张可比较切片均显示外侧混合MP10较低、Treg在全30类中的占比较高；但Treg在T/NK内的比例仅3/7为正。跨dataset也有类似分母依赖，不能概括成统一的Treg特异富集。
- 已有Epi-component分数仅覆盖Epi anchors，外侧邻域没有对应分数，故跨界面Epi-component比较全部不可评估；core/inner-edge内可用比较另表保留。空白不是零效应。
- Fibroblast全类型占比的外侧−内侧dataset中位差在7/7为正；EC为6/7，GSE238145为负。区域本身由Epi比例定义，其他谱系比例上升可能含组成闭合效应；模型丰度与谱系内分母表同时提供。

## 证据覆盖

| dataset_id | sections | positions | epi_regions | signed_positions | eligible_positions |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | 13 | 56357 | 319 | 55689 | 48203 |
| GSE171351 | 4 | 4086 | 16 | 4066 | 3293 |
| GSE238145 | 10 | 74264 | 465 | 69192 | 61512 |
| GSE285715 | 6 | 17241 | 53 | 16889 | 13510 |
| GSE319536 | 16 | 58041 | 739 | 56098 | 41966 |
| Inhousedata__UNVERIFIED | 7 | 31213 | 189 | 31151 | 28134 |
| NPJ__UNVERIFIED | 7 | 23800 | 135 | 23703 | 20920 |



身份核验：30张切片对应29个研究内来源明确患者；其余记录按既有修订政策分别列为匿名/假定样本或配对未解。高分辨率补充输入逐项见 additional_platform_coverage.tsv；已有HD表达/图像并不提供本轮可信上皮/目标注释，Xenium沿用architecture-only限制。



## 区域与边界可靠性

旧ROI均处于pending review，原标签完整保留。主区域使用冻结Epi≥0.25及真实六邻接，保留小岛、孔洞。原生CNV仅作伴随诊断：41片自动参考未经外部验证、19片混合基线未收敛、3片其他不确定。因此未制造CNV高低阈值或病理标签。逐区域CNV profile覆盖、原生aneuploid比例、MSQ可查cnv_region_evidence.tsv。

坐标沿用原图像像素；最大格点拟合误差0.268 d_NN，几何界面有这一定位误差。微米距离/物理面积不予声称。QC缺失、外包络与内部孔洞边独立标注，组织外缘与截断原因不够明确时保持unresolved。



## MP10 位置与混合影响

| dataset_id | metric | n_units_evaluable | n_positive | n_negative | median_effect |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| CanDis__UNVERIFIED | MP10_linear_epi_residual | 13 | 8 | 5 | 0.00224 |
| CanDis__UNVERIFIED | MP10_mixed_AUCell | 13 | 1 | 12 | -0.00621 |
| GSE171351 | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| GSE171351 | MP10_linear_epi_residual | 3 | 2 | 1 | 0.00316 |
| GSE171351 | MP10_mixed_AUCell | 3 | 2 | 1 | 0.00504 |
| GSE238145 | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| GSE238145 | MP10_linear_epi_residual | 9 | 5 | 4 | 0.00054 |
| GSE238145 | MP10_mixed_AUCell | 9 | 2 | 7 | -0.00539 |
| GSE285715 | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| GSE285715 | MP10_linear_epi_residual | 6 | 4 | 2 | 0.00582 |
| GSE285715 | MP10_mixed_AUCell | 6 | 1 | 5 | -0.0068 |
| GSE319536 | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| GSE319536 | MP10_linear_epi_residual | 16 | 11 | 5 | 0.00478 |
| GSE319536 | MP10_mixed_AUCell | 16 | 5 | 11 | -0.00237 |
| Inhousedata__UNVERIFIED | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| Inhousedata__UNVERIFIED | MP10_linear_epi_residual | 7 | 4 | 3 | 0.00031 |
| Inhousedata__UNVERIFIED | MP10_mixed_AUCell | 7 | 0 | 7 | -0.00779 |
| NPJ__UNVERIFIED | MP10_epi_component_existing | 0 | 0 | 0 | NA |
| NPJ__UNVERIFIED | MP10_linear_epi_residual | 7 | 6 | 1 | 0.00822 |
| NPJ__UNVERIFIED | MP10_mixed_AUCell | 7 | 0 | 7 | -0.00515 |

表为外侧邻域−内侧边缘；连续分数优先。raw MP10是位置混合表达；Epi成分分数仅使用既有有效估计，单独比较。线性残差是组成敏感性，不等同细胞内调控证据。core对edge结果、高低anchor频率及分母见配对效应和region_mp10_summary表。Epi随分带变化包含区域定义本身。



## T/NK 与 Treg

| dataset_id | metric | n_units_evaluable | n_positive | n_negative | median_effect |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | fine:Treg:T_NK_fraction | 13 | 8 | 5 | 0.00723 |
| CanDis__UNVERIFIED | fine:Treg:all_modeled_fraction | 13 | 12 | 1 | 0.01081 |
| CanDis__UNVERIFIED | major:T_NK:all_modeled_fraction | 13 | 13 | 0 | 0.04792 |
| GSE171351 | fine:Treg:T_NK_fraction | 3 | 1 | 2 | -0.02671 |
| GSE171351 | fine:Treg:all_modeled_fraction | 3 | 1 | 2 | -0.00156 |
| GSE171351 | major:T_NK:all_modeled_fraction | 3 | 3 | 0 | 0.00651 |
| GSE238145 | fine:Treg:T_NK_fraction | 9 | 3 | 6 | -0.00718 |
| GSE238145 | fine:Treg:all_modeled_fraction | 9 | 8 | 1 | 0.00831 |
| GSE238145 | major:T_NK:all_modeled_fraction | 9 | 9 | 0 | 0.0641 |
| GSE285715 | fine:Treg:T_NK_fraction | 6 | 1 | 5 | -0.00441 |
| GSE285715 | fine:Treg:all_modeled_fraction | 6 | 5 | 1 | 0.01005 |
| GSE285715 | major:T_NK:all_modeled_fraction | 6 | 6 | 0 | 0.07457 |
| GSE319536 | fine:Treg:T_NK_fraction | 16 | 10 | 6 | 0.00627 |
| GSE319536 | fine:Treg:all_modeled_fraction | 16 | 14 | 2 | 0.01038 |
| GSE319536 | major:T_NK:all_modeled_fraction | 16 | 15 | 1 | 0.07779 |
| Inhousedata__UNVERIFIED | fine:Treg:T_NK_fraction | 7 | 3 | 4 | -0.00397 |
| Inhousedata__UNVERIFIED | fine:Treg:all_modeled_fraction | 7 | 7 | 0 | 0.01709 |
| Inhousedata__UNVERIFIED | major:T_NK:all_modeled_fraction | 7 | 6 | 1 | 0.11843 |
| NPJ__UNVERIFIED | fine:Treg:T_NK_fraction | 7 | 2 | 5 | -0.01316 |
| NPJ__UNVERIFIED | fine:Treg:all_modeled_fraction | 7 | 7 | 0 | 0.00887 |
| NPJ__UNVERIFIED | major:T_NK:all_modeled_fraction | 7 | 7 | 0 | 0.06995 |

Treg两种分母在121个可评估section×contrast中有47个方向相反；这说明全组织组成和谱系内组成回答不同问题。它们是c2l模型估计，不是直接计数。



## Endothelial

| dataset_id | metric | n_units_evaluable | n_positive | n_negative | median_effect |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | major:Endothelial:all_modeled_fraction | 13 | 11 | 2 | 0.01547 |
| GSE171351 | major:Endothelial:all_modeled_fraction | 3 | 2 | 1 | 0.00292 |
| GSE238145 | major:Endothelial:all_modeled_fraction | 9 | 3 | 6 | -0.00542 |
| GSE285715 | major:Endothelial:all_modeled_fraction | 6 | 6 | 0 | 0.06187 |
| GSE319536 | major:Endothelial:all_modeled_fraction | 16 | 14 | 2 | 0.01996 |
| Inhousedata__UNVERIFIED | major:Endothelial:all_modeled_fraction | 7 | 7 | 0 | 0.024 |
| NPJ__UNVERIFIED | major:Endothelial:all_modeled_fraction | 7 | 7 | 0 | 0.0257 |



## Fibroblast 与 Mural

| dataset_id | metric | n_units_evaluable | n_positive | n_negative | median_effect |
| --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | major:Fibroblast:all_modeled_fraction | 13 | 12 | 1 | 0.03729 |
| CanDis__UNVERIFIED | major:Mural:all_modeled_fraction | 13 | 13 | 0 | 0.02016 |
| GSE171351 | major:Fibroblast:all_modeled_fraction | 3 | 2 | 1 | 0.0064 |
| GSE171351 | major:Mural:all_modeled_fraction | 3 | 2 | 1 | 0.01254 |
| GSE238145 | major:Fibroblast:all_modeled_fraction | 9 | 7 | 2 | 0.0166 |
| GSE238145 | major:Mural:all_modeled_fraction | 9 | 8 | 1 | 0.01866 |
| GSE285715 | major:Fibroblast:all_modeled_fraction | 6 | 6 | 0 | 0.09486 |
| GSE285715 | major:Mural:all_modeled_fraction | 6 | 6 | 0 | 0.04377 |
| GSE319536 | major:Fibroblast:all_modeled_fraction | 16 | 15 | 1 | 0.03638 |
| GSE319536 | major:Mural:all_modeled_fraction | 16 | 16 | 0 | 0.01251 |
| Inhousedata__UNVERIFIED | major:Fibroblast:all_modeled_fraction | 7 | 7 | 0 | 0.0234 |
| Inhousedata__UNVERIFIED | major:Mural:all_modeled_fraction | 7 | 7 | 0 | 0.01787 |
| NPJ__UNVERIFIED | major:Fibroblast:all_modeled_fraction | 7 | 7 | 0 | 0.0636 |
| NPJ__UNVERIFIED | major:Mural:all_modeled_fraction | 7 | 7 | 0 | 0.05556 |



## Niche 对照

同位置local与R1/R2标签的归一化互信息范围为0.050–0.473，两者不能直接互换；这是已赋值共同anchors上的描述性重叠，不是独立样本验证。详见niche_system_overlap_summary.tsv。

local niche、candidate meta-niche与既有R1/R2共识标签全部复用。niche_region_crosswalk保存区域/带交叉；local_R12_overlap保存同位置重叠，避免把不同算法编号当作等价。R1/R2只覆盖其原Epi-anchor宇宙，因而在外侧邻域通常不可评估。candidate meta-niche尚未冻结，不用于确认性结论。代表图并列三类标签、区域界面、MP10和实际细分谱系。



## 不一致与限制

GSE238145_B13的4,881个位置均未达到冻结Epi门槛，因此没有本轮Epi-rich区域或界面；该片仍保留完整输入、CNV、MP10和niche记录。Inhousedata_T4384654和GSE171351_D1无可评估core；不为补齐分带改阈值。

各dataset的方向、分布和可评估数量分别报告。患者层面结果只汇总可核实身份，匿名记录另表；不把spot数当独立生物重复，不计算探索性小P值。局部CNV与Epi不一致、阈值附近网格误差、孔洞和截断可影响边界；本轮只完成代表片2w敏感性，不按目标效应选带宽。单点岛和部分切片无core/无外侧支持时明确NA。

当前现象是可继续追踪的空间组成线索，不能由共定位推导通讯、免疫抑制、机制或预后。没有新增全量CNV、CellChat/COMMOT、空间注释或niche聚类，未改写scRNA全Epi定义。



代表片1w与2w敏感性在48个可比较的section×contrast×指标组合中，5个方向改变；结果见representative_bandwidth_direction.tsv。带宽结果未用于选择主设置。

## 读取接口

- regions/position_regions.tsv.gz：每位置原始候选/最终mask、区域与旧ROI标签。

- interfaces/position_interface_distances.tsv.gz：每位置最近有效界面的region_id/boundary_id、绝对和有符号距离、单位、分带及评价状态。主分析筛evaluation_status=evaluable。

- regions/masks/*.geojson.gz：逐片区域、有效观测支持及分类边界，图像像素坐标。

- input_review/position_features.tsv.gz：保留的MP10/CNV/组成门槛/niche字段，通过dataset_id+section_id+position_id连接。

- target_profiles/：区域/分带组成、连续曲线、section与verified patient效应、dataset分层描述。

- figure_index.tsv与review.zip提供轻量审阅；methods.md说明分母、几何、身份和复现。
