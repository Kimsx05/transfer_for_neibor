# Cervilla MP10–Treg / 内皮空间验证结果

完成日期：2026-11-14T01:59:53 (Asia/Shanghai)。这是单张切片内的空间支持分析，不是患者重复验证。全部授权可执行分支已运行；Treg 相关推断因身份信息不足明确记为不可评估。未运行通讯、配体受体、NicheNet、TF 或恶性判定。

## 完成状态与正式单元

| section | branch | execution_status | evaluability | association_status |
| --- | --- | --- | --- | --- |
| Cervilla | fixed_ROI_bin8_RCTD | completed | partially_evaluable | not_an_association |
| Cervilla | fixed_ROI_bin16_RCTD | completed | partially_evaluable | not_an_association |
| Cervilla | full_existing_cells_fine_reference | completed | partially_evaluable | not_an_association |
| Cervilla | Treg_identity_and_endpoints | evaluated_not_identifiable | not_evaluable | not_tested |
| Cervilla | Endothelial_and_MP10 | completed | partially_evaluable | density_no_clear_support; distance_exploratory_positive |
| Cervilla | niche50_and200 | completed | evaluable_composition; partial_MP10_coverage | K4_no_BHq_below_0.05; K2_200_sensitivity_support |
| D1_Johnson | light_target_registration_check | completed | partially_evaluable | not_tested |

本轮冻结的操作性主线为 **既有 cell 的保守谱系层用于位置与邻域 + 上皮限定的 64 μm 网格内连通小区域用于 MP10**。并非每个细胞都有可靠 MP10，也不意味着 cells 普遍优于 bins。全片直接细分参考预测覆盖93,593个原分析对象、60,670,091原始UMI；原158,144几何对象及QC/counts未改写。8个固定ROI的三路比较如下：

| unit | n | fitted | reliable_fit | Treg_candidates | Treg_multigene | Endo_multigene | Epi_score_eligible | mixed_Epi_fraction | selected |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cell | 2468 | 2466 | NA | 2 | 0 | 13 | 62 | 0.3053 | True |
| bin8 | 7642 | 5263 | 3935 | 1 | 0 | 23 | 6 | 0.08216 | False |
| bin16 | 1910 | 1658 | 1302 | 0 | 0 | 19 | 86 | 0.2996 | False |

8 μm 保留原始定位优势，但相同ROI仅6个Epi高贡献bin通过个体MP10测量资格；16 μm虽提高检测量，按原生selected-model规则，Epi高贡献bin中约30%仍有Epi/胶原共检出，不能用模型Epi比例恢复纯上皮RNA。cells约30.5%同类共检出也不被当作纯上皮评分。**修正后的比较不支持cells比16 μm混合更少，bin/cell没有明确通用赢家。** 保留已冻结cells操作性主线，是为核支持位置与Epi-only连通区域求和，不能宣称身份或评分准确率更高；主结果因此保持探索性。8 μm原始信号应配合使用，不按Treg数量、分配UMI或显著性重选。16 μm从全部297,675 supplied8 bins聚合成75,293格，64,528,669 UMI逐基因及全局守恒；74,730格通过本轮宽松表达QC。没有先删除低UMI8 μm格。

## 目标可信程度

- 冻结30类精细参考训练实际18,289个细胞；按 `orig.ident` 样本留出，不宣称患者独立。TIGIT不参与分类/RCTD判别。Treg内部参考精度在预计100/300/1000 UMI约0.61/0.80/0.89，且CD4/Tfh仍有混淆；这不是HD准确率。
- 上皮参考1800个细胞中，KRT5/KRT14/KRT6A/KRT6B/KRT17/COL17A1分别在352/165/239/64/1144/161个细胞检出，说明参考存在这些表达模式，但不是独立病理亚型验证。当前SingleR实际DE特征包含KRT7/KRT19，而上述基底/鳞状面板、EPCAM/UPK2未入选；它们仍进入原始RNA复核面板，不能把“基因存在”等同“分类器实际使用”。这限制对基底/鳞状上皮覆盖的强断言。Treg的FOXP3/IL2RA/CTLA4/IKZF2实际进入判别，LRRC32未入选，TIGIT明确排除。
- Cervilla有144个原始Treg预测候选，其中28个具有至少1个原始T-panel基因，9个有FOXP3或IL2RA；没有足够的成组Treg/T-lineage联合证据。候选组平均/检出图也不支持把全部候选当作Treg。8/16 μm RCTD、既有4 μm保守归属ROI均未提供可信补救。**Treg、Treg/CD4、Treg/T、Treg距离及反向Treg/Tconv/Tfh比较全部NA**，不是零Treg，也不是不存在联系。
- 204个Endothelial支持候选有参考和至少2个原始内皮marker；PECAM1/VWF等群体证据明显。108个同时有至少2个Epi或Fib marker，因此不解释为完整、纯净血管细胞图谱。固定排除这类共检出后剩96个位置用于一次敏感性分析。内皮亚型未作为正式终点。
- 38,027个支持上皮对象作为niche中心，30,009个未标记混合的上皮对象可参与Epi-only区域构建；仅1,394个单对象达到预定MP10支持条件。8,008个候选连通区域中1,209个可评分。所有20基因均在完整18,085基因输入中，未使用HVG/BANKSY表达或TME表达评分。
- 主MP10严格复用历史**平均并列秩**公式，K=905；与原SciPy实现1800个检查实例最大误差0。原生AUCell随机并列/排除零秩仅为敏感性。合格区域半深度排序ρ=0.716，预计3000UMI排序ρ=0.483；检测量调整后仍只能做区域层面解释。部分合格区域AUC=0表示已观察基因未进入高秩窗口，并不等于未测量；不合格区域不分Low。

## 50 / 200 μm 空间关系

| radius_um | epi_anchors | neighborhood_eligible | unresolved | training | predicted | score_regions | association_regions | distance_observable_regions | Endo_positive_regions | mean_known_fraction | median_known_neighbors | Endo_objects |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | 38027 | 37897 | 130 | 15102 | 22795 | 1209 | 1209 | 916 | 187 | 0.8124 | 39 | 204 |
| 200 | 38027 | 38027 | 0 | 15221 | 22806 | 1209 | 1202 | 916 | 989 | 0.7984 | 581 | 204 |

中心自身被排除，邻居池是全部原QC cells，包含周围Epi和全部可支持TME，Unknown从组成分子和条件分母同时排除。50/200为累计真实物理圆，不是环带、kNN或图跳数。类型贡献是核支持对象计数，不是真实无漏检细胞数。区域MP10没有复制成独立细胞观测。

主要条件结果（效应为标准化秩斜率；400 μm空间区块，499次）：

| radius_um | predictor | n | n_blocks | effect | low | high | P | BH_q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | Endothelial_density | 1209 | 85 | -0.0284 | -0.09867 | 0.03796 | 0.382 | 0.4366 |
| 50 | Endothelial_distance | 916 | 68 | 0.1269 | 0.04017 | 0.2318 | 0.008 | 0.016 |
| 200 | Endothelial_density | 1202 | 85 | -0.07134 | -0.172 | 0.01847 | 0.132 | 0.2112 |
| 200 | Endothelial_distance | 916 | 68 | 0.1434 | 0.0458 | 0.2544 | 0.008 | 0.016 |

**内皮支持密度未检出明确关联；距离有小幅正向关系，即较高MP10偏向离已识别内皮位置更远，而非内皮富集。** 最近距离本身不依赖邻域半径，两尺度共享916个可观察区域，不能当成两次独立复制。更大800 μm区块的距离检验P=0.048/0.040、敏感性家族q=0.192；排除混合Endo后距离q约0.099。加入已有BANKSY区域、边缘及可注释支持的条件分析方向仍相同，但这是同一片内的相关证据。因此总体结论为**探索性距离信号，未达到跨质量设定均稳定的强支持**。

| radius_um | predictor | n_target_sensitivity | n | effect | low | high | P | BH_q |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | Endothelial_density_unmixed | 96 | 1209 | -0.07386 | -0.1563 | -0.003288 | 0.074 | 0.09867 |
| 50 | Endothelial_distance_unmixed | 96 | 773 | 0.1521 | 0.02092 | 0.3145 | 0.05 | 0.09867 |
| 200 | Endothelial_density_unmixed | 96 | 1202 | -0.08893 | -0.2048 | 0.02309 | 0.15 | 0.15 |
| 200 | Endothelial_distance_unmixed | 96 | 773 | 0.1616 | 0.01997 | 0.3258 | 0.048 | 0.09867 |

距离超过最近观察边界时删失；缺测、未拟合与零值分开。现有组织区域条件分析只将区域质心附近<=8 μm的原bin域作为空间背景协变量，不转移细胞身份。未据其结果改写niche或MP10定义。

## niche 组成与 MP10

两个尺度均按冻结规则选K=4，并保留K=2。50 μm K4的种子最低ARI0.986、区块ARI0.972；200 μm分别1.000、0.952。空间连续性不是准确率。组成描述：

| radius_um | niche | n | composition_description | mean_known_fraction |
| --- | --- | --- | --- | --- |
| 50 | 1 | 14676 | Epithelial 88% / Fibroblast 10% / Mural 1% | 0.8351 |
| 50 | 2 | 12426 | Epithelial 98% / Myeloid 1% / Fibroblast 1% | 0.826 |
| 50 | 3 | 2045 | Fibroblast 65% / Epithelial 28% / Mural 3% | 0.6706 |
| 50 | 4 | 8750 | Epithelial 64% / Fibroblast 31% / Mural 3% | 0.7966 |
| 200 | 1 | 2857 | Epithelial 47% / Fibroblast 34% / Plasma 7% | 0.7943 |
| 200 | 2 | 18904 | Epithelial 90% / Fibroblast 8% / Mural 1% | 0.8228 |
| 200 | 3 | 12721 | Epithelial 67% / Fibroblast 28% / Mural 3% | 0.7977 |
| 200 | 4 | 3545 | Fibroblast 57% / Epithelial 38% / Mural 3% | 0.6743 |

这些niche主要解释上皮优势、上皮–间质梯度及较大尺度少量髓系/浆细胞成分；没有形成可验证的Treg niche。图像和原始Epi/胶原信号支持广义组织差异，但共享marker不是独立真值。两尺度编号不对应，交叉表在 `niche_scale_crosstab.tsv`。

K4的MP10主要家族没有q<0.05（最小q=0.056）。50 μm niche1较高、niche4较低和200 μm上皮优势niche2较高仅为方向性趋势。50 μm间质优势niche3只有1个可评分区域，正式P撤为NA；不能把该niche判为MP10低。K2的200 μm粗结构敏感性出现上皮优势区域较高MP10(q=0.016)，但不据此替换冻结K4主结果。H50/H10/H05与L50的数量、比例及描述性优势比完整保留，不选择最显著分组。

**解释范围**：稳定组成差异成立；细粒度niche–MP10关联未检出明确多重校正支持。评分资格主要偏向含足够可信上皮的区域，不能推广到全部38,027中心或上皮稀少niche。

## D1 轻量检查

原364,805 QC8 μm bins及裁剪原点保留。原始T-lineage多基因支持14,443 bins，内皮多基因支持2,993，Treg-panel与T-lineage联合原始线索232，均不等于已确认亚型。既有 `conflict_focus` 和一致对照ROI的H&E/原始Epi/T/Endo图已实际查看；可支持粗组织配准，不能凭图定Treg。3000个固定抽样Epi候选全部20个MP10基因在输入内，单bin检出成员中位数1、126个检出>=5，提示仍有单bin评分稀疏问题。没有对D1做新全片精细分类、RCTD、niche或关联。D1评估结论为**部分可评估，值得后续目标验证但本轮不扩展**。

## 空间复核与后续适用性

已查看固定T-lineage、vascular、dense、far-UMI及补充原始panel图。`vascular_focus` 的命名仅来自预先选中的内皮panel热点，不是病理确认的血管ROI；H&E显示核支持结构，但现有扩展边界仍可能跨越邻近细胞。D1冲突ROI有组织化Epi表达，T/Endo信号较稀疏。无输入覆盖范围不作为目标阴性区域。

可进入下一轮的是：可追溯的上皮限定区域评分、广义上皮–间质组成、保守内皮支持位置及其不确定性。**当前不具备Treg特异通讯推断条件**；内皮相关结果也需独立细胞/血管身份或蛋白验证、更多切片及深度稳健性后再解释机制。未检出明确支持不是证明无生物学联系，本轮不输出通讯增强结论。

## 完成、验证与限制

Python 53项针对性核验、R原对象实际重读均通过：ID、整数counts、坐标、原QC、16 μm逐基因守恒、Epi-only区域求和、邻域分母/中心排除及两尺度独立字段。重读对象生成 `R_loader_readback_niche50.png`。原counts及正式上游文件只读，没有保存新全量Seurat/H5AD。

过程修复透明记录：RCTD工作进程命名空间问题以进程内FORK适配解决；原生输出的位置权重列按实际第一/第二类型解码；最终又核对get_beta_doublet，对singlet采用第一类=1，只有certain doublet保留两类权重，未定贡献为NA。该修正撤回“cells混合低于16 μm”的初始理由；修正前汇总保留，原冻结cell主线/MP10/niche和关联未变，单元选择置信度明确为provisional；历史MP10平均秩在关联前纠正；全片邻域训练索引修复为位置索引。最终统计检查发现极小niche组不能借总样本量产生P，增加每比较侧>=20区域且>=8区块要求，50 μm niche3撤去P。后续只增加任务要求的组织背景和一次内皮混合敏感性，未改动主单元、评分、中心、K或主要效应。

没有未启动却标为完成的计算分支。未进行有效检验的目标是Treg及其相对比例/距离/反向检验；D1正式关联依范围未做；所有通讯分析依停止条件未做。主要剩余限制为参考跨平台适用性、Treg稀疏、细胞归属混合、内皮不完全检测、Epi评分覆盖偏倚、单片空间自相关及边界删失。

新增目录占用(含review)：0.481 GiB；review.zip：26.4 MiB。实际资源记录见RESOURCE_USAGE.json和logs；默认4 CPU，未升级共享环境。
