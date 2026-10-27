# MP10 区域与界面分析方法（v1）

输入为 spatial_ecology_top20_v2 冻结的 63 张 ordinary Visium、7 个 dataset、265,002 个 QC 保留组织位置。按 dataset_id、section_id、完整 position_id 合并，保留原始坐标、原分数及 CNV 标签。全部 30 类 c2l posterior mean 来自已完成的 v2_30type；不重新推断、不复制表达母对象。空间 Epithelial 是合并标签，不能由 scRNA 标签列表还原细分。Pericyte/SMC 独立属于 Mural。

## 证据与区域

旧 ROI 为 cnv_candidate_v1，8,257 个候选均无已填写 reviewed_label。原 region_auto、region_final、region_review_status、roi_id 保留。没有可对齐且已核实病理轮廓。63 张已有 CopyKAT 逐 spot 结果；41 张自动基线未经外部验证，19 张混合模型未收敛、3 张另有不确定标记。原生 aneuploid/diploid 和 MSQ 保留为诊断性证据；MSQ 定义为 native final profile 平方均值，不跨片合并、不移植 Yan 阈值。scRNA inferCNV 没有可靠空间映射，不能当 spot CNV。本轮 CNV-supported 状态为 Unknown，并非 CNV 阴性。

主区域沿用冻结 Epithelial / 全部30类 mean abundance ≥0.25 门槛，命名 Epi-rich region。区域参数在读取目标差异之前写入 frozen_rules.json。仅使用真实六邻接且实际间距≤1.35×每片中位最近邻距 d_NN；连通分量逐个编号，保留单点、孔洞与分离岛；不平滑、不填孔、不建跨间质凸包。Epi 丰度高是定义的一部分，不能作为独立验证。

几何在原图像像素框架中以真实格点拟合 affine，映射标准六角 Voronoi 单元，保留镜像/旋转；原位置 x/y 不改写。逐片记录残差，最大残差容限0.35 d_NN（冻结技术容限，非临床精度）。区域/支持 mask 是观测网格的派生覆盖，不是组织学分割或实际 spot 捕获圆面积。面积只报 grid image-pixel²，不报 μm²。输入坐标沿用已审计 lowres image reference，不再次套 scale factor。图像上视觉检查只排查明显错位。

## 边界与距离

两侧均有有效观测且属于不同 Epi mask 的共享单元边是有效 Epi-rich interface。完整原 tissue 网格用来区分 QC 缺失观测边；其余组织/捕获边保留原因未决。另按几何区分支持集外包络与内部孔洞。缺失组织边和可能截断边不参加区域界面距离。每个线段有稳定 boundary_id；轮廓同时保存 GeoJSON（本地图像坐标，非经纬度）。不从 H&E 自动猜测组织学身份。

计算点到全部候选线段的精确欧氏最短距离，不以最近线段中点近似。区域内只匹配自己的 region；外侧匹配同一有效组织连通分量中最近的 region，故每个位置只计一次。内部负、外侧正；distance_to_boundary 为绝对值，signed_distance 才携带符号。接口表 region_id 为所参照界面的区域，position_region_id 为位置本身区域，外侧为 Outside。无同分量有效边界则 NA，不填0。

w=d_NN（每片图像像素）。core d<-w；inner edge -w≤d≤0；outer neighborhood 0<d≤w。距任一缺失/组织/捕获支持边≤w 者保留距离但主统计排除；far tissue 仅在其支持边距离>w时报告。小区域无core即不可评估。代表片每dataset按 section_id 字典序首先选一张，不根据目标富集选择；主规则应用全片后，代表片另做2w敏感性。

## 分数、组成与统计

主 MP10 是冻结 mixed-spot raw Top20 AUCell，不称单个上皮细胞分数。可用的既有 Epi-component matched-model-universe 分数另列，严格沿用 epi_score_valid；不能与原 raw score 的数轴混合。原分数与 model_definition_id 保留。冻结输入没有本任务可直接复用的 raw-score high/low，本轮仅在冻结 Epi anchors 内按片中位数派生 relative_high（>中位数）、relative_low（<）、tied_unclassified（=）；非 anchor 不分组。high比例分母是该region/band所有Epi anchors，包含未分类相等者，并报告分母。所有原 scRNA 全Epi准入不变。

先汇总全部9 major，再列完整30个实际可识别类型（重点17种T/NK、EC、Fibroblast/Mural）。mean_model_abundance 是模型期望值，不是实测细胞数；mean_local_fraction 是每位置全30类比例的平均；pooled_modeled_fraction 是组内该类估计量和/全类型估计量和；mean_lineage_fraction 是该类/所属major在位置内的比例再平均。Treg/全部类型与Treg/T_NK分别报告。major 的谱系内比例恒为1，仅是分母定义。

先每section报告配对 inner_edge−inner_core、outer_neighborhood−inner_edge。患者仅来自修订来源身份表，同患者多片效应等权平均；匿名/假定独立记录单列敏感性，GSE285715未解配对不当6名独立患者。跨dataset只分层描述方向、分位数、覆盖，无spot级P值或机制推断。线性MP10~Epi丰度残差仅为片内组成敏感性，不能校正全部混合或证明因果。效应表保留左右位置数；连续曲线展示仅使用每section/bin≥10位置，原表保留全部覆盖，无平滑外推。

local niche、candidate meta-niche、已有R1/R2 accepted labels原样复用。meta-niche未确认，始终标candidate。R1/R2仅适用于原共同Epi-anchor宇宙；未被赋值不当作生物学缺失。crosswalk按同位置交叉，不强制编号一一对应。niche_target_summary为这些标签下当前中心位置组成，非重新计算ring组成。无新niche/CM聚类、无通讯推断。

## 复现

在本目录之外不写入。现有 /usr/bin/python3.12 + spatial_ecology_top20_v2/runtime 提供运行依赖；数值线程1，无GPU。

```bash
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 'mp10_cnv_regions_interface_v1/scripts/run_analysis.py'
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 'mp10_cnv_regions_interface_v1/scripts/summarize_and_validate.py'
```

以上命令在项目根执行；重跑只更新本任务独占输出。pilot用 --pilot，完成视觉校验后使用完整命令。正式图为透明矢量PDF（H&E保留原始位图）+300dpi透明PNG，少量白底预览在logs。review.zip只含核心图、小表和报告。


最终独立复核：运行 scripts/independent_checks.py，用独立Shapely点到线距离抽查62张有界面切片共744个位置，核查符号与mask内外；全部265,002位置原x/y和raw MP10逐值相符。最后运行 scripts/finalize_review.py，补充明确缺失说明、niche重叠和1w/2w方向表并封装review。
