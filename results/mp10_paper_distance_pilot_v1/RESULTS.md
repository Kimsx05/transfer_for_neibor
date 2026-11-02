# MP10 距离先导结果

7 个 dataset 各取 1 张原 v8 代表切片，共 31,024 个组织 spot，15,055 个冻结 Epi queries。7 片均为肿瘤研究标本，未纳入邻近正常对照。所有片都有正式 posterior mean、有效原始 MP10、图像与表达。R12 原标签覆盖 13,248 个 Epi anchors (88.0%)，全组织 dataset-local 标签覆盖 31,024 个 spot。匿名病例与未解析患者别名保留原身份。审阅可从 [先导总览](pilot_overview.pdf) 开始，其余逐片原图见 [图索引](figure_index.tsv)。

## 三条路线的距离梯度

下表为 P80 全组织 target、含同 spot 零距离的片内 Pearson r。正值表示 MP10 随着到目标的距离增加而升高，负值表示靠近目标的 Epi queries 中 MP10 较高。所有效应保留；每片 Pearson、Spearman、零距离数、可用距离范围和四种条件见 [距离表](tables/associations.tsv)。

| 切片 | CAF | TAM | Treg | 空间观察 |
|---|---:|---:|---:|---|
| CanDis Sample9 | +0.141 | +0.180 | +0.161 | 三条均为弱正梯度，去掉零距离后分别为 +0.096、+0.136、+0.141。 |
| GSE171351 B1 | +0.038 | +0.026 | +0.035 | 三条接近零，目标重叠并没有带来清楚的 MP10 梯度。 |
| GSE238145 B14 | −0.098 | −0.034 | +0.085 | 均弱；CAF 有较长的非重叠距离，正距离 r 仍为 −0.098。 |
| GSE285715 T2 | +0.128 | +0.232 | −0.110 | TAM 与 Treg 方向相反，去掉零距离后分别为 +0.158 和 −0.144。 |
| GSE319536 BC01 | +0.132 | +0.014 | +0.014 | 主要是弱 CAF 梯度。只取边界带 CAF 后为 −0.069。 |
| Inhousedata T4469479 | +0.192 | −0.151 | +0.152 | TAM 的负梯度保留于仅正距离数据；改用 P90 后反转为 +0.205。 |
| NPJ sample6 | +0.193 | +0.361 | +0.409 | TAM 与 Treg 最清楚。仅正距离分别为 +0.312、+0.367；P90 为 +0.338、+0.430。 |

**NPJ sample6 的现象主要见于候选内侧的 Epi queries。** 内侧 TAM、Treg 的 r 均约 +0.42，外侧分别为 −0.04、−0.02。原 R12 标签 1 (旧图 niche2) 内，两者分别为 +0.39、+0.47；标签 0 的对应值为 −0.01、+0.06。原始参考 R1 profile 中，标签 1 的上皮和 GC B 平均相对组成为 59.1% 和 15.0%，用于说明该邻域的背景。进一步限于 v8 E0 时，NPJ 的 TAM、Treg 主梯度仍为 +0.372、+0.454。全组织 TAM、Treg targets 中 68% 和 74% 位于候选外侧。只保留边界带 target 后，700 个 TAM targets 剩 38 个，700 个 Treg targets 剩 102 个，r 分别变为 −0.077、−0.063。这张片支持与全组织目标分布相关的梯度，不能据此定位为 200 μm 边界带特异效应。[分层效应](tables/context_distance.tsv) · [目标定位](tables/target_localization.tsv)

GSE285715 T2 的 TAM 梯度在候选内、外侧均为正 (+0.23、+0.24)；Treg 的负梯度在内侧更明显 (−0.19，外侧 −0.02)，R12 标签 1 中为 −0.16。Inhousedata 的 TAM 在内侧为 −0.21、外侧为 +0.16，R12 标签 2 和 3 中分别为 −0.22 和 −0.31。该片 P90 的方向翻转，说明其 TAM 关系对所选富集区域敏感。CAF 在 7 片中的 6 片为正，但多数较弱；NPJ 的 CAF 全片 r 为 +0.19，分开内外侧后分别为 −0.11、+0.07，提示分区构成本身足以影响整体方向。[区域构成](tables/region_composition.tsv) · [R12 与 local niche 覆盖和 MP10](tables/region_niche.tsv)

## target 的含义和覆盖

三条路线在 7 片均可定义，21 个主比较均有可用距离。各片每条路线入选 204–2,087 个点，占合格组织 20.00%–20.06%。各片 CAF 的 P80 阈值为全部参考量的 3.70%–35.12%，TAM 为 1.72%–10.60%，Treg 仅为 1.66%–3.80%。相应 target 内 posterior mean 丰度的片内中位数范围为 CAF 0.37–5.79、TAM 0.22–2.49、Treg 0.115–0.803。固定 P80 适合描述片内相对富集；不同切片的 target 不具有相同的绝对量。特别是 Treg，七片都属于低估计丰度的探索性结果。TAM 仅含 C1QC TAM 和 SPP1 TAM，FOLR2 无独立可评估因子。[定义、阈值及实际丰度](tables/targets.tsv)

245 个 query–target 组合在同一组织块无 target，距离为 NA；另有 20,323 个组合缺精确病灶边界距离，但有可用的目标距离。二者均按原定义处理。T2 有约 17%–19% 的 query 与 P80 target 同 spot；去零后其 CAF、TAM 正相关减弱。NPJ 的 TAM 和 Treg 仅正距离梯度保留，但仍可能受到上皮比例、组织结构和多个细胞共享 spot RNA 的影响。所有相关都是单片描述，没有将 spot 当作独立患者或将 7 片合并作显著性推断。[距离覆盖](tables/distance_coverage.tsv) · [实距分组及点数](tables/distance_bins.tsv)

## 功能程序

在各片完整 RNA 轴 (15,966–33,538 个基因) 上计算 50 个 Hallmark、3 个 Reactome 和 2 个预设 C7，385 个片–集组合均可评分。两种胆固醇程序与 Epi queries 的原 MP10 均为正相关：Hallmark 胆固醇稳态 r 约 +0.11 至 +0.63，Reactome 胆固醇生物合成约 +0.10 至 +0.65。两者作为代谢程序一致性参照有价值，但共享转录数据和基因成员。[逐片程序一致性](tables/cholesterol_MP10.tsv) · [成员匹配](tables/gene_set_coverage.tsv)

**下一批最值得保留的是 EMT–CAF 距离和胆固醇一致性。** EMT spot 程序在 7 片均靠近 CAF 时较高 (r 约 −0.47 至 −0.09)；Reactome ECM 在 6 片同向，其中前五片为 −0.47 至 −0.31。缺氧、E2F/G2M 和炎症程序适合作为逐片背景，缺氧–CAF 的 r 在 T2 为 −0.28、NPJ 为 +0.19。首两片里 UCell 与 AddModuleScore 的全组织 Spearman ρ，EMT 为 0.965、0.922，缺氧为 0.815、0.685，胆固醇生物合成为 0.787、0.615；三程序 × 三目标 × 两片的 18 个梯度方向均一致。这些是混合 spot 内的转录程序，不能据此归因于单个 CAF。[程序距离效应](tables/program_distance.tsv) · [两方法参照](tables/method_comparison.tsv)

两个 C7 程序分别保留 DN (activated 相对高表达) 和 UP (resting 相对高表达)。全组织 C7 DN 与 Treg 组成的 r 约 −0.25 至 +0.01，CTLA4 程序约 −0.36 至 +0.12；没有可用作每个 Treg 功能强度的稳定正对应。建议将两组原始程序与 Treg 估计量并列保留作空间背景，而不据此定义抑制功能。[全组织 Treg 程序与组成](tables/Treg_program_composition.tsv)

## 与作者实现的关系

本批复用了 Oral/cSCC 的目标中心到 query 中心最小欧氏距离定义、cSCC 的连续程序散点及 Lung 的区室并列展示。片内 P80、同组织块搜索、200 μm 目标限制、项目冻结的 anchor/MP10、R12 和标准基因集 UCell 均为本项目适配。原脚本与适配脚本分开保存。本批计算止于这 7 张切片。[方法和运行入口](METHODS.md)

