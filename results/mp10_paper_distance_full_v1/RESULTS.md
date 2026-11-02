# 63 张 Visium 切片的 MP10 空间距离结果

冻结 manifest 登记的 **63/63 张切片、7 个 dataset、265,002 个组织 spot** 全部完成逐点连接、分区组成及标准程序评分。共有 143,461 个冻结 Epi queries。原 7 张 pilot 各出现一次；P80、正距离、P90 和 margin 条件的 84 个 pilot Pearson r 与原结果逐项相同。全量中的 niche 中位数使用修正后的统计集合。[名单](full_sections.tsv) · [分来源覆盖](tables/cohort_coverage.tsv) · [总览](full_overview.pdf)

## 实际覆盖

| 模块 | 肿瘤研究来源 (60 张) | 邻近正常 (3 张) |
|---|---:|---:|
| 登记并整合的 spot | 259,374 | 5,628 |
| 主 Epi queries | 141,152 | 2,309 |
| 至少有一个有效 R12 标签的切片 | 58 | 3 |
| 全片 P80、仅正距离、P90 三路线均可估 | 58 | 3 |
| P80 margin 三路线均可估 | 55 | 0 |
| 固定调整三路线均可估 | 58 | 3 |
| 同配对集模型评分层三路线均可估 | 57 | 3 |
| 全部 55 个 UCell 程序可评分 | 60 | 3 |
| 全部 55 程序 × 3 距离效应可估 | 58 | 3 |

GSE238145_B13 无冻结 Epi query；GSE171351_D1 有 3 个 query，每条路线仅 2 个同块有效配对。两片的 Epi-query 距离相关记 NA，仍有全组织组成和功能图。GSE238145_B5 的主 MP10 可用，但没有同时满足两模型评分层有效性与距离条件的配对。三张邻近正常没有可见病灶外缘；其 9 个 margin 比较均不可用。CanDis_2026_Sample11、GSE319536_BC15 和 Inhousedata_T4463209 分别没有 CAF margin target，另外两条路线仍可估。margin 合计可报告 171/189 个逐片路线效应。分层或查询子集少于 20 对时不展示 r；niche 单元中少于 5 个有效 Epi queries 时不展示中位数。[不可估计效应明细](tables/unavailable_effects.tsv) · [逐片和分层覆盖](tables/distance_coverage.tsv)

## 三条路线的片内距离

全部 63 片均可定义三类 P80 富集目标。入选约占各片合格组织的 20.00%–20.10%。Treg 的 P80 阈值约为全部 30 因子组成的 0.20%–7.81%，各片目标内 posterior mean 丰度中位数为 0.023–1.451，其中 59/63 片低于 1。阈值和估计量逐片保留。主条件下可用逐片路线的零距离比例中位数为 11.3%，范围 0–53.1%。2,507 个 query–target 组合因同一组织块无目标而缺距离；另有 252,232 个组合缺精确病灶边界距离但可计算目标距离。[目标阈值与丰度](tables/targets.tsv) · [实距分箱](tables/distance_bins.tsv)

下表限肿瘤研究来源。每个 dataset 先取其可用切片 r 中位数，再对实际有结果的 dataset 等权作 Fisher-z 汇总。正值表示距离越远、分数越高。匹配模型两行使用完全相同的配对点。末列是可评估切片数，顺序固定为 CAF、TAM、Treg。

| 条件 | CAF | TAM | Treg | 切片数 |
|---|---:|---:|---:|---|
| Mixed-spot MP10，P80 全组织含零 | +0.126 | +0.065 | +0.109 | 58 / 58 / 58 |
| P80 全组织仅正距离 | +0.096 | +0.016 | +0.097 | 58 / 58 / 58 |
| P90 全组织含零 | +0.122 | +0.069 | +0.095 | 58 / 58 / 58 |
| P80 margin 含零 | +0.093 | +0.070 | +0.078 | 55 / 58 / 58 |
| 调整完整集合上的未调整 r | +0.126 | +0.065 | +0.109 | 58 / 58 / 58 |
| 固定组成、深度与位置调整后 | +0.057 | +0.020 | +0.037 | 58 / 58 / 58 |
| Matched-model spot MP10 | +0.107 | +0.055 | +0.096 | 57 / 57 / 57 |
| Matched-model Epi-component MP10 | +0.029 | +0.046 | +0.028 | 57 / 57 / 57 |

主 P80 的各 dataset 切片中位数范围分别为 CAF −0.072 至 +0.253，TAM +0.024 至 +0.135，Treg +0.033 至 +0.262。固定调整后，三条路线的数据集等权描述值均比同集合未调整值小；各片实际方向、同分区或 R12 内的关系及 Spearman 结果保留在表中。[逐片四条件](tables/associations.tsv) · [调整](tables/adjusted.tsv) · [模型层](tables/matched_layers.tsv) · [分层](tables/context_distance.tsv)

**身份覆盖影响汇总的对象。** 登记的肿瘤研究标本有 56 个既有单位，其中 29 个为核验患者；主距离可估切片中，55 张有可用单位，归入 54 单位、6 个 dataset。以 dataset 等权汇总各 dataset 的等权单位均值，CAF、TAM、Treg 分别为 +0.120、+0.062、+0.125。限核验患者时覆盖 28 张、27 位患者、3 个 dataset，对应 +0.031、+0.047、+0.074。P1 内 B5 和 B18 各贡献一个片内主效应，等权合并一次。GSE285715 缺少既有单位 ID，其切片进入逐片和 dataset 描述。全研究未使用 spot 数作为患者或单位权重。[单位及成员切片](tables/unit_effects.tsv.gz) · [直接等权单位结果](tables/unit_equal.tsv) · [各 dataset 分布](tables/dataset_effects.tsv.gz) · [数据集等权结果](tables/dataset_equal.tsv)

**邻近正常参照单列。** N1、N2、N3 的 P80 全组织 r 按 CAF/TAM/Treg 顺序分别为 (+0.337, +0.168, +0.114)、(+0.374, +0.181, +0.399)、(+0.210, +0.214, +0.130)。均保留原身份未解析状态。来源记录没有提供 T/N 配对，不据编号推断配对。

## 分区、niche 与功能

R12 覆盖 126,075/143,461 个 Epi queries；dataset-local 覆盖 265,002 个有效组织 spots。两种版本分别按各自合格且有标签的点计算分区构成。region × niche 的主 MP10 中位数统一改为该单元中有效 Epi queries 的中位数；表里同时给出组成点数、n_MP10_valid 和该 niche 的区域分布。以原 7 片的 dataset-local 单元核对，59 个原本可显示的中位数发生数值变化，另有 60 个原来显示的单元因有效 query 不足 5 个而改为 NA。v3.4 无可见外缘、无外侧路径及未定位的状态分别保留。[区域组成](tables/region_composition.tsv) · [niche 分母与中位数](tables/region_niche.tsv)

63 × 55 = **3,465 个 section–program 组合全部完成 UCell 评分**。各片完整 RNA universe 为 13,929–33,538 个基因；各片各集匹配比例 45.0%–100%，中位数 96.5%。区域表分别对有效组织和 Epi queries 报告 N、Q1、中位数及 Q3。Epi-query 程序距离共 10,395 个逐片路线组合，10,065 个有可估相关，其余 330 个来自上述两张查询不足切片。[匹配覆盖](tables/gene_set_coverage.tsv) · [程序距离](tables/program_distance.tsv) · [分区程序](tables/region_programs.tsv)

在肿瘤研究来源，Hallmark 胆固醇稳态和 Reactome 胆固醇生物合成与冻结 MP10 的数据集等权一致性 r 分别为 +0.396 和 +0.382 (59 片、7 个 dataset)。两种分数共享表达信息及部分基因。全组织中 resting-vs-activated Treg 的 DN、UP 和 Reactome CTLA4 程序与 Treg 组成的对应 r 分别为 −0.070、−0.029 和 −0.069 (60 片、7 个 dataset)。两个 C7 方向独立保留。其余程序均保留对应表和热图。报告限于混合 spot 转录程序、组成和空间效应的描述。[胆固醇程序](tables/cholesterol_MP10.tsv) · [Treg 组成参照](tables/Treg_program_composition.tsv)

