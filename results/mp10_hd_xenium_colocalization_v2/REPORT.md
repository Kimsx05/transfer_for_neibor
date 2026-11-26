# Visium HD 与 Xenium 的 MP10–Treg 空间关系

Provisional / Codex-reviewed，人工待签署。以下按切片报告。

## 执行与输入

原 8 µm 表达、成员、全 RCTD 权重、Xenium 原分割、参考转移及兼容 UCell 分数均已复用。运行核查包含 21 个状态/尺度组合。记录到的计算 scope 峰值为 3.71 GB，实际最高 7 个进程或线程。任务使用单 CPU 亲和性和单线程数值库，硬上限为每次 7 tasks / 165 GB，重计算阶段不重叠。峰值计入子进程与 cgroup 文件缓存。各阶段秒数及历史报错见 [资源表](reports/resources.tsv)。

## Treg 初审

| sample | candidates | supported | possible | unresolved | unsupported | reliable_T | old_supported | new_supported | human_requested |
|---|---|---|---|---|---|---|---|---|---|
| Cervilla | 906 | 55 | 76 | 593 | 182 | 1157 | 0 | 55 | 30 |
| D1_Johnson | 6772 | 782 | 1430 | 3820 | 740 | 18033 | 18 | 764 | 30 |
| UC_SCC | 6266 | 194 | 395 | 4747 | 930 | 12232 | 2 | 192 | 30 |
| A1R2 | 440 | 47 | 38 | 163 | 192 | 866 | 23 | 24 | 12 |
| A1R4 | 179 | 4 | 4 | 35 | 136 | 347 | 1 | 3 | 12 |
| A1R5 | 684 | 207 | 73 | 136 | 268 | 1084 | 170 | 37 | 12 |
| A1R6 | 182 | 5 | 5 | 26 | 146 | 300 | 1 | 4 | 12 |
| A2R2 | 89 | 13 | 14 | 20 | 42 | 501 | 10 | 3 | 12 |

三个 HD 数量为成分 bin，五个 Xenium 数量为分割细胞。低支持切片的个别目标可改变符号。[逐片留一表入口](reports/leave_one_summary.tsv)

| sample | old_Treg_not_top1 | old_missing_anchor | old_old_combination | old_parent_not_TNK |
|---|---|---|---|---|
| Cervilla | 43 | 3 | 0 | 52 |
| D1_Johnson | 631 | 23 | 0 | 663 |
| UC_SCC | 150 | 2 | 0 | 180 |

旧门控原因可以重叠。[全部入选入口](reports/candidate_routes.tsv) · [全部流失原因](reports/gate_loss.tsv)

## 新 16/32 µm MP10 与 Xenium 三状态

| sample | res | state | n | high | q95 | high_multigene | QC |
|---|---|---|---|---|---|---|---|
| Cervilla | 16 | MP10 | 37898 | 1893 | 0.09005 | 1 | low_support |
| Cervilla | 32 | MP10 | 10705 | 535 | 0.09818 | 1 | low_support |
| D1_Johnson | 16 | MP10 | 21831 | 1091 | 0.1222 | 1 | low_support |
| D1_Johnson | 32 | MP10 | 5526 | 277 | 0.1225 | 1 | evaluable |
| UC_SCC | 16 | MP10 | 64595 | 3230 | 0.09967 | 1 | evaluable |
| UC_SCC | 32 | MP10 | 16634 | 831 | 0.103 | 1 | evaluable |
| A1R2 | cell | MP10_covered | 17432 | 872 | 0.3427 | 1 | evaluable |
| A1R2 | cell | Chol_biosynthesis_8 | 17432 | 872 | 0.4372 | 1 | evaluable |
| A1R2 | cell | H05_transfer | 17432 | 872 | 0.7018 | — | evaluable |
| A1R4 | cell | MP10_covered | 4069 | 204 | 0.4006 | 1 | evaluable |
| A1R4 | cell | Chol_biosynthesis_8 | 4069 | 204 | 0.5762 | 1 | evaluable |
| A1R4 | cell | H05_transfer | 4069 | 204 | 0.5367 | — | evaluable |
| A1R5 | cell | MP10_covered | 21636 | 1080 | 0.3789 | 1 | evaluable |
| A1R5 | cell | Chol_biosynthesis_8 | 21636 | 1082 | 0.5276 | 1 | evaluable |
| A1R5 | cell | H05_transfer | 21636 | 1082 | 0.6668 | — | evaluable |
| A1R6 | cell | MP10_covered | 805 | 41 | 0.3061 | 1 | evaluable |
| A1R6 | cell | Chol_biosynthesis_8 | 805 | 41 | 0.4408 | 1 | evaluable |
| A1R6 | cell | H05_transfer | 805 | 41 | 0.5814 | — | evaluable |
| A2R2 | cell | MP10_covered | 3970 | 199 | 0.259 | 1 | evaluable |
| A2R2 | cell | Chol_biosynthesis_8 | 3970 | 199 | 0.3204 | 1 | evaluable |
| A2R2 | cell | H05_transfer | 3970 | 199 | 0.6371 | — | evaluable |

Cervilla 和 D1_Johnson 的 16 µm high 在 80% 稀释后的集合重叠低于冻结界限。连续排名和 high 稳定性分别解释。D1_Johnson 的 32 µm 分组达到稳定性门槛。三套 Xenium 状态在同片共享合格上皮。逐片重叠及相关见 [Jaccard](reports/high_overlap.tsv) 与 [得分相关](reports/score_correlations.tsv)。

## 50 µm 主检验

每行的密度与 Treg/T 差值都按 high 减 low，后者为比例单位。q 是平台内 BH 校正结果。

| sample | state_definition | targets | score_QC | effect_Treg/T | effect_density | FDR_Treg/T | FDR_density |
|---|---|---|---|---|---|---|---|
| A1R2 | Chol_biosynthesis_8 | 47 | evaluable | 0.008882 | -0.2947 | 0.9099 | 0.8178 |
| A1R2 | H05_transfer | 47 | evaluable | 0.04194 | -0.3826 | 0.0948 | 0.6072 |
| A1R2 | MP10_covered | 47 | evaluable | 0.002868 | -2.786 | 0.8662 | 0.1123 |
| A1R4 | Chol_biosynthesis_8 | 4 | evaluable | -0.0005607 | -0.1054 | 0.9674 | 0.787 |
| A1R4 | H05_transfer | 4 | evaluable | -0.000561 | -0.1054 | 0.787 | 0.6072 |
| A1R4 | MP10_covered | 4 | evaluable | -0.0005613 | -0.1054 | 0.9099 | 0.8178 |
| A1R5 | Chol_biosynthesis_8 | 207 | evaluable | -0.07761 | -2.889 | 0.4307 | 0.1273 |
| A1R5 | H05_transfer | 207 | evaluable | -0.1052 | -3.142 | 0.0948 | 0.0948 |
| A1R5 | MP10_covered | 207 | evaluable | -0.09713 | -3.014 | 0.108 | 0.0948 |
| A1R6 | Chol_biosynthesis_8 | 5 | evaluable | 0.143 | 4.96 | 0.1658 | 0.372 |
| A1R6 | H05_transfer | 5 | evaluable | 0.09492 | 10.65 | 0.1273 | 0.0948 |
| A1R6 | MP10_covered | 5 | evaluable | -0.01225 | -2.463 | 0.787 | 0.8178 |
| A2R2 | Chol_biosynthesis_8 | 13 | evaluable | 0.04452 | 5.602 | 0.1658 | 0.5673 |
| A2R2 | H05_transfer | 13 | evaluable | 0.03106 | 6.224 | 0.5576 | 0.4307 |
| A2R2 | MP10_covered | 13 | evaluable | 0.02145 | 3.966 | 0.5576 | 0.5576 |
| Cervilla | MP10 | 55 | low_support | -0.005806 | -0.06847 | 0.92 | 0.6249 |
| D1_Johnson | MP10 | 782 | low_support | -0.002249 | -3.134 | 0.4208 | 0.0336 |
| UC_SCC | MP10 | 194 | evaluable | 0.002262 | -0.3452 | 0.92 | 0.4208 |

通过平台内主检验 FDR 0.05 的项目：D1_Johnson MP10 density Δ=-3.13，q=0.0336。D1_Johnson 的二分类稳定性限制其解释。Xenium 的 A1R5 三定义原始密度差为负，A2R2 三定义为正，A1R6 随定义变化。逐片信号不能建立跨队列普遍正向富集。

## 暴露、距离与分割边界近接触

| sample | state_definition | density Δ | Treg/T Δ | Epi→Treg Δ µm | ≤2 µm 接触概率 Δ | Treg中心 high 暴露−零期望 | Treg→high 距离−零期望 |
|---|---|---|---|---|---|---|---|
| Cervilla | MP10 | -0.06847 | -0.005806 | -0.02284 | — | -0.03016 | 9.983 |
| D1_Johnson | MP10 | -3.134 | -0.002249 | -1.66 | — | -0.01475 | 8.631 |
| UC_SCC | MP10 | -0.3452 | 0.002262 | 2.76 | — | -0.008969 | 3.388 |
| A1R2 | MP10_covered | -2.786 | 0.002868 | 4.225 | 0.0005946 | -0.01478 | 1.415 |
| A1R2 | Chol_biosynthesis_8 | -0.2947 | 0.008882 | -0.2308 | -0.000597 | -0.0006962 | -2.24 |
| A1R2 | H05_transfer | -0.3826 | 0.04194 | 5.851 | -0.0005999 | 0.0007699 | 7.096 |
| A1R4 | MP10_covered | -0.1054 | -0.0005613 | -0.4344 | -0.0002618 | -0.03774 | -1.186 |
| A1R4 | Chol_biosynthesis_8 | -0.1054 | -0.0005607 | -0.4881 | -0.0002617 | -0.03887 | 17.92 |
| A1R4 | H05_transfer | -0.1054 | -0.000561 | -4.537 | -0.0002616 | -0.1774 | -1.14 |
| A1R5 | MP10_covered | -3.014 | -0.09713 | 17.3 | -0.0002937 | -0.0207 | 4.797 |
| A1R5 | Chol_biosynthesis_8 | -2.889 | -0.07761 | 17.33 | -0.0002937 | -0.006173 | 11.33 |
| A1R5 | H05_transfer | -3.142 | -0.1052 | 15.92 | -0.0002938 | -0.02135 | 14.81 |
| A1R6 | MP10_covered | -2.463 | -0.01225 | -14.32 | -0.001397 | 0 | -7.436 |
| A1R6 | Chol_biosynthesis_8 | 4.96 | 0.143 | -25.62 | 0.02778 | 0.09919 | -9.356 |
| A1R6 | H05_transfer | 10.65 | 0.09492 | -4.672 | 0.02439 | 0.4058 | -7.413 |
| A2R2 | MP10_covered | 3.966 | 0.02145 | -4.608 | 0.01717 | 0.00732 | -6.66 |
| A2R2 | Chol_biosynthesis_8 | 5.602 | 0.04452 | -6.775 | 0.006593 | -0.006361 | -6.08 |
| A2R2 | H05_transfer | 6.224 | 0.03106 | -13.78 | 0.001302 | 0.009236 | -1.044 |

负距离差表示更近。反向暴露与最近 high 距离的零模型保留每层 high 数量。近接触仅限 Xenium 的有效原边界。不同几何读数有不同有效分母；零分母及无可达目标见逐中心文件。[完整终点与区间](reports/endpoints.tsv.gz) · [有效覆盖与 T 分母](reports/denominators.tsv)

## 三个空间工具

Squidpy 全局标签打乱衡量固定图上的互斥身份聚集；HD 的独立成分属性通过专项统计保留。CRAWDAD 固定焦点、改变局部标签，spatstat 在上皮内局部分层移动状态。各工具对应不同零假设。方向与校验见 [逐片工具对照](reports/tool_comparison.tsv)。

## 三个 niche 方法

| sample | BuildNicheAssay | CellCharter all | CellCharter masked | MENDER |
|---|---|---|---|---|
| A1R2 | 0.9486 | 0.6199 | 0.6262 | 0.82 |
| A1R4 | 0.866 | 0.7751 | 0.7542 | 0.6212 |
| A1R5 | 0.9985 | 0.9987 | 1 | 0.8854 |
| A1R6 | 0.5909 | 0.6679 | 0.5808 | 0.679 |
| A2R2 | 0.8533 | 0.5268 | 0.9689 | 0.9321 |
| Cervilla | 0.8776 | 1 | 0.9979 | 0.846 |
| D1_Johnson | 0.5586 | 1 | 0.4715 | 0.6933 |
| UC_SCC | 0.9936 | 0.9192 | 0.9591 | 0.5412 |

| sample | MENDER 各区的首位身份及其占比 |
|---|---|
| A1R2 | 0 Epithelial 67%, 1 Epithelial 32%, 2 Epithelial 89%, 3 Unknown 46% |
| A1R4 | 0 Epithelial 46%, 1 Unknown 80%, 2 Epithelial 82%, 3 Unknown 72% |
| A1R5 | 0 Epithelial 87%, 1 Unknown 53%, 2 Epithelial 50%, 3 Unknown 27% |
| A1R6 | 0 Unknown 84%, 1 Unknown 63%, 2 Epithelial 78%, 3 Unknown 92% |
| A2R2 | 0 Unknown 77%, 1 Unknown 61%, 2 Unknown 93%, 3 Epithelial 58% |
| Cervilla | 0 Epithelial 95%, 1 Unknown 52%, 2 Unknown 96%, 3 Mixed 47%, 4 Epithelial 65%, 5 Matrix_CAF 35% |
| D1_Johnson | 0 Epithelial 88%, 1 Mixed 63%, 2 Unknown 95%, 3 Mixed 47%, 4 Epithelial 44%, 5 Matrix_CAF 60% |
| UC_SCC | 0 Mixed 42%, 1 Epithelial 96%, 2 Epithelial 48%, 3 Unknown 61% |

Mixed/Unknown 主导的区首先反映身份覆盖；深度与组织边缘也随 niche 变化。niche 在揭示代谢状态之前冻结，随后逐片叠加 high 份额、密度和 Treg/T。可在 [niche 覆盖与叠加](reports/niche_overlay.tsv)、[组成相似性](reports/niche_similarity.tsv) 与统一长表的 within MENDER niche 家族审阅。

## 人工复核与计算状态

首批 [浏览入口](review/index.html) 含三片 HD 各 30 个和五片 Xenium 各 12 个，共 150 个。优先复核 HD D1_Johnson 的非 T/NK 主要身份但有 Treg 成分者，以及 A1R4、A1R6 稀疏目标的边界和竞争线索。Xenium 的 Venous EC、SPP1 TAM、SMC 缺少原 anchor，各片按 panel-limited 保留。H05 为转移支持。现有数据支持优先校准 Treg 身份及高分尾部稳定性，之后检验候选组织区域中的具体关系。

[计算状态](reports/completion.tsv) · [六工具逐片状态](reports/tools.tsv) · [图件索引](reports/figures.tsv) · [方法与复现](METHODS.md)

[效应区间、三个状态、工具核验与 niche 的逐片判读](reports/INTERPRETATION.md) · [首批便携复核包](review/first_batch.zip)
