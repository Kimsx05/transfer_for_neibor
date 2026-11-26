# 结果判读与方法校验

## 主检验与跨片方向

D1_Johnson 的 50 µm high−low 密度差为 -3.13 个支持 bin/mm²，区组重采样 95% 区间 [-6.20, 0.24]。受限置换的中心为 0.81，双侧 P=0.0056，HD 族 q=0.0336。置换检验衡量相对分层零模型的偏离，区间描述切片内区组重采样的原始差值。上皮 high 的降深稳定性为低支持。

Xenium 平台 30 个主检验中最小 BH q=0.0948。A1R5 的三种 high 在 50 µm 均呈较低 Treg 密度，200 µm 内到固定 Treg 集合的限制平均距离均较长。A2R2 三种 high 的密度与 Treg/T 都较高、距离较短、≤2 µm 分割边界近接触比例较高，但全片主定义目标仅 13 个。A1R4 有 4 个，A1R6 有 5 个；A1R6 的 H05 高转移支持有局部正向信号，合成读数的密度方向可随单个 Treg 的剔除翻转。

## 三个 Xenium 状态

下表为五片共同合格上皮上的范围。连续秩相关与极端尾部交集分别报告。

| 状态 A | 状态 B | Jaccard 最小 | Jaccard 最大 | Spearman 最小 | Spearman 最大 |
|---|---|---|---|---|---|
| Chol_biosynthesis_8 | H05_transfer | 0.03109 | 0.1109 | 0.3273 | 0.4982 |
| MP10_covered | Chol_biosynthesis_8 | 0.2556 | 0.3815 | 0.6624 | 0.8332 |
| MP10_covered | H05_transfer | 0.025 | 0.1031 | 0.1483 | 0.465 |

## 原生工具、零模型与校验

| sample | state | density_delta | density_q | Z | OE | CRAWDAD_Z_100 | CRAWDAD_Z_800 | curve_q | spatstat_outside_global_band |
|---|---|---|---|---|---|---|---|---|---|
| A1R5 | MP10_covered | -3.014 | 0.0948 | -16.11 | 0.002966 | -2.727 | -8.11 | 0.008 | True |
| A1R5 | Chol_biosynthesis_8 | -2.889 | 0.1273 | -16.28 | 0.005921 | -1.876 | -7.765 | 0.008727 | True |
| A1R5 | H05_transfer | -3.142 | 0.0948 | -15.54 | 0 | -3.163 | -7.796 | 0.006 | True |

Squidpy 的目标 Z 对比切片全局身份打乱。CRAWDAD 在去重的焦点缓冲区内计算上皮状态相对局部打乱的 Z，允许不同类型坐标组成共同变化。冻结图的主检验只交换合格上皮内状态。数值和显著性不能互换。

spatstat::rlabel 的独立 4,999 次联合随机标记对 36 个主读数进行重算。观测差值最大绝对误差 1.16e-13。新随机化样本的 P 和 BH 另表保留。[原生受限置换核对](spatstat_primary_check.tsv)

CRAWDAD 的原生窗口与稀疏全片实现比较了 96 个类别/尺度，70 个有有限 Z。有限结果的最大绝对误差为 5.33e-15。原生 Z 遇到任一置换的合并比例为 0 或 1 时缺失，另列可定义置换次数和按有限值条件计算的 Z。[完整校验](CRAWDAD_kernel_checks.tsv)

## Niche 的生物学可解释范围

Cervilla 与 D1_Johnson 均有 Matrix_CAF 占首位的 MENDER 区，也有上皮、Mixed 和 Unknown 占首位的区。UC_SCC 及五片 Xenium 的首位类别常为上皮、Unknown 或 Mixed。身份未解析率与原始深度在不同区内变化。可复现的是这些组成差异，当前不足以给所有同号区赋予统一生物学命名。

| sample | state | niches | Epi_high_min | Epi_high_max | supported_Treg_min | supported_Treg_max |
|---|---|---|---|---|---|---|
| A1R2 | Chol_biosynthesis_8 | 4 | 0.0452 | 0.06767 | 3 | 24 |
| A1R2 | H05_transfer | 4 | 0.007519 | 0.05808 | 3 | 24 |
| A1R2 | MP10_covered | 4 | 0.03972 | 0.07519 | 3 | 24 |
| A1R4 | Chol_biosynthesis_8 | 4 | 0.0277 | 0.1728 | 0 | 4 |
| A1R4 | H05_transfer | 4 | 0.03047 | 0.1574 | 0 | 4 |
| A1R4 | MP10_covered | 4 | 0.04848 | 0.06996 | 0 | 4 |
| A1R5 | Chol_biosynthesis_8 | 4 | 0.005089 | 0.06769 | 0 | 135 |
| A1R5 | H05_transfer | 4 | 0 | 0.0689 | 0 | 135 |
| A1R5 | MP10_covered | 4 | 0.002545 | 0.06621 | 0 | 135 |
| A1R6 | Chol_biosynthesis_8 | 4 | 0.01899 | 0.1176 | 0 | 4 |
| A1R6 | H05_transfer | 4 | 0 | 0.09821 | 0 | 4 |
| A1R6 | MP10_covered | 4 | 0 | 0.1716 | 0 | 4 |
| A2R2 | Chol_biosynthesis_8 | 4 | 0.04773 | 0.129 | 0 | 6 |
| A2R2 | H05_transfer | 4 | 0 | 0.06325 | 0 | 6 |
| A2R2 | MP10_covered | 4 | 0.0452 | 0.08065 | 0 | 6 |
| Cervilla | MP10 | 6 | 0.007605 | 0.0608 | 1 | 40 |
| D1_Johnson | MP10 | 6 | 0.003906 | 0.05468 | 2 | 298 |
| UC_SCC | MP10 | 4 | 0.01741 | 0.0641 | 3 | 118 |

niche 内共有 84 个可评主定义次要检验，BH q<0.05 的项目为 0 个。high 分布和 Treg 所在 niche 的重叠仍属描述。

CellCharter 按固定预算迭代，38/48 次拟合在 40 轮内满足原生收敛判据。未收敛的结果作为有限配置的补充，解释还需结合种子 ARI 和屏蔽目标基因的敏感性。[收敛表](CellCharter_convergence.tsv)

## 最少必要人工任务

首批 150 个对象均保留原始表、独立坐标和局部图。16 个对象有 Codex 的逐图记录，其余候选为按冻结规则逐行完成的 table_only 判读。人工签署数为零。优先判断 D1_Johnson 混合 bin 中 FOXP3 与 T 谱系的归属，以及 A1R4/A1R6 稀少目标的边界与活化 T 竞争。Xenium 的六个 pan-T 候选仅 CD3E、CD3G 实测。完成这些身份核对后，以独立测量复核胆固醇合成程序；A1R5 和 A2R2 可作为方向相反的组织区域进行后续验证。[便携复核包](../review/first_batch.zip)
