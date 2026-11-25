# MP10–SREBF–mTOR–NECTIN 虚拟扰动 v2

8,000 个原上皮细胞、38 位记录患者、39 个标本、5 个来源。两种模型使用同一观察性队列。scTenifoldKnk 完成 38 个正式 WT，C9 因仅 29 个细胞、2 位患者而记 NE。CellOracle 已完成小试、主种子、替代种子和五组按来源留一的簇级模型。正式簇任务共 63 个成功、7 个 NE。

## 去除核糖体和线粒体模块后的前排基因

在 P_full，五个主扰动的前 50 名中有 21–22 个细胞质核糖体结构蛋白及 0–8 个线粒体基因组基因。P_clean 精确去掉 39 个面板成员并重新建网。五个主扰动仍有 47 个基因全部进入前 50，任意两者的全排名 Spearman 为 0.979–0.998，前 50 两两重合 48–50 个。

| 扰动源 | scTenifold 前五（距离） | CellOracle 前五（患者等权绝对 Δ） |
| --- | --- | --- |
| SREBF1 | PRDX1、FOXN3、NME2、TPI1、KCNQ1OT1 | KRT13、KRT20、HSPA5、CENPF、MSMO1 |
| SREBF2 | NME2、PRDX1、TPI1、PPIA、EEF1G | HMGCS1、S100A7、GAPDH、SAA1、AQP3 |
| RPTOR | PRDX1、NME2、FOXN3、TPI1、KCNQ1OT1 | NE（非 TF） |
| MTOR | PRDX1、NME2、TPI1、FOXN3、KCNQ1OT1 | NE（非 TF） |
| RICTOR | NME2、PRDX1、TPI1、PPIA、EEF1G | NE（非 TF） |

PRDX1、NME2 均进入五个主源和全部 25 个匹配对照的前 50。S100A14 和 CALR 同时进入 SREBF2 与 RICTOR 的前 50；在此阈值未见仅某一主源独有的前排候选。scTenifold 适合提示共有响应，单靠其排名难以选出主源特异靶。

## A　SREBF 与 MP10/脂质转录程序

两种模型均给出 MP10 成员总体向下的预测，支持程序方向；特异性和成员层面证据有差别。scTenifold 主模型中，两源均为 20/20 向下，20 次患者子抽样均维持这一比例，位次会漂移。两组匹配对照的程序负向比例中位数也都是 100%。

| 源 | 程序 Full | CellOracle 可评估/冻结 | 负向成员 | 排名中位分位 | 可评估对照负向中位 | 五次留一负向范围 |
| --- | --- | --- | --- | --- | --- | --- |
| SREBF1 | MP10 | 20/20 | 85.0% | 8.2% | 40.0% | 85.0%–100.0% |
| SREBF1 | CHOL_SYNTH | 17/27 | 76.5% | 51.6% | 39.9% | 61.1%–88.9% |
| SREBF1 | SCENIC_SREBF1 | 8/21 | NE | NE | NE | NE |
| SREBF1 | SCENIC_SREBF2 | 28/36 | 75.0% | 34.9% | 64.3% | 71.4%–92.9% |
| SREBF2 | MP10 | 19/20 | 94.7% | 16.1% | 51.3% | 70.0%–100.0% |
| SREBF2 | CHOL_SYNTH | 10/27 | NE | NE | 50.0% | 41.2%–80.0% |
| SREBF2 | SCENIC_SREBF1 | 6/21 | NE | NE | NE | NE |
| SREBF2 | SCENIC_SREBF2 | 25/36 | 88.0% | 20.1% | 52.0% | 80.8%–96.4% |

覆盖门槛为至少 3 个成员且达到冻结集合的一半。NE 程序仍公开每个可达成员。M20 结果及成员见程序表；MP10 的 M20 无成员。

SREBF1 对 SCD 的 CellOracle 患者等权有符号 Δ 为 -0.0993，绝对变化排名 12/705。9 个主模型簇保留该直接推断边，其中 9 个系数为正。

SREBF2 对 HMGCS1 的 CellOracle 患者等权有符号 Δ 为 -0.2682，绝对变化排名 1/673。9 个主模型簇保留该直接推断边，其中 9 个系数为正。

## B　RPTOR、MTOR、RICTOR 与脂质程序

三者与 SREBF 相关脂质程序的负向转录轮廓相容。对 C1 特异性的支持不足。下表在五源共同可评估的非源 MP10 成员上比较，低排名分位代表更靠前。

| 源 | MP10 负向 | 成员排名中位分位 | 五对照排名中位分位 | 20次重拟合排名中位分位范围 |
| --- | --- | --- | --- | --- |
| SREBF1 | 100.0% | 6.4% | 6.5% | 4.9%–28.6% |
| SREBF2 | 100.0% | 6.6% | 6.6% | 5.3%–27.4% |
| RPTOR | 100.0% | 6.5% | 7.6% | 5.4%–29.6% |
| MTOR | 100.0% | 6.5% | 7.9% | 5.6%–38.2% |
| RICTOR | 100.0% | 7.0% | 7.1% | 5.5%–29.1% |

RPTOR 和 MTOR 未显示区别于 RICTOR 的一致程序选择性。mTORC1 Hallmark 含 200 个成员，各主源与各自对照的排名也接近。MTOR 同时参与两种复合体。LPIN1、RPS6KB1/2、EIF4EBP1 的基因行均在模块表；它们的 RNA 读出不能确定底物磷酸化。CellOracle 标准 TF 模型不评估这三个非 TF 源。

## C　SREBF 与 NECTIN/PVR

NECTIN2 在两种主模型中的有符号预测均向下。scTenifold 距离未通过 BH，CellOracle 的路径、幅度、对照和留一结果需要分开看。

| 源 | 终点 | sc 排名/2193 | sc q | CellOracle Δ | CO 排名 | 可达簇 | 低于主源的有效 TF 对照 | 留一 Δ<0/可评估 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SREBF1 | NECTIN2 | 178 | 1.000 | -0.0066 | 368 | 6/9 | 1/3 | 5/5 |
| SREBF1 | NECTIN1 | 415 | 1.000 | -0.0055 | 352 | 6/9 | 0/2 | 4/5 |
| SREBF1 | NECTIN3 | 196 | 1.000 | -0.0039 | 546 | 1/9 | 0/2 | 4/5 |
| SREBF1 | NECTIN4 | 367 | 1.000 | -0.0203 | 199 | 6/9 | 3/3 | 5/5 |
| SREBF1 | PVR | 877 | 1.000 | NE | NE | 0/9 | NE | 0/0 |
| SREBF2 | NECTIN2 | 171 | 0.883 | -0.0231 | 148 | 5/9 | 2/2 | 5/5 |
| SREBF2 | NECTIN1 | 479 | 1.000 | 0.0053 | 300 | 6/9 | 0/2 | 3/5 |
| SREBF2 | NECTIN3 | 208 | 1.000 | -0.0151 | 345 | 2/9 | 0/0 | 3/5 |
| SREBF2 | NECTIN4 | 355 | 1.000 | -0.0036 | 582 | 6/9 | 0/2 | 5/5 |
| SREBF2 | PVR | 949 | 1.000 | NE | NE | 0/9 | NE | 0/0 |

SREBF1 的 NECTIN2 在 scTenifold 原生、回归测序深度、再去 PC1 后的方向分数为 -0.00407、-0.00011、0.00058；距离超过 3/5 个匹配对照。CellOracle 在 38 位患者有可评估细胞，0 个主簇保留直接边，可达簇最短路径为 2、3 步。

SREBF2 的 NECTIN2 在 scTenifold 原生、回归测序深度、再去 PC1 后的方向分数为 -0.00855、-0.00299、-0.00125；距离超过 4/5 个匹配对照。CellOracle 在 38 位患者有可评估细胞，0 个主簇保留直接边，可达簇最短路径为 2、3 步。

映射至面板的官方启动子先验未列出 SREBF1/2→NECTIN2 直接边。两源的 NECTIN2 在 20 次 scTenifold 患者子抽样和 5 次来源留一中均保留原生负方向，均无 q<0.05。SREBF1 去除 PC1 后方向翻转。合并考虑，NECTIN2 获得弱的计算方向线索，直接调控和源特异性未获确证。SREBF2–NECTIN1 在两种主模型间方向相反；表中其余成员按各自覆盖与数值解释。

## 诊断、分层和可复核边界

SREBF1 的经典 manifold 对 NECTIN2 排名 175/2193，q=1.000；C7 独立 heat 模型（698 个细胞、28 位患者）排名 329/1203，方向 up。CellOracle 主簇合并后，2 步 -0.0043, 3 步 -0.0066, 4 步 -0.0071；官方随机图 3 步为 NE。各步和随机图的患者分母可不同。

SREBF2 的经典 manifold 对 NECTIN2 排名 142/2193，q=1.000；C7 独立 heat 模型（698 个细胞、28 位患者）排名 624/1203，方向 up。CellOracle 主簇合并后，2 步 -0.0205, 3 步 -0.0231, 4 步 -0.0181；官方随机图 3 步为 NE。各步和随机图的患者分母可不同。

冻结分层的 SREBF1/2–NECTIN2 配对 High 5%、10%、50% 与 Low 50% 分别有 23、34、35 位共同合格患者。两源各分层 High−Low 距离差的患者中位数均为正，个别患者为负。该距离反映固定 WT 加上层内表达分布的条件读出。

CellOracle 空扰动最大 |Δ| 为 0。主运行保留未额外裁剪的输出，实际越出 WT 范围和取零计数见 celloracle_simulation_qc。人工添加边为零，P/q 空白表示该算法未提供该统计量。

## 后续验证优先级

1. 优先检查 SREBF1–SCD、SREBF2–HMGCS1。在独立材料中区分转录响应与蛋白/脂质输出。这两对有较靠前的 CellOracle 位次及跨簇同号拟合边。

2. 若比较 RPTOR 与 RICTOR，应同步检查核内成熟 SREBP、复合体底物磷酸化及脂质读出。当前 RNA 模型无法区分 C1 与 C2 特异性。

3. NECTIN2 优先级低于上述脂质读出。若继续，分别测 RNA 与细胞表面蛋白，并检验直接调控证据。本数据未检验 NECTIN2–TIGIT 结合、Treg 抑制、胆固醇通量或药敏。

## 文件导航

主结果在 tables/。knk_all_gene_index.tsv 和 celloracle_cluster_index.tsv 定位每次真实拟合的全基因记录；celloracle_all_genes.tsv.gz 为跨簇患者等权排名。control_specificity、celloracle_specificity、cross_method_summary、mtor_module_summary 以及两种方法的 program 与 stability 表分别回答特异性、跨法、程序和重拟合问题。图见 FIGURE_INDEX.tsv，输入、范围和复现入口见 METHODS.md。
