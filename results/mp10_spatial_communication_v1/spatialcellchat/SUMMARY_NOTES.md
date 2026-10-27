# SpatialCellChat 两片试点重点解读

本说明仅比较各片内部的 H/L 方向和谱系排序；两片不是原 scRNA 五队列的复刻，不比较不同切片的绝对数值。GSE171351_A1 与 Inhousedata_T4469479 各一片，重复数不足以作患者差异检验。以下 mean 均为总 mean strength：所有纳入且经位置 QC 的 LR 逐个 F_l/K_l 后相加。

## Treg：主 mean 两片均下降，谱系份额与 flow 不完全同向

|切片|方向|Low mean|High mean|Low flow|High flow|T/NK排名 L→H|T/NK内份额 L→H|全TME份额 L→H|
|---|---|---:|---:|---:|---:|---|---|---|
|GSE171351_A1|Epi→TME|0.383295|0.379406|66.879|61.119|2/5→2/5|17.70%→16.85%|2.41%→2.34%|
|GSE171351_A1|TME→Epi|0.386759|0.382181|68.088|63.689|2/5→2/5|17.67%→16.85%|2.41%→2.33%|
|Inhousedata_T4469479|Epi→TME|0.204514|0.147416|106.585|120.582|4/5→5/5|13.71%→14.05%|2.43%→2.24%|
|Inhousedata_T4469479|TME→Epi|0.215975|0.156735|115.411|133.793|4/5→5/5|13.87%→14.22%|2.50%→2.30%|

两片两向 Treg 的主 mean 都是 High<Low。A1 中 Treg 在可辨 T/NK 层级排名第二，但不属于 high 状态下增量突出者；其 T/NK 内份额下降。Inhouse 的 high 状态 Treg mean 排名5/5，尽管绝对 mean 下降，T/NK 内份额反而略升——该片其他 T/NK 成分下降幅度不同。因此不应把“份额上升”写成“Treg绝对通讯增强”。T_NK_unresolved 汇合了三个不可独立辨识标签，保留于分母和排名只是组成事实，不能命名为新亚型。

Inhouse Treg 的 flow 与 mean 反向：Epi→Treg flow 106.585→120.582，而 mean 0.204514→0.147416；Treg→Epi flow 115.411→133.793，而 mean 0.215975→0.156735。这是 F 与每条LR的支持链接分母 K 口径不同所致，不是两套独立H/L网络。总 mean 是 sum_l(F_l/K_l)，一般不等于 sum_l F_l / sum_l K_l；表里的 summed_LR_supported_pair_denominators 仅作覆盖汇总，不能当主 mean 的单一分母。A1 的 mean 和 flow 两向都下降，故仅按 mean 的两片同向不能升级为总量的重复下降。

## EC、Fibroblast 与 Mural：哪些方向重复，哪些相反

在两个方向上都重复的主 mean 方向：Arterial EC、Venous EC、Lymphatic EC 两片均 High<Low；iCAF 两片均 High>Low；Pericyte 两片均 High<Low。

在两个方向上均出现切片反例：Capillary EC、vCAF、Matrix_CAF、rCAF、SMC 在 A1 中 High>Low，在 Inhouse 中 High<Low。没有把这些不一致条目剔除或称为稳定高状态候选。

EC 的绝对排序也不固定：A1 中 Venous EC 最大，Inhouse 中 Lymphatic EC 最大；Capillary EC 的正向结果只在首片出现。内皮内份额中，Venous EC 两片两向均下降；Arterial/Lymphatic 在 A1 下降而 Inhouse 上升，Capillary相反。因此绝对 mean 与谱系内份额必须并列看。参考中未独立建模 Tip/Cycling EC，本分支不制造其列。

Fibroblast 中 vCAF 在两片两向均为绝对 mean 最大，但其 H−L 方向不重复；iCAF 是这两片唯一两向均 mean 增加的 CAF 标签，且两向谱系内份额均上升，但只可作为两片试点的探索性重复。Matrix_CAF 的谱系份额在 A1 下降、Inhouse 上升，即使其绝对方向相反。rCAF 的份额变化很小或因方向不同，不据此升级生物学结论。

Mural 按原参考保留独立 major。A1 中 SMC 高于 Pericyte，Inhouse 中 Pericyte 高于 SMC；SMC 的 mean 方向跨片相反，Pericyte 两片均下降，但 Pericyte 在 Mural 内份额 A1 下降、Inhouse 上升。它们没有被强行改名为 CAF。

|谱系/亚群|方向|A1 Low→High mean|A1 Δ|Inhouse Low→High mean|Inhouse Δ|
|---|---|---|---:|---|---:|
|Endothelial / Arterial EC|Epi→TME|0.23301→0.22792|-0.00509|0.07976→0.06293|-0.01683|
|Endothelial / Arterial EC|TME→Epi|0.23584→0.23091|-0.00494|0.08351→0.06596|-0.01755|
|Endothelial / Capillary EC|Epi→TME|0.45335→0.49496|+0.04161|0.13292→0.09730|-0.03562|
|Endothelial / Capillary EC|TME→Epi|0.46330→0.51796|+0.05466|0.13879→0.10226|-0.03653|
|Endothelial / Lymphatic EC|Epi→TME|0.29645→0.28680|-0.00965|0.21072→0.17613|-0.03459|
|Endothelial / Lymphatic EC|TME→Epi|0.29644→0.28555|-0.01089|0.21092→0.17646|-0.03446|
|Fibroblast / Matrix_CAF|Epi→TME|0.47218→0.47923|+0.00705|0.43612→0.36526|-0.07086|
|Fibroblast / Matrix_CAF|TME→Epi|0.48329→0.49618|+0.01290|0.46653→0.39037|-0.07616|
|Mural / Pericyte|Epi→TME|0.90591→0.79903|-0.10688|0.24842→0.23918|-0.00924|
|Mural / Pericyte|TME→Epi|0.90520→0.80542|-0.09978|0.25062→0.24203|-0.00858|
|Mural / SMC|Epi→TME|1.86621→1.96842|+0.10221|0.13971→0.08458|-0.05513|
|Mural / SMC|TME→Epi|1.93138→2.01912|+0.08774|0.15367→0.09302|-0.06066|
|Endothelial / Venous EC|Epi→TME|1.63740→1.58626|-0.05114|0.11859→0.08715|-0.03144|
|Endothelial / Venous EC|TME→Epi|1.63295→1.59321|-0.03975|0.12371→0.09052|-0.03319|
|Fibroblast / iCAF|Epi→TME|0.20542→0.22216|+0.01674|0.11670→0.12084|+0.00414|
|Fibroblast / iCAF|TME→Epi|0.20870→0.22326|+0.01457|0.11752→0.12280|+0.00529|
|Fibroblast / rCAF|Epi→TME|0.21008→0.22571|+0.01563|0.10440→0.08752|-0.01688|
|Fibroblast / rCAF|TME→Epi|0.21235→0.22646|+0.01411|0.11036→0.09163|-0.01873|
|Fibroblast / vCAF|Epi→TME|2.00284→2.17043|+0.16759|1.02076→0.66096|-0.35980|
|Fibroblast / vCAF|TME→Epi|2.04261→2.20119|+0.15857|1.05960→0.69714|-0.36246|

## 首片有限范围敏感性

固定1037个LR、相同Q与参数，仅将 native diffusible 范围200→300μm；该版本ECM与分泌共用此范围，接触保持105μm。70个major/fine方向比较68项符号保留。两项接近零的翻转为 Epi→Mural（−0.004674→+0.000984）和 Epi→Classical Mono（+0.000334→−0.000174）。Treg两向mean仍下降；其余可辨EC、Fibroblast、Mural fine方向保留。每条Contact位置概率与主运行完全一致，已作独立数值核对。

所有结论均为 c2l 比例加权的跨spot通讯归属，依赖去卷积与混合表达；不能确定某个受体转录本一定来自同spot中的Treg，不能将这些数值称为免疫抑制强度、功能验证或因果机制。保留参考identifiability caution。另见 METHOD_NOTES.md、LR_category_summary.tsv、run_summary.tsv、numerical_checks.tsv。
