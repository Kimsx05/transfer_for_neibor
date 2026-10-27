# MP10 round3：全TME总通讯、上皮状态比较与Treg定位

## 1. 全TME的大群与亚群总通讯格局

本轮完成模块0–4，主结果覆盖全部TME，不以B/Plasma、T/NK或旧AX列表限定。复用v2全部12个可用共享状态对象：五队列与ALL、major/fine各一套。实际纳入40个肿瘤标本、39个derived/provisional记录患者ID；Inhouse为12标本、11个记录ID。两类epi是每个标本内部的相对MP10状态，共享同一份TME。未重新评分、推断网络或跑丰度模型，v1/v2保持只读。

**全局总输出最突出的谱系是Fibroblast；Treg并不是全局总收发最强的节点。** 全局incoming、outgoing、与epi的强度、MP10状态差和通讯份额是不同问题，下面分开回答。W为全部内部P<0.05的已计算LR概率之和，C为支持LR记录数量，R为不按P筛选的原始和。这些总量没有患者差异显著性P。

### 1.1 读取对象与实际覆盖

| scope | network_level | n_nodes | n_computed_LR | n_specimens | n_recorded_patient_IDs |
|---|---|---|---|---|---|
| ALL | fine | 47 | 1874 | 40 | 39 |
| CNP0000460 | fine | 34 | 1771 | 7 | 7 |
| GSE326225 | fine | 32 | 1723 | 6 | 6 |
| HRA003620 | fine | 30 | 1744 | 7 | 7 |
| InhouseData | fine | 41 | 1815 | 12 | 11 |
| PRJNA662018 | fine | 40 | 1862 | 8 | 8 |
| ALL | major | 9 | 1774 | 40 | 39 |
| CNP0000460 | major | 9 | 1612 | 7 | 7 |
| GSE326225 | major | 9 | 1570 | 6 | 6 |
| HRA003620 | major | 9 | 1632 | 7 | 7 |
| InhouseData | major | 9 | 1743 | 12 | 11 |
| PRJNA662018 | major | 9 | 1710 | 8 | 8 |


major固定9个显示节点（两个epi＋七TME），fine固定47个显示节点（两个epi＋45个TME标签）。队列中不可计算节点显示NE，排名不把它们当0纳入。对象只覆盖其已计算LR维度，不代表所有可能分子联系。完整表保存TME↔TME、epi↔TME、epi↔epi及自环。

### 1.2 大群主格局

ALL的七个TME大群结果如下；排名仅在七个可计算TME内进行：

| node | incoming_all_W | rank_incoming_all_W | outgoing_all_W | rank_outgoing_all_W | self_W |
|---|---|---|---|---|---|
| B | 3.2417 | 4 | 0.85549 | 6 | 0.17218 |
| Plasma | 1.0789 | 7 | 0.48783 | 7 | 0.0027472 |
| T/NK | 3.5319 | 3 | 0.90467 | 4 | 0.28977 |
| Myeloid | 3.8025 | 2 | 1.5691 | 3 | 0.43659 |
| Fibroblast | 2.4182 | 5 | 10.883 | 1 | 1.2506 |
| Endothelial | 4.1628 | 1 | 6.1562 | 2 | 1.1679 |
| Mast | 2.1408 | 6 | 0.89308 | 5 | 0.0062981 |


Fibroblast的ALL outgoing W为10.8828，Endothelial为6.1562，Myeloid为1.5691。**这三个大群在五个dataset中均分别排outgoing第1、2、3**，比单一pooled图更有重复覆盖。incoming方面，Endothelial在四队列第1，HRA由T/NK居第1。B/Plasma/Mast等仍保留完整结果，没有因为主线改变而删除。

自环不能忽视：ALL Fibroblast的incoming中约51.7%来自大群自环；这是同一大群内不同细胞/亚型也被合并后的联系，不能称为单细胞自我作用。所有IO另给去自环版本。major和fine各自独立推断，不能把细群和当作大群，也不能把两层绝对W直接比较。

### 1.3 细群格局

ALL incoming前列的代表项（完整45群见表/图）：

| node | incoming_all_W | rank_incoming_all_W | n_matched_recorded_IDs |
|---|---|---|---|
| CD8 Exhausted | 48.976 | 1 | 29 |
| SPP1 TAM | 36.731 | 2 | 19 |
| Cycling EC | 36.456 | 3 | 2 |
| Tip EC | 36.36 | 4 | 8 |
| Arterial EC | 35.978 | 5 | 5 |
| NK / gamma-delta-boundary | 35.006 | 6 | 28 |
| dataset_bias_myeloid | 33.458 | 7 | 11 |
| FOLR2 TAM | 33.139 | 8 | 29 |


ALL outgoing前列的代表项：

| node | outgoing_all_W | rank_outgoing_all_W | n_matched_recorded_IDs |
|---|---|---|---|
| dataset_bias_fibro | 77.298 | 1 | 3 |
| Matrix_CAF | 76.787 | 2 | 11 |
| Cycling stromal | 70.842 | 3 | 1 |
| Pericyte | 65.702 | 4 | 15 |
| vCAF | 56.125 | 5 | 2 |
| iCAF | 55.554 | 6 | 18 |
| iCAF-like | 48.37 | 7 | 3 |
| rCAF | 45.117 | 8 | 3 |
| Cycling EC | 42.296 | 9 | 2 |
| Tip EC | 40.93 | 10 | 8 |


CD8 Exhausted的ALL incoming W为48.9758，排1/45；SPP1 TAM和多种内皮亚群也处于前列。outgoing则主要由成纤维/基质类标签占据。dataset_bias_fibro虽数值最高，但其实际细群网络只在Inhouse可计算、配对覆盖3个ID，且有既有偏倚标签，不能把这个排名升级为跨队列生物学候选。Cycling stromal、vCAF等也存在稀少覆盖；Matrix_CAF、iCAF、Pericyte更适合结合各自队列范围继续阅读。

全量circle保留所有正边与实际自环；fine图较密，精读时使用配套完整source×target热图。主图入口是聚合网络，不是单LR图。

## 2. 两类上皮对全部TME的总输出、总输入与群体差异

对同一对象内的同一组可计算TME，U_out只汇总epi→TME、U_in只汇总TME→epi，均排除epi↔epi及epi自环。以下fine总量每一行只在对应对象内比较：

| scope | direction | n_common_TME | U_low | U_high | U_difference |
|---|---|---|---|---|---|
| ALL | epi_to_TME | 45 | 11.346 | 14.319 | 2.972 |
| ALL | TME_to_epi | 45 | 17.782 | 26.08 | 8.297 |
| CNP0000460 | epi_to_TME | 32 | 10.415 | 13.994 | 3.5786 |
| CNP0000460 | TME_to_epi | 32 | 13.193 | 17.12 | 3.9273 |
| GSE326225 | epi_to_TME | 30 | 11.274 | 15.067 | 3.7929 |
| GSE326225 | TME_to_epi | 30 | 9.1043 | 17.064 | 7.9594 |
| HRA003620 | epi_to_TME | 28 | 12.485 | 12.263 | -0.22175 |
| HRA003620 | TME_to_epi | 28 | 6.6783 | 6.5491 | -0.12923 |
| InhouseData | epi_to_TME | 39 | 11.085 | 13.945 | 2.8601 |
| InhouseData | TME_to_epi | 39 | 8.7972 | 18.869 | 10.072 |
| PRJNA662018 | epi_to_TME | 38 | 9.1602 | 11.802 | 2.6423 |
| PRJNA662018 | TME_to_epi | 38 | 17.732 | 24.163 | 6.4308 |


ALL中，epi输出总W由11.3465增至14.3185（+26.2%），从TME接收的总W由17.7825增至26.0795（+46.7%）。五队列里CNP、GSE、Inhouse、PRJNA两个总方向均偏high；HRA两个方向的总体W均略偏low。major也呈同样的队列级方向格局，但绝对值不能跨层比较。

全部TME节点的对象内方向计数：

| scope | direction | n_computed | n_high_greater | n_low_greater | n_zero |
|---|---|---|---|---|---|
| ALL | TME_to_epi | 45 | 45 | 0 | 0 |
| ALL | epi_to_TME | 45 | 45 | 0 | 0 |
| CNP0000460 | TME_to_epi | 32 | 31 | 0 | 1 |
| CNP0000460 | epi_to_TME | 32 | 31 | 0 | 1 |
| GSE326225 | TME_to_epi | 30 | 29 | 0 | 1 |
| GSE326225 | epi_to_TME | 30 | 30 | 0 | 0 |
| HRA003620 | TME_to_epi | 28 | 18 | 10 | 0 |
| HRA003620 | epi_to_TME | 28 | 6 | 22 | 0 |
| InhouseData | TME_to_epi | 39 | 38 | 0 | 1 |
| InhouseData | epi_to_TME | 39 | 37 | 1 | 1 |
| PRJNA662018 | TME_to_epi | 38 | 38 | 0 | 0 |
| PRJNA662018 | epi_to_TME | 38 | 38 | 0 | 0 |


**ALL细群的两个方向均为45/45节点high增强。** 因而“某群high联系更强”在pooled结果中本身缺乏相对特异性，需要再看份额、状态差排名和队列覆盖。不同对象的原始W没有被直接平均成患者效应，ALL也不计作第六个独立队列。

ALL epi→TME的最大W差出现在Tip EC（+0.20844）、Capillary EC（+0.20380）、Cycling EC（+0.19325）、Arterial EC（+0.17044）等；Treg差为+0.03993，排37/45。反方向最大差包括Matrix_CAF（+0.85653）、dataset_bias_fibro（+0.80091）、Cycling stromal（+0.78106）、Pericyte（+0.69100）；Treg为+0.01887，排27/45。Tip EC只有Inhouse一个细群网络可计算，Cycling EC主要单队列覆盖，不能仅凭ALL差最大就称为最稳健。

## 3. Treg的全网络地位、与两类上皮的聚合联系和份额

### 3.1 Treg总体收发是否强？

同一shared对象只有一个Treg节点，下面全局值只有一套；没有high组Treg/low组Treg两套总量。排名是当前对象可计算TME中的名次：

| scope | incoming_all_W | outgoing_all_W | rank_incoming_all_W | rank_outgoing_all_W | n_TME_ranked | self_W | rank_incoming_other_W | rank_outgoing_other_W | n_matched_specimens | n_matched_recorded_IDs |
|---|---|---|---|---|---|---|---|---|---|---|
| ALL | 18.179 | 8.024 | 36 | 34 | 45 | 0.066342 | 34 | 34 | 34 | 33 |
| CNP0000460 | 7.8065 | 2.4394 | 22 | 28 | 32 | 0.010147 | 22 | 28 | 7 | 7 |
| GSE326225 | 19.192 | 14.95 | 24 | 25 | 30 | 0.083057 | 24 | 22 | 4 | 4 |
| HRA003620 | 9.8369 | 8.4021 | 18 | 27 | 28 | 0.022751 | 18 | 27 | 6 | 6 |
| InhouseData | 25.099 | 11.821 | 27 | 31 | 39 | 0.25057 | 23 | 31 | 12 | 11 |
| PRJNA662018 | 12.263 | 3.3495 | 25 | 33 | 38 | 0.0023781 | 25 | 32 | 5 | 5 |


ALL Treg incoming W=18.1792（36/45）、outgoing W=8.0240（34/45）；去自环后incoming=18.1128（34/45）、outgoing=7.9576（34/45）。自环W=0.06634，只占incoming的0.365%、outgoing的0.827%，所以较低排名主要不是自环造成的。其global C分别为823和430，是跨伙伴细胞群对累计的支持LR记录数，不是823/430个独立机制。

Inhouse Treg incoming=25.0988（27/39）、outgoing=11.8206（31/39）；去自环incoming排23/39、outgoing仍31/39。各队列都不能据此称Treg为“全局通讯最强”或“免疫抑制最强”。

### 3.2 Epi-high与Treg的总联系是否更强？

| scope | direction | low_W | high_W | delta_W | low_C | high_C | delta_R | W_R_direction_agree | rank_delta_W |
|---|---|---|---|---|---|---|---|---|---|
| ALL | epi_to_TME | 0.18372 | 0.22366 | 0.039933 | 14 | 18 | 0.039933 | True | 37 |
| ALL | TME_to_epi | 0.019563 | 0.038435 | 0.018871 | 5 | 7 | 0.023763 | True | 27 |
| CNP0000460 | epi_to_TME | 0.23676 | 0.32068 | 0.083916 | 10 | 12 | 0.083916 | True | 23 |
| CNP0000460 | TME_to_epi | 0.026933 | 0.035385 | 0.0084521 | 3 | 3 | 0.0086758 | True | 25 |
| GSE326225 | epi_to_TME | 0.33203 | 0.40836 | 0.076332 | 11 | 16 | 0.093466 | True | 24 |
| GSE326225 | TME_to_epi | 0.016449 | 0.029643 | 0.013194 | 4 | 6 | 0.053536 | True | 26 |
| HRA003620 | epi_to_TME | 0.24586 | 0.22404 | -0.021825 | 20 | 21 | -0.017619 | True | 23 |
| HRA003620 | TME_to_epi | 0.041799 | 0.043688 | 0.0018896 | 9 | 10 | 0.0025289 | True | 16 |
| InhouseData | epi_to_TME | 0.28556 | 0.33465 | 0.04909 | 10 | 14 | 0.035952 | True | 27 |
| InhouseData | TME_to_epi | 0.014872 | 0.041749 | 0.026876 | 3 | 5 | 0.047383 | True | 31 |
| PRJNA662018 | epi_to_TME | 0.20441 | 0.26352 | 0.059118 | 10 | 15 | 0.059118 | True | 24 |
| PRJNA662018 | TME_to_epi | 0.021 | 0.034399 | 0.013399 | 4 | 5 | 0.016406 | True | 25 |


epi→Treg在CNP、GSE、Inhouse、PRJNA偏high，HRA偏low；反方向Treg→epi在五队列均偏high。Treg这12个对象/方向比较的W和辅助R符号一致，但内部支持阈值会改变增量大小，不能把W的增加自动当作统计显著。

ALL epi→Treg：0.18372→0.22366，支持LR数14→18；Treg→epi：0.01956→0.03843，支持LR数5→7。Inhouse分别为0.28556→0.33465（10→14条）和0.01487→0.04175（3→5条）。HRA epi→Treg的W下降，但C从20增至21，直接说明“数量更多”和“总强度更高”不能混用。

### 3.3 增强是否相对偏向Treg？

| scope | direction | low_share_pct | high_share_pct | rank_low_W | rank_high_W | rank_delta_share |
|---|---|---|---|---|---|---|
| ALL | epi_to_TME | 1.6192 | 1.562 | 31 | 33 | 23 |
| ALL | TME_to_epi | 0.11001 | 0.14738 | 38 | 35 | 8 |
| CNP0000460 | epi_to_TME | 2.2733 | 2.2916 | 22 | 22 | 18 |
| CNP0000460 | TME_to_epi | 0.20414 | 0.20668 | 27 | 27 | 11 |
| GSE326225 | epi_to_TME | 2.945 | 2.7103 | 16 | 19 | 20 |
| GSE326225 | TME_to_epi | 0.18067 | 0.17372 | 25 | 26 | 17 |
| HRA003620 | epi_to_TME | 1.9692 | 1.8269 | 25 | 26 | 25 |
| HRA003620 | TME_to_epi | 0.62589 | 0.66709 | 24 | 24 | 17 |
| InhouseData | epi_to_TME | 2.576 | 2.3997 | 20 | 22 | 22 |
| InhouseData | TME_to_epi | 0.16905 | 0.22125 | 34 | 32 | 20 |
| PRJNA662018 | epi_to_TME | 2.2315 | 2.2328 | 22 | 23 | 19 |
| PRJNA662018 | TME_to_epi | 0.11843 | 0.14236 | 30 | 30 | 10 |


ALL epi分配给Treg的输出份额由1.6192%降至1.5620%，强度排名由31/45变为33/45；所以虽然绝对联系增强，**并没有相对更偏向Treg**。该输出份额只有CNP、PRJNA微升，GSE、HRA、Inhouse下降。Inhouse对应2.5760%→2.3997%。

反方向，Treg占epi全部TME输入来源的份额由0.1100%增至0.1474%，份额差排名8/45；五队列有四个上升、GSE下降。它在这个组成指标上比全局排名更突出，但绝对份额仍很小，而且与数据库反向分子记录相关，不能直接称为Treg特异反馈。

### 3.4 Treg总量主要来自谁、流向谁？

ALL四类来源/去向分解：

| direction | category | W | C | percent_of_direction_total |
|---|---|---|---|---|
| incoming | Epi_MP10_high | 0.22366 | 18 | 1.2303 |
| incoming | Epi_MP10_low | 0.18372 | 14 | 1.0106 |
| incoming | Treg_self | 0.066342 | 6 | 0.36493 |
| incoming | other_TME | 17.705 | 785 | 97.394 |
| outgoing | Epi_MP10_high | 0.038435 | 7 | 0.479 |
| outgoing | Epi_MP10_low | 0.019563 | 5 | 0.24381 |
| outgoing | Treg_self | 0.066342 | 6 | 0.82679 |
| outgoing | other_TME | 7.8996 | 412 | 98.45 |


Treg incoming的97.39%、outgoing的98.45%来自/流向其他TME；两类epi合计只占其incoming的2.24%、outgoing的0.72%。因此全局Treg总量与它的epi状态边不是同一个量。所有单独伙伴均在`treg_partner_breakdown.tsv`和`Treg_all_partners.pdf`中保留，没有仅展示强伙伴。

配对覆盖为34个标本、33个记录ID（CNP7、GSE4、HRA6、Inhouse11、PRJNA5个ID）；网络整体仍是40标本/39 ID。覆盖差异没有被当作删除全谱系网络的新门槛，也没有用这个数冒充患者级总W的重复。

## 4. 其他值得讨论的群体与反向/不一致队列

下表从全量结果中摘出便于对照的群体，不是另建加权候选分数；全部群仍在完整表和固定顺序方向矩阵中：

| node | direction | n_computable_datasets | n_high_greater | n_low_greater | n_zero_difference | n_datasets_at_least_two_matched_IDs | n_high_share | n_low_share |
|---|---|---|---|---|---|---|---|---|
| C1QC TAM | TME_to_epi | 5 | 5 | 0 | 0 | 5 | 4 | 1 |
| C1QC TAM | epi_to_TME | 5 | 4 | 1 | 0 | 5 | 2 | 3 |
| CD8 Exhausted | TME_to_epi | 5 | 4 | 1 | 0 | 5 | 2 | 3 |
| CD8 Exhausted | epi_to_TME | 5 | 5 | 0 | 0 | 5 | 2 | 3 |
| Capillary EC | TME_to_epi | 4 | 4 | 0 | 0 | 2 | 1 | 3 |
| Capillary EC | epi_to_TME | 4 | 4 | 0 | 0 | 2 | 3 | 1 |
| FOLR2 TAM | TME_to_epi | 5 | 5 | 0 | 0 | 5 | 3 | 2 |
| FOLR2 TAM | epi_to_TME | 5 | 4 | 1 | 0 | 5 | 2 | 3 |
| Matrix_CAF | TME_to_epi | 3 | 3 | 0 | 0 | 2 | 1 | 2 |
| Matrix_CAF | epi_to_TME | 3 | 3 | 0 | 0 | 2 | 3 | 0 |
| Pericyte | TME_to_epi | 4 | 4 | 0 | 0 | 3 | 2 | 2 |
| Pericyte | epi_to_TME | 4 | 4 | 0 | 0 | 3 | 4 | 0 |
| SPP1 TAM | TME_to_epi | 5 | 4 | 1 | 0 | 3 | 3 | 2 |
| SPP1 TAM | epi_to_TME | 5 | 4 | 1 | 0 | 3 | 1 | 4 |
| Tip EC | TME_to_epi | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| Tip EC | epi_to_TME | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| Treg | TME_to_epi | 5 | 5 | 0 | 0 | 5 | 4 | 1 |
| Treg | epi_to_TME | 5 | 4 | 1 | 0 | 5 | 2 | 3 |
| dataset_bias_fibro | TME_to_epi | 1 | 1 | 0 | 0 | 1 | 0 | 1 |
| dataset_bias_fibro | epi_to_TME | 1 | 1 | 0 | 0 | 1 | 1 | 0 |
| iCAF | TME_to_epi | 5 | 4 | 1 | 0 | 3 | 1 | 4 |
| iCAF | epi_to_TME | 5 | 4 | 1 | 0 | 3 | 4 | 1 |


可以分别提出以下具体问题：

- **Fibroblast/CAF是全局主要输出谱系。** 大群outgoing五队列均第1；Matrix_CAF反向TME→epi的ALL差最大，但细群只在3队列计算、其中2队列有至少2个匹配ID。iCAF更广泛，HRA存在反向；Pericyte在4个可计算队列中两个方向均high增强，配对覆盖较完整的有3队列。它们值得研究的理由与Treg不同，不能统称免疫抑制。
- **CD8 Exhausted是强接收节点。** ALL global incoming第1、epi→该群在五队列均high增强，五队列均有至少两个匹配ID；但输出份额只在两队列增加、三队列下降，ALL份额也下降，仍需区分绝对增量和相对偏好。
- **C1QC TAM/FOLR2 TAM是覆盖较完整的双向对照。** 五队列都可计算且有至少两个匹配ID；其→epi总量五队列均high增强，而epi→它们在HRA反向。这些是可以与Treg共同比较的细胞群，而非沿用v1的单LR shortlist。
- **Capillary EC有较大epi状态差。** 四个可计算队列中均偏high，但只有CNP、PRJNA达到至少两个同标本匹配ID。Tip EC、Cycling EC及部分CAF/dataset_bias标签更偏单队列或低覆盖，保留数值而不宣称普遍性。
- **B、Plasma、Mast、其他T/NK及未解析标签完整保留。** 例如B/Plasma并未成为本轮入口限制；其队列内结果与所有群使用相同W/C/R、排名和覆盖口径。

HRA是结构性不一致队列：epi总输出/输入均略低，但部分节点（例如CD8 Exhausted的epi接收方向、C1QC/FOLR2的→epi方向）仍偏high。不能把HRA整队列简单删除，也不应把缺失内皮/CAF节点记成“反向”。

## 5. 单条LR如何解释Treg的总量变化

每条贡献是支持阈值后的高低差，不是未过滤概率差；所有已计算LR都进入分解。ALL的正贡献代表项：

| direction | LR_id | supported_low | supported_high | contribution | raw_difference |
|---|---|---|---|---|---|
| epi_to_TME | NECTIN2_TIGIT | 0.010031 | 0.021017 | 0.010986 | 0.010986 |
| epi_to_TME | APP_CD74 | 0.022905 | 0.029928 | 0.0070233 | 0.0070233 |
| epi_to_TME | LAMA5_CD44 | 0 | 0.0036834 | 0.0036834 | 0.0036834 |
| epi_to_TME | MIF_CD74_CD44 | 0.062836 | 0.065836 | 0.0030003 | 0.0030003 |
| epi_to_TME | MIF_CD74_CXCR4 | 0.038513 | 0.040402 | 0.001889 | 0.001889 |
| epi_to_TME | JAM1_ITGAL_ITGB2 | 0.0020904 | 0.003876 | 0.0017857 | 0.0017857 |
| epi_to_TME | MDK_ITGA4_ITGB1 | 0.0084387 | 0.01014 | 0.001701 | 0.001701 |
| TME_to_epi | TIGIT_NECTIN2 | 0.010031 | 0.021017 | 0.010986 | 0.010986 |
| TME_to_epi | ADGRE5_CD55 | 0 | 0.0051231 | 0.0051231 | 0.00275 |
| TME_to_epi | NAMPT_INSR | 0.0050141 | 0.0060199 | 0.0010058 | 0.0010058 |
| TME_to_epi | CD46_JAG1 | 0 | 0.00086612 | 0.00086612 | 0.00086612 |
| TME_to_epi | CD96_NECTIN1 | 0.0023551 | 0.002779 | 0.00042396 | 0.00042396 |
| TME_to_epi | SEMA4D_PLXNB2 | 0.00080167 | 0.0011035 | 0.00030185 | 0.00030185 |
| TME_to_epi | CD6_ALCAM | 0.0013612 | 0.0015257 | 0.00016451 | 0.00016451 |


ALL epi→Treg有18条正贡献、无负贡献，总差+0.03993；NECTIN2_TIGIT贡献+0.01099，约27.51%，随后有APP_CD74 +0.00702、LAMA5_CD44 +0.00368、MIF_CD74_CD44 +0.00300等。NECTIN2在这个pooled对象是最大单条贡献，但在CNP、GSE、Inhouse的正贡献排名分别是5、6、3，不是每队列都第一。

ALL Treg→epi的TIGIT_NECTIN2贡献也是+0.01099，约占该方向净差58.21%。数据库的正反记录共享分子对，当前原始概率也相同；主总量按标准DB维度保留，但这不能当作两条独立重复或已证实双向功能反馈。另一个ADGRE5_CD55贡献+0.00512，其中low内部P=0.9未支持、high支持，故W贡献大于其原始概率差；这正是必须区分W差与R差的例子。

HRA的负贡献并未丢弃：

| LR_id | contribution |
|---|---|
| MDK_NCL | -0.004699 |
| COL4A5_CD44 | -0.004575 |
| APP_CD74 | -0.0044543 |
| MDK_ITGA4_ITGB1 | -0.0032138 |
| LGALS9_CD45 | -0.0023781 |


在HRA，NECTIN2_TIGIT仍增加约0.00149，但MDK_NCL、COL4A5_CD44、APP_CD74等负贡献更大，合计epi→Treg总W下降0.02183。**一条轴增强不能替代总通讯结论。** 完整LR分解和“其余正/负贡献”图已交付，贡献和与总差的误差处于浮点精度范围。

## 6. NECTIN2–TIGIT小型补充

此补充放在`supplementary/NECTIN2_TIGIT`，不替代全网络主图。真实单LR值如下：

| scope | epi_state | raw_probability | internal_P | display_strength | status |
|---|---|---|---|---|---|
| ALL | Epi_MP10_low | 0.010031 | 0 | 0.010031 | detected |
| ALL | Epi_MP10_high | 0.021017 | 0 | 0.021017 | detected |
| CNP0000460 | Epi_MP10_low | 0.017422 | 0 | 0.017422 | detected |
| CNP0000460 | Epi_MP10_high | 0.025811 | 0 | 0.025811 | detected |
| GSE326225 | Epi_MP10_low | 0.0087672 | 0 | 0.0087672 | detected |
| GSE326225 | Epi_MP10_high | 0.014679 | 0 | 0.014679 | detected |
| HRA003620 | Epi_MP10_low | 0.009771 | 0 | 0.009771 | detected |
| HRA003620 | Epi_MP10_high | 0.011259 | 0 | 0.011259 | detected |
| InhouseData | Epi_MP10_low | 0.0078834 | 0 | 0.0078834 | detected |
| InhouseData | Epi_MP10_high | 0.016552 | 0 | 0.016552 | detected |
| PRJNA662018 | Epi_MP10_low | 0.010469 | 0 | 0.010469 | detected |
| PRJNA662018 | Epi_MP10_high | 0.022158 | 0 | 0.022158 | detected |


所有五队列和ALL的low/high两端均有内部支持，圈图两条边都画，没有制造“high有、low无”。内部P=0是100次置换下的经验输出，不是患者差异P，也不说明真实概率精确为0。

同标本表达复核：

| dataset | n_specimens | n_recorded_IDs | n_positive | n_negative | n_zero | mean_expression_difference | median_expression_difference | LOO_min | LOO_max | LOO_direction_stable |
|---|---|---|---|---|---|---|---|---|---|---|
| CNP0000460 | 7 | 7 | 7 | 0 | 0 | 0.21492 | 0.17706 | 0.18525 | 0.2304 | True |
| GSE326225 | 4 | 4 | 4 | 0 | 0 | 0.16715 | 0.12406 | 0.088348 | 0.21722 | True |
| HRA003620 | 6 | 6 | 4 | 2 | 0 | 0.056719 | 0.068634 | 0.036092 | 0.075286 | True |
| InhouseData | 12 | 11 | 11 | 0 | 0 | 0.14149 | 0.14612 | 0.12554 | 0.1533 | True |
| PRJNA662018 | 5 | 5 | 5 | 0 | 0 | 0.18871 | 0.17453 | 0.14208 | 0.22405 | True |


共34个合格标本/33个记录ID；先把同ID重复标本的表达差平均后，31个ID的epi NECTIN2高状态均值较高，2个较低（均在HRA）。五队列的留一ID平均差均保持正向。检测支持分别记录：34标本中low端32、high端33个满足该状态NECTIN2及共享Treg TIGIT各≥10%，未要求两端同时通过才保留。保留真实零和低表达，不以细胞作为患者重复。

CD8 Exhausted及其他实际TIGIT阳性T/NK接收群的同口径表达/网络简洁对照在`supplementary_TIGIT_positive_receiver_comparison.tsv`。这说明分子轴可供下一轮具体核查；**没有验证Treg全网络总强度，也没有证明Treg特异接触或抑制。**

## 7. 当前结论与下一步具体问题

本轮支持的描述是：全网络主要输出来自成纤维/基质谱系；部分髓系、内皮和CD8 Exhausted具有较强接收或上皮联系；Treg整体收发排名偏后，但epi状态边在多数队列增强。这个增强通常伴随广泛TME增强，并不等于epi更偏向Treg。Treg输出份额和输入来源份额可出现不同方向，HRA保留反向。内部支持总和没有自动获得患者级差异显著性。

下一步应把问题拆开：首先在明确标本/患者单位下，检验CAF/基质输出、内皮/髓系与epi的广泛联系；若聚焦Treg，检验其相对于CD8 Exhausted、C1QC/FOLR2等完整对照的邻近性与边界位置，以及NECTIN2/TIGIT是否在对应细胞和空间关系中匹配。HRA的总量下降但单轴增加是优先解释的不一致，不宜通过改阈值或删队列消除。

本轮没有从pooled对象反推出患者总W，没有把LR表达乘积命名为患者CellChat总通讯，没有新空间/注释/CNV/预后分析。实际主结论已有五队列Treg覆盖说明，不需要为绘图或单一预期重推网络；稀少其他节点按覆盖限制描述。

阅读入口为[逐队列图索引](FIGURE_INDEX.md)。主circle是全LR聚合图A/B/C，W与C分别输出；完整热图、IO含/不含自环、角色散点、全群高低比较与份额、Treg面板均按major/fine及五队列/ALL输出。透明PDF和300dpi PNG保存在完整目录，`review.zip`仅含核心PDF、主要表、README、方法和本报告；不含旧对象、表达矩阵或PNG。
