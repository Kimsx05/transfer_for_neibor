# MP10 round1 v2 中文分析报告

## 本轮结论

本轮已按全上皮定义完成模块0–3；v1保持只读。**标本级 MP10 与 B/Plasma 丰度没有得到经多重校正的关联证据，但标本内 epi 状态分析给出可以进一步验证的 LR 线索。** 目前优先考虑 epi→B/Plasma 的 **APP–CD74**；B/Memory B→epi 的 **ADGRE5–CD55** 是方向相反的通讯路线候选，Plasma→epi 的 **CD46–JAG1** 保留为更有限的候选。它们不是已证实的接触、激活或抑制机制。

APP–CD74 主状态结果在四队列偏 epi-high，HRA003620 明确偏 epi-low；共同 B/Plasma 标本复核保留 CNP、Inhouse、PRJNA 的 high 方向及 HRA 的反向，但 GSE 只剩一个记录 ID且方向改变。不能把四队列初始同向直接表述为稳健四队列验证。Memory/Naive B 的公共队列覆盖比 GC B 更适合下一轮具体空间问题；Inhouse IFN-B 是保留标签的单队列探索对象，不能等同于公共队列 B 亚群。

主状态与标本分组是不同层次：前者比较同一标本内两类上皮与共享 TME；后者比较不同标本的全上皮平均 MP10。两者不一致不自动否定任一结果，也不能相互冒充独立验证。

## 1. 实际输入、样本和评分变化

有效细胞 **374,329**，全部上皮 **141,881**；58个标本、42个记录患者ID、5个dataset。15个Normal、43个Tumor。所有患者ID仍是derived/provisional，未经生物学身份核实。Inhouse有13个肿瘤标本、12个记录ID，其中12标本符合MP10和双epi状态条件；同一记录ID的两个标本在网络保留、表达复核先平均患者差，比例主模型采用固定一标本回退。

| dataset | 标本数 | 记录患者ID数 | 肿瘤标本 | 全上皮MP10合格 | 双epi状态合格 | 补算细胞 |
|---|---|---|---|---|---|---|
| CNP0000460 | 12 | 7 | 7 | 7 | 7 | 6827 |
| GSE326225 | 8 | 8 | 8 | 6 | 6 | 617 |
| HRA003620 | 14 | 7 | 7 | 7 | 7 | 1597 |
| InhouseData | 13 | 12 | 13 | 12 | 12 | 4556 |
| PRJNA662018 | 11 | 8 | 8 | 8 | 8 | 4585 |


123,699 个上皮复用历史分数，18,182 个按相同90,011-feature counts秩公式补齐。当前全部上皮均有可信有限分数，**无缺失补0**。60个跨队列重叠细胞 counts 完全一致，分数最大绝对误差小于1e-12。三份未进入标本 MP10 模型的肿瘤标本不是缺分数，而是epi少于30：GSE326225_BC5为22、GSE326225_BC7为19、T4466585为18；其细胞与其他适用描述/条件分析仍保留。

v1有27个MP10合格肿瘤标本，v2有40个；**新增13个，原合格标本6个二分标签翻转**。阈值是在各dataset新的合格肿瘤集合上一次固定，后续不随覆盖或分期重划。全部40个合格肿瘤标本也满足两类epi各≥30。

新增资格：

| dataset | orig.ident | old_n_malignant | n_epi_total | MP10_sample_allEpi | sample_MP10_group_allEpi |
|---|---|---|---|---|---|
| CNP0000460 | CNP0000460_P02T | 0 | 5443 | 0.2502 | high |
| CNP0000460 | CNP0000460_P04T | 2 | 7963 | 0.2073 | low |
| HRA003620 | HRA003620_BC-01 | 2 | 369 | 0.2572 | high |
| HRA003620 | HRA003620_BC-11 | 2 | 362 | 0.1967 | low |
| InhouseData | T4398343 | 14 | 940 | 0.2252 | high |
| InhouseData | T4450398 | 0 | 2038 | 0.175 | low |
| InhouseData | T4471438 | 9 | 825 | 0.1635 | low |
| PRJNA662018 | PRJNA662018_PRJNA662018_P1T1 | 25 | 1206 | 0.3264 | high |
| PRJNA662018 | PRJNA662018_PRJNA662018_P2T2 | 5 | 1779 | 0.1997 | low |
| PRJNA662018 | PRJNA662018_PRJNA662018_P3T3 | 3 | 722 | 0.2543 | high |
| PRJNA662018 | PRJNA662018_PRJNA662018_P4T4 | 1 | 1386 | 0.1313 | low |
| PRJNA662018 | PRJNA662018_PRJNA662018_P5T5 | 10 | 5238 | 0.1802 | low |
| PRJNA662018 | PRJNA662018_PRJNA662018_P6T6 | 0 | 10023 | 0.1635 | low |


原合格标本的分组翻转：

| dataset | orig.ident | old_MP10_sample | MP10_sample_allEpi | old_MP10_group | sample_MP10_group_allEpi |
|---|---|---|---|---|---|
| GSE326225 | GSE326225_BC3 | 0.1398 | 0.1398 | high | low |
| GSE326225 | GSE326225_BC4 | 0.1258 | 0.1406 | low | high |
| HRA003620 | HRA003620_BC-05 | 0.3083 | 0.306 | low | high |
| HRA003620 | HRA003620_BC-24 | 0.3603 | 0.1645 | high | low |
| InhouseData | T4091568 | 0.09461 | 0.2266 | low | high |
| PRJNA662018 | PRJNA662018_PRJNA662018_P8T8 | 0.2498 | 0.2496 | low | high |


原始均值、描述中位数、旧新阈值、计数、状态覆盖和排除原因均在 `tables/old_vs_new_samples.tsv` 与 `samples_v2.tsv`。图见 `figures/overall/old_vs_new_eligibility.pdf`。

## 2. B/Plasma 丰度与 B 内部组成

全局主分析使用asin√原比例；ALL调整dataset。40个标本经重复记录ID固定回退后实际拟合39标本；分期完整模型调整前后同为25标本。B parent>0 对应37标本，分期完整25；B≥30分期完整对应19标本，两种转换使用完全相同的这19标本。

以下为major/TME，连续效应按每1原始MP10单位报告，区间为95% CI：

| subtype | comparison | coefficient | CI_low | CI_high | q | n_specimens |
|---|---|---|---|---|---|---|
| B | sample_MP10_group_allEpi | -0.01018 | -0.1055 | 0.08512 | 0.8324 | 39 |
| Plasma | sample_MP10_group_allEpi | 0.01015 | -0.0863 | 0.1066 | 0.8324 | 39 |
| B | MP10_continuous | -0.06942 | -0.7222 | 0.5834 | 0.8382 | 39 |
| Plasma | MP10_continuous | -0.2383 | -0.8947 | 0.4181 | 0.8382 | 39 |
| B | MP10_stage_complete_unadjusted | -0.1954 | -0.717 | 0.3261 | 0.5337 | 25 |
| Plasma | MP10_stage_complete_unadjusted | -0.2655 | -0.9889 | 0.4579 | 0.5337 | 25 |
| B | MP10_stage_adjusted | -0.2702 | -0.8412 | 0.3008 | 0.3938 | 25 |
| Plasma | MP10_stage_adjusted | -0.4129 | -1.181 | 0.3556 | 0.3888 | 25 |


B/TME的sample high−low系数约−0.0102、q=0.832；Plasma/TME约+0.0102、q=0.832。连续、分期调整和all_cells分母均未支持确定性关联；全体谱系的ALL主MP10模型也无q<0.05结果。非显著不等于证明无关系；当前样本量和队列差异限制精度。

**Inhouse与公共队列没有统一丰度方向。** Inhouse B的二分与连续系数均正；CNP、GSE、HRA均负，PRJNA为正。Plasma在Inhouse二分/连续均负；公共队列二分有正有负，连续只有PRJNA为正，其余负。逐队列结果均未提供B/Plasma的校正显著关联。原始点、低B标记及连续散点已逐队列输出。

B内部组成的分期调整敏感性：

| subtype | analysis_variant | coefficient | CI_low | CI_high | q | n_specimens |
|---|---|---|---|---|---|---|
| GC B | primary_asin | -0.1802 | -0.7868 | 0.4265 | 0.968 | 25 |
| IFN-stimulated B | primary_asin | -0.1538 | -0.4593 | 0.1516 | 0.9666 | 25 |
| Memory B | primary_asin | -0.4149 | -0.9105 | 0.08071 | 0.8123 | 25 |
| Naive B | primary_asin | 0.4988 | -0.05869 | 1.056 | 0.8123 | 25 |
| GC B | parent30_asin | -0.2877 | -0.767 | 0.1916 | 0.881 | 19 |
| IFN-stimulated B | parent30_asin | -0.1188 | -0.5623 | 0.3248 | 0.9024 | 19 |
| Memory B | parent30_asin | -0.2932 | -0.7668 | 0.1804 | 0.881 | 19 |
| Naive B | parent30_asin | 0.4157 | -0.04973 | 0.8811 | 0.881 | 19 |
| GC B | parent30_logit | 0.7342 | -5.519 | 6.988 | 0.9927 | 19 |
| IFN-stimulated B | parent30_logit | 0.3137 | -4.525 | 5.153 | 0.9927 | 19 |
| Memory B | parent30_logit | 0.7171 | -3.988 | 5.423 | 0.9558 | 19 |
| Naive B | parent30_logit | 4.149 | -0.2421 | 8.541 | 0.9558 | 19 |


主asin中Memory B负、Naive B正；限制B≥30后仍如此，但换成相同标本的logit，Memory B、GC B和IFN-B等部分方向改变。所有区间均跨0且q不支持确定性结论。这说明低/稀疏组成和转换口径影响解释，不能择取更好看的方向。Naive B正向也只能是待验证趋势。二分、连续及分期完整调整前/后的三种分析完整见 `B_sensitivity_effects.pdf` 和模型表，不能用敏感性替代主分析。

### B标签的可比性核查

| dataset | GC B | IFN-stimulated B | Memory B | Naive B |
|---|---|---|---|---|
| CNP0000460 | 217 | 0 | 1400 | 1041 |
| GSE326225 | 24 | 1 | 445 | 1562 |
| HRA003620 | 468 | 2 | 4004 | 1103 |
| InhouseData | 24 | 2753 | 3 | 5 |
| PRJNA662018 | 78 | 0 | 1141 | 1011 |


Inhouse 2,785个B中2,753个标为IFN-stimulated B；公共队列该标签仅GSE1个、HRA2个。Inhouse Memory B仅3个、Naive B仅5个、GC B24个。标签本身高度队列相关，无法靠dataset调整完全恢复同一亚群的可比生物学分布。

冻结marker面板保留B身份支持；在Inhouse IFN-B中可见B身份表达，但ISG15/STAT1/MX1相对其他标签/队列的模式并不构成充分、独立的跨队列一致IFN定义。公共队列那1–2个IFN-B不能当作验证集。marker点图已标细胞和标本数；基因是否在母RNA存在、逐队列是否有非零证据另表保存。此核查既不证明标签错误，也不证明其比例差是生物学差异；沿用冻结注释并保留原技术警告。

## 3. 主状态通讯与患者表达复核

主状态有12个独立推断对象（5dataset＋ALL，major/fine各一），次级条件有56个，共68个；另有6个共同B/Plasma覆盖敏感性对象。大群从原始大群抽样重推，细群独立抽样重推。两类epi逐标本等量、共享TME只出现一次。ALL仅概览；未支持的已算概率仍保存原值，LR筛掉、节点不足或测量不明均不能当0。

下面是用于阅读的10条代表轴，不是“显著top10”。同向dataset数要求每队列≥2个表达支持记录ID、网络原始差与配对表达中位方向一致，且该队列至少一状态有内部支持。该数仍只是覆盖描述，不是差异P。

| axis_id | network_level | TME_type | direction | LR_id | n_datasets_same_direction | n_supported_patients | support_datasets |
|---|---|---|---|---|---|---|---|
| AX01 | major | B | epi_to_TME | APP_CD74 | 4 | 24 | CNP0000460;GSE326225;InhouseData;PRJNA662018 |
| AX02 | major | B | TME_to_epi | ADGRE5_CD55 | 3 | 23 | CNP0000460;HRA003620;InhouseData |
| AX03 | major | Plasma | epi_to_TME | APP_CD74 | 4 | 26 | CNP0000460;GSE326225;InhouseData;PRJNA662018 |
| AX04 | major | Plasma | TME_to_epi | CD46_JAG1 | 2 | 23 | CNP0000460;PRJNA662018 |
| AX05 | fine | Memory B | epi_to_TME | APP_CD74 | 3 | 14 | CNP0000460;GSE326225;PRJNA662018 |
| AX06 | fine | Memory B | TME_to_epi | ADGRE5_CD55 | 2 | 14 | CNP0000460;HRA003620 |
| AX07 | fine | Naive B | epi_to_TME | APP_CD74 | 3 | 7 | CNP0000460;GSE326225;PRJNA662018 |
| AX08 | fine | Naive B | TME_to_epi | ADGRE5_CD55 | 1 | 7 | CNP0000460 |
| AX09 | fine | IFN-stimulated B | epi_to_TME | APP_CD74 | 1 | 9 | InhouseData |
| AX10 | fine | IFN-stimulated B | TME_to_epi | COL4A3_SDC4 | 1 | 9 | InhouseData |


GC B未列为候选：合格同标本覆盖只见HRA和PRJNA各1个记录ID，没有任何dataset达到2个；两者方向还相反，**不能据此排序**。其完整状态和表达表仍保留。

四条大群轴的逐队列患者方向与网络差如下。患者正/负计数基于全部可测、足覆盖的配对表达，`n_supported_patients`另计≥10%检测支持；不把表达辅助乘积当CellChat概率。

| axis_id | dataset | n_patients | n_supported_patients | n_positive | n_negative | raw_delta | supported_median_delta | coverage_concordant |
|---|---|---|---|---|---|---|---|---|
| AX02 | CNP0000460 | 5 | 5 | 5 | 0 | 0.003742 | 0.07303 | True |
| AX02 | GSE326225 | 2 | 2 | 2 | 0 | 0 | 0.03309 | False |
| AX02 | HRA003620 | 5 | 5 | 4 | 1 | 0.002224 | 0.07131 | True |
| AX02 | InhouseData | 9 | 8 | 9 | 0 | 0.005815 | 0.1748 | True |
| AX02 | PRJNA662018 | 3 | 3 | 1 | 2 | 0.001877 | -0.006443 | False |
| AX01 | CNP0000460 | 5 | 5 | 5 | 0 | 0.04794 | 0.3282 | True |
| AX01 | GSE326225 | 2 | 2 | 2 | 0 | 0.05556 | 0.2871 | True |
| AX01 | HRA003620 | 5 | 5 | 1 | 4 | -0.009535 | -0.1983 | True |
| AX01 | InhouseData | 9 | 9 | 8 | 1 | 0.0469 | 0.368 | True |
| AX01 | PRJNA662018 | 3 | 3 | 3 | 0 | 0.01506 | 0.603 | True |
| AX04 | CNP0000460 | 4 | 4 | 4 | 0 | 0.0002552 | 0.006223 | True |
| AX04 | GSE326225 | 3 | 3 | 2 | 1 | 0 | 0.003832 | False |
| AX04 | HRA003620 | 5 | 3 | 1 | 4 | -0.0003072 | -0.01434 | True |
| AX04 | InhouseData | 9 | 7 | 7 | 2 | 0 | 0.007947 | False |
| AX04 | PRJNA662018 | 6 | 6 | 5 | 1 | 7.02e-06 | 0.005857 | True |
| AX03 | CNP0000460 | 4 | 4 | 4 | 0 | 0.01821 | 0.08342 | True |
| AX03 | GSE326225 | 3 | 3 | 2 | 1 | 0.03358 | 0.1384 | True |
| AX03 | HRA003620 | 5 | 5 | 1 | 4 | -0.002995 | -0.05578 | True |
| AX03 | InhouseData | 9 | 9 | 8 | 1 | 0.01661 | 0.0773 | True |
| AX03 | PRJNA662018 | 6 | 5 | 6 | 0 | 0.007994 | 0.227 | True |


**APP–CD74：** 比较的是同标本epi APP表达差及共享B/Plasma CD74端。B、Plasma和B亚群共享同一个epi APP差，因此这些多个接收群不是多次独立验证。HRA的低状态增强在配对表达中也出现，必须保留。Memory B与Naive B的主要公共队列支持来自CNP、GSE、PRJNA；HRA Memory B反向。Inhouse IFN-B只在该队列有足够覆盖。

**ADGRE5–CD55：** 反方向比较B/Memory B的共享ADGRE5端与epi CD55差。B大群在CNP、HRA、Inhouse具较完整同向支持；GSE网络无状态差，PRJNA网络正而患者表达中位方向负，不应称为跨五队列一致。Naive B反方向只有CNP满足同向覆盖，等级更低。

**Plasma→epi：** CD46–JAG1在CNP和PRJNA具high方向覆盖，HRA则偏low；它作为独立于APP接收路线的有限候选保留。Plasma的ADGRE5–CD55即使某些原始概率差为正，也不能在该队列内部支持不足时算作强支持；完整表保留它，而不强行与B大群共用最佳反向轴。

IFN-B的COL4A3–SDC4仅为Inhouse单队列探索；标签与测量可比性限制依然适用。所有候选都须回到具体发送/接收基因及空间位置验证，不能由通讯方向推断激活、抑制或直接接触。

### 一次共同B/Plasma标本限制复核

两类epi、B和Plasma各≥30后，共21标本、20记录ID。该复核比单目标覆盖更严格，GSE只剩1个ID，结果只能描述。大群代表轴结果：

| axis_id | scope | raw_delta | n_specimens | n_patients |
|---|---|---|---|---|
| AX01 | ALL | 0.01976 | 21 | 20 |
| AX02 | ALL | 0.004822 | 21 | 20 |
| AX03 | ALL | 0.008003 | 21 | 20 |
| AX04 | ALL | 0.0004325 | 21 | 20 |
| AX01 | CNP0000460 | 0.06284 | 4 | 4 |
| AX02 | CNP0000460 | 0.005338 | 4 | 4 |
| AX03 | CNP0000460 | 0.0248 | 4 | 4 |
| AX04 | CNP0000460 | 0 | 4 | 4 |
| AX01 | GSE326225 | -0.00788 | 1 | 1 |
| AX02 | GSE326225 | 0.002618 | 1 | 1 |
| AX03 | GSE326225 | -0.00423 | 1 | 1 |
| AX04 | GSE326225 | 0 | 1 | 1 |
| AX01 | HRA003620 | -0.0085 | 5 | 5 |
| AX02 | HRA003620 | 0.001994 | 5 | 5 |
| AX03 | HRA003620 | -0.002661 | 5 | 5 |
| AX04 | HRA003620 | 0 | 5 | 5 |
| AX01 | InhouseData | 0.05303 | 8 | 7 |
| AX02 | InhouseData | 0.0006176 | 8 | 7 |
| AX03 | InhouseData | 0.01959 | 8 | 7 |
| AX04 | InhouseData | 0 | 8 | 7 |
| AX01 | PRJNA662018 | 0.02084 | 3 | 3 |
| AX02 | PRJNA662018 | 0.003164 | 3 | 3 |
| AX03 | PRJNA662018 | 0.01135 | 3 | 3 |
| AX04 | PRJNA662018 | 5.475e-05 | 3 | 3 |


APP–CD74在CNP、Inhouse、PRJNA仍偏high、HRA仍偏low；GSE变为low方向且仅一个ID，提示不能忽视未匹配epi对原pool的贡献。B反方向ADGRE5–CD55总体仍偏high，但是否有内部支持和样本表达一致性需一起读。Plasma反方向CD46–JAG1在共同覆盖限制后，仅PRJNA保留非零状态差，其余队列差为已计算的0，因此不宜升级为稳健多队列主线。不同敏感性对象重新归一化，不能把绝对差大小变化解释为效应放大倍数。

## 4. 与次级标本分组、组织状态和T分期的关系

| axis_id | scope | epi_state_within_specimen | sample_MP10_group_allEpi |
|---|---|---|---|
| AX01 | CNP0000460 | 0.04794 | -0.05502 |
| AX01 | GSE326225 | 0.05556 | 0.02173 |
| AX01 | HRA003620 | -0.009535 | -0.01145 |
| AX01 | InhouseData | 0.0469 | 0.05415 |
| AX01 | PRJNA662018 | 0.01506 | -0.05281 |
| AX02 | CNP0000460 | 0.003742 | -0.007172 |
| AX02 | GSE326225 | 0 | 0.004059 |
| AX02 | HRA003620 | 0.002224 | 0.007804 |
| AX02 | InhouseData | 0.005815 | 0.01042 |
| AX02 | PRJNA662018 | 0.001877 | 0.005145 |
| AX03 | CNP0000460 | 0.01821 | -0.0307 |
| AX03 | GSE326225 | 0.03358 | -0.01085 |
| AX03 | HRA003620 | -0.002995 | -0.01275 |
| AX03 | InhouseData | 0.01661 | 0.01753 |
| AX03 | PRJNA662018 | 0.007994 | -0.02737 |
| AX04 | CNP0000460 | 0.0002552 | 0.0002481 |
| AX04 | GSE326225 | 0 | 0 |
| AX04 | HRA003620 | -0.0003072 | -0.00171 |
| AX04 | InhouseData | 0 | 0.001173 |
| AX04 | PRJNA662018 | 7.02e-06 | -0.000179 |


表内两列分别是标本内epi状态差和不同标本MP10分组差，不能混称high/low。例如APP–CD74的B大群ALL主状态差为+0.0102，而标本MP10组差为−0.0205；CNP和PRJNA也由状态high增强转为标本组high较弱，Inhouse与GSE同为正、HRA同为负。APP–CD74可在同标本状态中偏high，而在不同标本均值分组时反向或不一致；这反映比较单位、细胞组成、受体背景及独立对象尺度均不同，不能声称“同一效应被验证/推翻”。每个队列的热图和独立对象circle还给出Normal/Tumor和T1/T2，未计算的LR为NE。

只有CNP同时提供T1和T2的队列内条件；其他队列缺一条件或缺明确T，不能制造组间比较。ALL的T比较也受dataset和样本构成限制。Normal/Tumor有配对记录ID补充，但Normal来源和治疗资料缺失，解释应保留这一限制。

Treg、FOLR2/C1QC/SPP1 TAM、CAF、内皮及其他T/NK/髓系均保留原比例模型、独立网络、IO与完整数值边。`all_groups_evidence_summary.tsv`与`all_TME_state_LR_descriptive.tsv`提供适用汇总；后者只是每对象内数值差较大的LR，不是完成了同等患者表达验证的新生物学优先排名。本轮不继承v1的Treg/FOLR2 shortlist，也不因某一pooled强边覆盖B/Plasma研究重点。

## 5. 下一轮可检验的空间问题

本轮只读取现有accepted SAMPLE_R12_LEIDEN标签与轻量组成摘要，未做空间重分析。既有粗粒度标签是GC B、Naive/Memory B、Plasma，不能直接等同于本轮细分Memory、Naive或IFN-B。读取到的相关摘要为：

| accepted_label | niche | celltype | mean_unit_fraction | mean_unit_difference | n_units |
|---|---|---|---|---|---|
| 1 | niche1 | GC B | 0.1105 | 0.05315 | 39 |
| 1 | niche1 | Naive/Memory B | 0.02982 | -0.01112 | 39 |
| 1 | niche1 | Plasma | 0.01406 | -0.006101 | 39 |
| 2 | niche2 | GC B | 0.04827 | -0.01326 | 55 |
| 2 | niche2 | Naive/Memory B | 0.05615 | 0.03291 | 55 |
| 2 | niche2 | Plasma | 0.04214 | 0.0262 | 55 |


因此可把 **niche2的Naive/Memory B与Plasma相关环境** 作为APP–CD74路线的对接入口；niche1虽有GC B相关组成，但本轮scRNA GC B覆盖不足，不宜作为已获支持的首选GC轴。这里的B/Plasma富集来自既有组成定义，不是新的独立发现；需要下一轮确认标签版本、切片/患者单位及映射。

建议下一轮具体检验：在控制B/Plasma总量、切片与患者后，MP10高/低上皮与这些既有niche是邻接、重叠还是回避；APP高的epi是否位于CD74阳性B/Plasma邻近边界；反方向候选的ADGRE5阳性B与CD55阳性epi、CD46阳性Plasma与JAG1阳性epi是否在相应位置匹配。HRA应预设为反向待解释队列，GSE不能当共同覆盖稳健重复；IFN-B单独作为Inhouse探索层，不与公共Naive/Memory标签直接合并。

## 6. 交付与阅读顺序

先读本报告及`methods.md`，再看`FIGURE_INDEX.md`。每个dataset都包含主状态和原条件图、B/Plasma覆盖/比例图、真实LR圈图和配对表达图，含负向、无支持及NE，不只收录ALL或阳性图。透明PDF与300dpi PNG留在完整目录；轻量`review.zip`只装核心PDF、关键表、方法和报告，不含原始表达、母对象或大型输入。完整数值网络和重现脚本留在v2目录，v1未改写。

本轮结论止于探索性关联及空间可检验假设，没有患者层面CellChat差异P、空间接触证据或功能因果结论；不自动开始下一轮空间模块。
