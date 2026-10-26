# MP10 第一轮亚群筛选报告

已完成模块 0–3，停止在候选汇总；未运行新的空间、CNV、niche、CM 或 COMMOT/NICHES 分析。

## 核心结论

跨队列主要比例分析未发现经过 BH 校正的 MP10 关联大群或亚群。不能据此宣称 MP10-high 存在已证实的某类细胞扩增。通讯与样本表达核查支持两个优先探索候选：**Treg、FOLR2 TAM**；**C1QC TAM** 为证据较弱的次级候选。所有通讯差异均为 pooled 模型方向与样本表达支持，不是患者层面差异显著性或因果关系。

## 实际输入与覆盖

保留 374,329 个有效细胞、58 个标本、5 个 dataset。患者映射共 42 个记录 ID，均为 derived/provisional，不能冒称已临床核实的 unique patients；当前 verified 映射数为 0。InhouseData 完整保留 13 标本、12 个记录患者 ID。原始注释映射中已剔除的污染细胞本轮不恢复；为训练空间参考而标记不纳入的其他标签没有被丢弃。

| dataset | specimens | patient_IDs_recorded | normal | tumor | MP10_eligible |
| --- | --- | --- | --- | --- | --- |
| CNP0000460 | 12 | 7 | 5 | 7 | 5 |
| GSE326225 | 8 | 8 | 0 | 8 | 6 |
| HRA003620 | 14 | 7 | 7 | 7 | 5 |
| InhouseData | 13 | 12 | 0 | 13 | 9 |
| PRJNA662018 | 11 | 8 | 3 | 8 | 2 |

实际对象、字段与绝对路径见 [input_sources.tsv](tables/input_sources.tsv)。核心表达对象为 `final_annotation_v2_pre_reference.qs` 的 RNA:data；比例注释为 frozen v2 `major_celltype/subtype_v2`。使用既有 `malignant_like_main` 判定与 Top20_score，按标本目标上皮均值汇总。冻结 Top20 全部一致，未重选基因、改评分或借用其他工作结果。

主要 MP10 覆盖为 27 个肿瘤标本、26 个记录患者 ID；Inhouse 其中 9 个标本。低覆盖均值仍保存，但不参与阈值或主要 MP10 模型。各 dataset 内中位数定义相对 high/low，不是统一绝对代谢阈值。

## 哪些比较可估计

- Normal/Tumor：58 个标本，15 Normal、43 Tumor。CNP0000460、HRA003620、PRJNA662018 可作队列内对照；GSE326225 与 Inhouse 只有 Tumor，单组描述保留。另有按已有映射的配对结果，解释以映射正确为前提。
- 精确 T1/T2：17 个 Tumor 标本。只有 CNP0000460 同时覆盖两组，跨队列调整模型的主要组间信息来自此队列；HRA003620 仅 T1，Inhouse 仅 T2，不能单独作 T1/T2 对照。GSE326225、PRJNA662018 缺少可靠 T 分类，没有将 MIBC/NMIBC 换算为 T 分期。
- MP10 二分及连续：27 标本；部分 parent 缺失使有效数为 26。PRJNA662018 只有 2 个合格标本，不能拟合其单队列比例模型或当作有重复样本的通讯验证。
- 分期完整的连续 MP10：19 标本，Fibroblast parent 有效 18。未调整与分期调整严格使用相同子集。没有明确 pT/cT 信息，记录为 T 分类但来源前缀未注明。
- 恒定比例、缺父群、单组、设计不满秩或残差自由度不足均明确保留为不可估计；完整原因见 [model_status_summary.tsv](tables/model_status_summary.tsv)。没有通过删掉必要 dataset/分期项制造显著性。

## 大群先解释，再看亚群

Normal/Tumor 存在清楚的组成差异。例如 Fibroblast 在 TME 中的份额下降，同时其内部 Pericyte、Matrix_CAF 等相对份额上升，iCAF 相对份额下降。这是组织状态背景，不直接等于 MP10 证据。下表原始差值是未调整的样本均值百分点差；模型系数是 natural-logit 尺度，二者不能混称。

| subtype | denominator_type | coefficient | raw_difference_pp | P | q |
| --- | --- | --- | --- | --- | --- |
| Pericyte | parent | 3.564 | 22.02 | 4.361e-13 | 1.875e-11 |
| Cycling TAM/macrophage | parent | 2.274 | 3.898 | 9.037e-10 | 1.943e-08 |
| vCAF | parent | 2.498 | 4.175 | 6.516e-09 | 9.34e-08 |
| Matrix_CAF | parent | 2.507 | 15.72 | 3.56e-07 | 3.827e-06 |
| rCAF | parent | 2.121 | 0.5111 | 1.351e-06 | 1.162e-05 |
| Capillary EC | parent | 1.897 | 22.11 | 2.015e-06 | 1.444e-05 |
| dataset_bias_fibro | parent | 2.103 | 3.269 | 2.403e-06 | 1.476e-05 |
| Fibroblast | TME | -2.438 | -6.857 | 7.386e-06 | 5.17e-05 |
| iCAF | parent | -1.979 | -45.48 | 1.084e-05 | 5.827e-05 |
| Venous EC | parent | -1.389 | -34.41 | 2.175e-05 | 0.0001039 |
| Classical Mono | parent | -1.04 | -31.11 | 3.704e-05 | 0.0001593 |
| Cycling stromal | parent | 2.064 | 2.068 | 5.123e-05 | 0.0002002 |

T1/T2、MP10-high/low、连续 MP10 及其分期调整在主要口径均无 q<0.05 结果。连续 MP10 的较前列结果如下（完整表包含全部受检群体）：

| subtype | denominator_type | coefficient | CI_low | CI_high | P | q | n_specimens |
| --- | --- | --- | --- | --- | --- | --- | --- |
| iCAF-like | parent | -4.757 | -9.294 | -0.2203 | 0.04025 | 0.7985 | 26 |
| Pericyte | parent | 4.803 | -0.1107 | 9.717 | 0.05518 | 0.7985 | 26 |
| dataset_bias_myeloid | parent | 3.54 | -0.13 | 7.211 | 0.05822 | 0.7985 | 27 |
| Plasma | TME | -6.723 | -13.79 | 0.3474 | 0.06181 | 0.4327 | 27 |
| Naive B | parent | 2.772 | -0.8915 | 6.435 | 0.1316 | 0.7985 | 26 |
| excluded_myeloid | parent | -2.837 | -6.996 | 1.322 | 0.1747 | 0.7985 | 27 |
| Endothelial | TME | 4.511 | -2.322 | 11.34 | 0.19 | 0.665 | 27 |
| IFN-stimulated B | parent | 4.715 | -2.644 | 12.07 | 0.1987 | 0.7985 | 26 |
| rCAF | parent | 2.933 | -1.819 | 7.685 | 0.2207 | 0.7985 | 26 |
| NK | parent | 2.235 | -1.488 | 5.958 | 0.2315 | 0.7985 | 27 |

iCAF-like 的连续方向为负，Pericyte 为正，但均未通过多重校正；Pericyte 在分期调整后支持进一步减弱。Plasma 的二分与连续方向均偏负，也不显著。每个记录患者仅保留一个标本的敏感性分析仍无显著 MP10 比例结果。重复肿瘤患者很少时，相关结构借用同一分母全图谱的估计，方法及假设已单列，不能掩盖映射尚未核实的问题。

最终候选的比例证据同样为非显著，不能把通讯候选重新包装成比例候选。Treg 的二分系数略负、连续系数略正，不能称为方向一致；FOLR2 TAM、C1QC TAM 的二分与连续系数均正，但均不显著，分期调整后也没有转为显著：

| subtype | comparison | coefficient | P | q | n_specimens |
| --- | --- | --- | --- | --- | --- |
| C1QC TAM | MP10_low_high | 0.1107 | 0.6445 | 0.9029 | 27 |
| FOLR2 TAM | MP10_low_high | 0.4135 | 0.1424 | 0.7465 | 27 |
| C1QC TAM | MP10_continuous | 0.8012 | 0.6069 | 0.7985 | 27 |
| FOLR2 TAM | MP10_continuous | 1.627 | 0.3842 | 0.7985 | 27 |
| C1QC TAM | MP10_stage_complete_unadjusted | 1.086 | 0.578 | 0.8936 | 19 |
| FOLR2 TAM | MP10_stage_complete_unadjusted | 1.277 | 0.5525 | 0.8936 | 19 |
| C1QC TAM | MP10_stage_adjusted | 0.641 | 0.6868 | 0.9973 | 19 |
| FOLR2 TAM | MP10_stage_adjusted | 1.005 | 0.6163 | 0.9973 | 19 |
| Treg | MP10_low_high | -0.1152 | 0.6791 | 0.9029 | 27 |
| Treg | MP10_continuous | 0.4798 | 0.794 | 0.8576 | 27 |
| Treg | MP10_stage_complete_unadjusted | 0.652 | 0.7525 | 0.9517 | 19 |
| Treg | MP10_stage_adjusted | 0.7109 | 0.7502 | 0.9973 | 19 |

## CellChat 候选与样本核查

完成 28 个可用条件网络，包含总体与全部队列；使用 CellChat 1.6.1 human DB、triMean、population.size=FALSE、100 次内部置换、每群≥30 细胞，并在 sample×subtype 内封顶 100 细胞。大群概览是细分网络的聚合。原始样本表达支持采用未抽样的 RNA:data 平均表达和检出率，不把 pooled 内部 P 当患者差异 P，不比较独立单样本对象的未校准 probability。

以下“有重复支持”指低、高组各至少 2 个记录患者 ID 有≥30 个 source/target 细胞，且配体及所有受体亚基检出率均≥10%；这只是表达覆盖标准，不是统计显著性门槛。表达乘积及删一患者检查只是描述性支持，不是新的 CellChat probability。

### 优先探索 1：Treg

**Malignant_like_epithelial → Treg → NECTIN2–TIGIT（Cell–Cell Contact）**。CNP0000460、GSE326225、InhouseData 的网络均为 high>low，样本表达支持同向；这三个队列删去任一个记录患者 ID 后，表达支持方向均保留。总体 LR 支持强度为 low 0.008875、high 0.031679（共享图形标尺，不能当患者效应量）。HRA003620 high 只有一个可评估标本，表达方向略反向，不计为重复验证；PRJNA662018 缺少双方有效组覆盖。

目前最符合“通讯状态变化线索”，不是已证实的 Treg 比例变化。T1/T2 中该轴方向与 MP10 对照不同，不能称为所有对照一致增强。下一轮空间问题：MP10 较高的恶性样上皮是否更常直接邻接 TIGIT 表达的 Treg；这种接触富集能否在控制 Treg 数量和组织区域后保留。

### 优先探索 2：FOLR2 TAM

**Malignant_like_epithelial → FOLR2 TAM → MIF–(CD74+CD44)（Secreted Signaling）**。CNP0000460、HRA003620、InhouseData 均有重复样本表达支持，网络及表达支持均 high>low；删一患者后同方向比例分别为 80%、100%、87.5%。总体 low 0.099946、high 0.152394。GSE326225 low 仅一个可评估标本，表达支持方向反向，不计为重复验证。MIF 对多个髓系亚群有支持，不能称为 FOLR2 TAM 特异通路。

仍是通讯状态线索：FOLR2 TAM 的主要比例模型不显著。下一轮空间问题：MP10 较高、MIF 表达的恶性样上皮与 CD74/CD44 支持的 FOLR2 TAM 是否形成可重复邻近，并在匹配其他 TAM 密度后仍保留该关系。

### 次级探索：C1QC TAM

**C1QC TAM → Malignant_like_epithelial → HBEGF–(EGFR+ERBB2)（Secreted Signaling）**。CNP0000460、HRA003620、InhouseData 均有双方及复合物表达覆盖，网络均 high>low；但 HRA003620 的样本表达支持方向为负且删一后仍为负。CNP 与 Inhouse 表达支持为正。总体 low 0.001807、high 0.002626，变化较小。这种明确不一致使其优先级低于前两者。

下一轮空间问题：HBEGF 表达的 C1QC TAM 是否紧邻 EGFR/ERBB2 支持的 MP10 上皮区域；空间分布能否解释 HRA003620 中表达与 pooled 推断不一致。

每个候选的完整 Normal/Tumor、T1/T2、MP10 二分/连续/分期调整、parent/TME 分母、队列与有效样本证据见 [candidate_evidence.tsv](tables/candidate_evidence.tsv)、[candidate_axis_dataset_support.tsv](tables/candidate_axis_dataset_support.tsv) 和逐亚基 [selected_LR_gene_sample_support.tsv](tables/selected_LR_gene_sample_support.tsv)。

## 未优先推荐及反向结果

CAF、SPP1 TAM、Neutrophil、Endothelial 和全部其他已注释亚群都进入了统计和可评估网络。没有因大群 P 值过滤小群。CAF 的组织状态差异没有转化为本轮充分的 MP10 比例或多队列样本通讯支持；SPP1 TAM 部分轴的样本表达与网络方向不一致；Neutrophil 的 MIF–CD74/CXCR4 负向支持主要来自 Inhouse；Capillary EC→上皮 COL4A1–SDC1 为负向探索，但重复样本表达覆盖主要在 CNP0000460。它们并未被删掉或补为零，详见完整候选表与各队列图。WNT5A–MCAM 未被预指定为结果。

## 图形入口与交付

- 图 A：`figures/overall/stacked_MP10_low_high_Fibroblast.pdf`、`stacked_MP10_low_high_Myeloid.pdf`、`stacked_MP10_low_high_T_NK.pdf`；其余组织状态、T1/T2、大群和每个 parent 均在同目录及 `figures/by_dataset/<dataset>/`。总体柱使用共同 dataset 等权，绝非 pooled 细胞比例。
- 图 B：`figures/overall/circle_MP10_low_high_Treg_NECTIN2_TIGIT.pdf`、`circle_MP10_low_high_FOLR2_TAM_MIF_CD74_CD44.pdf`、`circle_MP10_low_high_C1QC_TAM_HBEGF_EGFR_ERBB2.pdf`。两条件使用相同节点与边宽最大值；[NE] 明确表示不可评估。正常组没有恶性样上皮，相关面板保持不可评估，不能解释为 Normal→Tumor 同类细胞通讯从零增大。
- `effects_*.pdf`、`candidate_points_*.pdf`、`candidate_continuous_*.pdf`、`network_major_*.pdf`、`pathway_bubbles_*.pdf`、`sample_support_*.pdf` 为配套结果；透明 PDF 均另有 300 dpi PNG。
- `review/` 为轻量阅读包；`review.zip` 不含原始数据、完整表达矩阵、Seurat 或完整 CellChat 输入对象。完整脚本、参数、表格与必要数值缓存留在本目录。

## 局限与完成边界

患者映射、部分正常来源、治疗状态与 pT/cT 前缀未得到临床核实。沿用的是既有 malignant-like 主集合而非独立恶性真值；低恶性上皮覆盖会限制 MP10 样本代表性。Inhouse 部分亚群带既有基因命名/技术偏差标签，已原样保留；未整体更换注释或补换基因。CellChat 各对象存在内部最大表达缩放，pooled 网络差值仅探索，候选的样本表达支持不能替代独立患者差异检验。两个主要候选值得空间验证，但本轮没有已确证的 MP10 亚群机制。

本轮停止于上述候选和空间问题，不自动进入下一轮，也未发布或上传到外部仓库。
