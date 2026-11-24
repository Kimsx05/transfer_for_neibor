# 方法与复现

## 项目、身份与计数

实际根目录为 `/home/data/t070721/codex_workspace/Bladder Metabolism`。指定的 `/home/jinshengxi/bladder_metabolism` 不存在。用 rg 检索当前工作区及隐藏文件，定位到 v1/v2 报告、原表和脚本，未找到项目 AGENTS.md 或指定 JSON。结果目录启动时尚不存在。本次只运行 B。

读取原 `mp10_cellchat_high_definition_sensitivity_v1/scripts/01_prepare.R`、机制 v1 的 pseudobulk、模型与患者映射代码，以及跟进 v1/v2 的评分脚本。复用 MP10 Top20（成员见 registration.md，保留 ACAT2）和 MP10_cell（原全 RNA 90,011 基因的平均并列秩）。标本内未舍入的 R type-7 分位数确定归属。High 主为 >Q95，Low 为 ≤Q50。H10/H50 分别为 >Q90、>Q50，限于既有敏感性；边界并列不拆分。未将固定比例作为跨标本暴露。

主身份为原 c2l 最终 Treg。逐标本保留原 subtype、state 组成和 CD3D/E、CD4、FOXP3、IL2RA、CTLA4、IKZF2 等标记的 UMI、阳性细胞数。原 pseudobulk.qs 同时包含全 RNA 的整数计数、逐基因阳性细胞数与 profile 信息。核实每份 Treg profile 的总 UMI；受者每标本一份。载体来源统计只计 Epi_all 与其余互斥细胞类型，按标本先计算 RNA 份额，再跨主患者代表取中位数。

由 specimen_exposure_eligibility、patient_relationship_audit、sample_metadata 三表核对 specimen、来源、原始 ID、记录 ID 与重复关系。另一份 resistance v1 映射中 39/40 份标本有同一记录 ID，T4450398 在该辅助表缺失。两来源沿用部分相同的命名规则，均不能替代临床身份验证。T4463209a/b 为已确认重复，主取 a，预设替换取 b，每个模型每个记录单位最多一行。其余保留 derived/provisional。没有引入 LLPS 的患者分母。

## 公开扰动与基因注释

缓存为 `cache/mp10_public_resources_v1/<accession>/snapshot_20261124/`。先写临时文件，核验 gzip、TAR、SOFT 或表头后原子改名。resource_manifest.tsv 记录 URL、HTTP 状态（000 表示没有获得 HTTP 响应）、字节、缓存 SHA256、实际列或成员。GPL 直接获取官方平台 SOFT 后 gzip 压缩，清单校验值对应压缩后的文件。GSM 与原始列、元数据及排除原因见 tables/sample_mapping.tsv。最初的 GPL23038 family 入口下载范围过大，中止于未取得最终 HTTP 状态，改用 platform self；该残缺文件未参与分析。

| 比较 | 实际样本 | 处理 |
|---|---|---|
| SCAP，B16 Treg | GSM4505100–03 WT，4505104–06 KO | 匹配 WT 4、KO 3，log2 RMA |
| FASN，B16 Treg | GSM4505107–09 WT，4505110–12 KO | 另一批 WT 3、KO 3，log2 RMA |
| 一般活化，外周淋巴结 WT Treg | GSM5029123–27 为 0 h，5029128–32 为 8 h | 各 5，log2 RMA，非配对 |
| CD36，Yumm1.7 肿瘤 Treg | GSM4138181–83 WT，4138184–86 KO | 各 3，逐样本 TXT，FPKM |
| SREBF2，VAT Treg | GSM7777514–16 WT，7777517–19 KO | 各 3，仅 VAT，RPKM |

GSE149573 核实 28,846 个转录簇、25 列；GSE165257 为 41,345 行、10 列。两平台分别建模，不再次 log。GSM4505104 的标题、B16 source、Foxp3-YFP⁺、tumor infiltrating Treg 和 Scap KO 等结构化字段，与原论文 B16 对照/KO 的 4/3 设计一致。OT-I 自由描述在该系列全部 25 份样本中相同。保留原字段并按结构化信息纳入。RNA 样本数按 GSM 计，rep 序号不作为供体配对；各 GSM 到动物／混池的完整关系未知。[原论文](https://www.nature.com/articles/s41586-021-03235-6) · [GSM4505104](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4505104)

GPL23038 精确连接转录簇 ID，从 mrna_assignment 的 RefSeq/UCSC 描述提取 RNA／转录本说明前括号内的 symbol。GPL16570 精确连接数字探针 ID，解析 gene_assignment 各命中的 symbol 字段。仅保留能解析为唯一 symbol 的特征，未解析或歧义者另列。两平台均在差异分析前按 symbol 对同一样本的多个探针取中位数。CD36 按 TAR 成员的 GSM 与 gene_id 合并，核对六份 gene_name。VAT 用完整 GEO title 对应相同的组织、基因型与 rep 列名。两个丰度矩阵均使用 log2(x+1)，同名特征按样本取中位数。常量或非有限行排除；VAT 未列入检测表的基因保留缺失。脾脏 GEO 为 2 WT/3 KO，而文件列为 3 WT/2 KO，整支暂停。[CD36 数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139325) · [VAT 数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243021)

limma 拟合截距和组别，芯片使用标准 empirical Bayes，RNA 丰度使用强度趋势。全表保存效应、95% CI、moderated t、P 和每个对比全基因 BH q。KO_UP/DOWN 对应 KO−WT；WT-dependent 取负 t，普通活化取 8 h−0 h 正 t。严格集合 q<0.05。SCAP 涉及 SREBP1/2、多支脂质及甲羟戊酸后果。CD36 的脂质适应没有胆固醇特异性，VAT SREBF2 用于跨组织参照。GSE149573/165257 来自同一研究。

先做小鼠平台内差异，再用冻结的 MGI 人鼠同源类（2026-11-24，快照 SHA256 见清单）。同类每个物种各一个 symbol–Entrez 配对，且两端 symbol 各自只出现在一个类中。记录多对多、缺物种、缺精确 symbol 与常数排除。没有跨平台表达合并。[MGI 同源表](https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt)

## 排名、标本模型与异质性

新关联前登记见 registration.md。旧摄取 ssGSEA 的 6,998 基因背景来自 Treg ≥30 参考 profile 的单组 edgeR filterByExpr。原 M0046 DE 的 13,407 基因按实际设计筛选，沿用 sign(logFC)×sqrt(max(F,0))。新评分从全基因 Treg counts 计算 edgeR TMM，固定在主队列中至少四分之一标本 CPM≥1 的 22,224 基因作为背景。使用 TMM-logCPM（prior.count=2）在每份标本全背景内计算平均并列秩百分位 (rank−0.5)/N。每种外部扰动的有向 t 在可用的一对一同源共有基因中转成中心化百分位权重 w。投影为 Σ(w×标本秩百分位)/Σ|w|，以 ≥50、main-a 的 32 份标本固定均值和 SD 标准化。分数越高表示越接近 WT 依赖或 8 h 方向。弱扰动参照的估计误差不会因此消失。正式摄取、合成和 G2M 使用背景内集合平均百分位，门槛为 ≥5 且 ≥50% 成员。原 ssGSEA 和新分数的绝对值不跨方法比较。

在既有受者全基因排序上计算双排序 Spearman。严格正／负集合另给其相对其余共有基因的竞争性秩位置 Mann–Whitney AUC−0.5。空集保留 NA；单成员只解释成单基因位置。基因不是患者重复。标本推断使用 标准化程序 ~ dataset + 原 Epi Q50_z，常规 OLS 残差 t 区间及双侧 P。主单位为原 Q50 一个固定 SD（0.078525），Q95 SD 为 0.147472。M0047–49 为 ≥30、Q95、换 b；M0050–54 为逐一删除数据集。另给各数据集内斜率。至少 4 个记录单位、残差 df≥3 且满秩。加入 TCR 或 G2M 分别比较条件关联。

每个设置、掩码及调整分开对 SCAP/FASN/TCR 三项、CD36/SREBF2 两项 BH；四个正式过程另属次要家族。Full、去 Top20、再去 SCAP/SREBF1/SREBF2/FASN/CD36 均重新计算背景和目标。另行实施事后非重叠诊断，以严格 TCR 成员自身背景的有向投影，或 G2M 成员平均 logCPM 作控制，代谢投影背景及成员移除对应控制全部成员、Top20 和直接节点。每控制对四项 BH。此诊断不替换登记主分析。

H05 发送—接收另按原 ≥50 规则，将既有 H05 发送端 ssGSEA 与同一标本的 Treg 摄取读出连接。沿用 HC3 和残差 t 区间，单位杠杆触发主模型 NE，原去 singleton 敏感性两路线一组 BH。公开扰动、空间邻近、载体 RNA 份额各自评价其实际环节。

## 运行和文件

从实际项目根执行 `bash results/mp10_treg_metabolic_mechanisms_v1/code/run.sh` 可在日期子目录重现分析。R --vanilla，BLAS/OpenMP、data.table、RcppParallel 单线程，pseudobulk 读到 R 进程结束全程持有共享 heavy_object.lock。没有加载或复制全量 Seurat/H5AD。使用已有环境，无新增依赖。R、包及绘图版本见 logs。jins-routine-plotting 的公共 Arial 主题和 jin_export 用于全部图；该技能无对应的森林、散点、热图及证据表模板，因此使用定制几何图层。每图含 PDF、300 dpi PNG 及 TSV，带图例者另存图例与主体。图注见 FIGURE_INDEX.tsv。
