# 方法与复现

项目根为 `/home/data/t070721/codex_workspace/Bladder Metabolism`。旧根在本机不存在，已在当前工作区检索 AGENTS.md、上游脚本、报告和轻量结果。计算、输入缓存及图表只写入本交付目录，共享输入和虚拟任务目录只读。执行脚本按顺序启动单线程 R 或 Python，BLAS、OpenMP、data.table、BiocParallel 及数值表达式线程均限制为 1。未新增依赖。

## 样本与测量

原 TCGA-BLCA PanCancer Atlas cBioPortal batch-normalized RSEM 在匹配临床的原发、结尾为 -01 的标本中，每位记录患者按排序取一份。逐行 log2(RSEM+1)，映射同一基因的行取均值；冻结评分背景按原脚本只执行 PVRL2 → NECTIN2，去常量及非有限行。ABSOLUTE 纯度按原键 array 归并。逐变量 Full/M20 复算后使用原冻结分数。旧 10 配对的 program/RNA × Full/M20 × raw/purity subset/adjusted Pearson、原 P/q 和区间保存在 `tables/TCGA_original_Pearson.tsv`，与新九项定义的 q 分列。

GSE32894、GSE31684 和 GSE13507 以官方 GEO 系列及补表为单位，逐 GSM 决定纳排。GSE32894 包含瑞典研究的 308 份肿瘤，GSE31684 为膀胱切除队列，GSE13507 为韩国原发肿瘤。GSE32549 与 GSE32894 有部分技术重合风险，未另计队列。各自使用本地样本标识与官方表，跨研究没有可链接的全球患者标识，不能逐人排除重合。样本级表保留了原始组织、临床注释及纳排理由，未记录治疗时保留 unknown。

GSE31684 使用提交的 GCRMA 尺度，GSE13507 使用提交的 Illumina 归一化表达。逐探针芯片注释来自 GPL570、GPL6102、GPL6947；处理版本存入校验表。剔除多基因歧义，按明确符号合并 PVRL1–4 与 NECTIN1–4，多探针对每个样本取中位数。受体 RNA 另用 org.Hs.eg.db 和平台唯一 Entrez ID 核对。TIGIT 对应 201633，GPL6947 的 VSTM3 与 GPL6102 的 VSIG9 均有明确匹配。未据关联强弱选择探针。CD226 与 PVRIG 亦输出逐探针核对。不同平台内独立计算，单基因保留队列内提交单位。

GSE32894 的发布矩阵用于单基因，程序改用官方背景校正强度。从技术表按列名匹配正常及重复，单个 UC_0169_1 技术列板号冲突时只用其余可匹配列。提取方案逐探针将 AllPrep 中位数调至 RNeasy，P1/P2/P5/P6 中位数调至 P4。缩放只使用有限正参照及组中位数，保留可计算的探针；分位数归一化后取 log2，再在可映射且强度可取对数的探针中保留平均强度上半部。对肿瘤样本按基因取探针中位数，去常量。程序使用未中心化矩阵。预处理辅助样本从未进入患者统计。该流程是根据公开强度与技术表重建，单个冲突列以及缺少原 BASE 中间对象限制逐数值等价性。

## 程序、重叠与统计

冻结定义及各集合成员见 `config/freeze.md`、`tables/variable_dictionary.tsv` 和每个队列的 `score_coverage.tsv`。MP10 为原 Top20；SREBF1/2 为纠正过的 regulon；mTORC1 为项目冻结的 Hallmark；附加 Treg 与总体 T 细胞代理源于同项目 HD 的 11/5 成员表。各集合至少覆盖 5 个、且覆盖当前集合的 50%。使用完整可用基因背景的 GSVA ssGSEA，alpha 0.25，normalize=FALSE。M20 保留原 MP10，其他下游集合和背景去除 Top20。另在 Full 背景中从一对程序的两侧去除共享成员，或从上游去除单基因靶点。无共享时该敏感性与原值一致。RNA 层中的 SREBF1/2 为转录本，主层为程序，MTOR 的单基因值单独提供。单细胞 H05、L50 遵循原标本内 >Q95、≤Q50 成员定义及并列规则，bulk 暴露为连续值。

主分析双侧 Pearson，逐配对有效 n，Fisher z 95% CI，t 检验。外部 9×3 BH 为主 q；另给队列内 9 项 BH。家族 20×3、Treg 3×3 和补充配对各自校正，TCGA 自成参考层。敏感性按队列、清单及模型分别校正。重复展示的 NECTIN2 对保留同一 analysis_id。敏感性包括秩相关、临床分层、单个分期或研究亚型/cluster 协变量、免疫成分、TCGA 的 ABSOLUTE 纯度、M20 与去重叠。在完整 T-lineage 未达门槛时另列 CD3D 单基因补充代理。Treg 丰度代理作为结局时不再调总体 T 细胞成分。偏相关将两变量分别对含截距协变量作残差化，以剩余自由度计算近似区间和 P；秩相关的区间亦为近似。原始 r、偏相关、Spearman ρ 和回归 β 分列。原 Pearson 的删一变化用完整离差和解析计算；不按影响删除个体。有研究侵袭性标记的 GSE13507 首选该标记，并另列纯 T 分类分层，其余外部队列按 T 分类，Tx 缺失。由于平台输入及覆盖有差异，仅汇总同一 Full 定义、同一相关类型。外部 Fisher z 方差 1/(n−3)，REML 估 τ² 后给正常近似区间、Cochran Q 与 I²。3 个队列时对异质性和均值区间均应谨慎解释。

## 扰动与受者

GSE236742 仅用 HCV29Rawcounts 工作表，Ensembl ID 经 org.Hs.eg.db 唯一映射到符号，同符号计数相加。GSE143803 使用两个细胞系各自提交的 gene readCounts，按列名与 RNA GSM 标题核对，ChIP 排除。计数须为非负整数，使用 edgeR filterByExpr、TMM、limma-voom 与经验 Bayes。GSE267018 只用人 HeLa FPKM，合并同基因表达后取 log2(FPKM+1)，去常量，limma-trend 用六组设计拟合；两条件各自 KO−WT 及两个差中之差均保留。每个对比的基因区间为 moderated t 区间，q_genome 全基因 BH，q 为预定家族 5 个基因 BH。完整排序可在 genome.tsv.gz 按 t 排序。程序先去除五个家族成员，再以每份 RNA 为单位独立拟合程序分数，普通残差 t 区间，按每个对比的 7 个程序 BH。G2M、凋亡和全基因响应比例用于描述广泛状态改变；无直接活力测量。额外凋亡及胆固醇稳态集合为 MSigDB 2024.1.Hs。没有单细胞或小鼠样本参与人源干预模型。独立的 RNA 提交重复数不等于已核实的独立细胞克隆数。

配对比较限定原 main_patient、≥50 个原注释 Treg，T4463209a/b 按原规则合为一个记录患者单位。原 Treg（全部、Full、standard）OXPHOS 为唯一受者结局；逐一加入一般 TCR、G2M 或 CHOL_SYNTH，不使用 NC/CY 子集作为独立重复。去掉控制集合与 OXPHOS 重叠基因后，TCR 沿用已有严格方向签名的秩加权读数，另两者按基因平均 TMM logCPM。控制读数再标准化。Epi 暴露为原 Q50_z 或 High 内 NECTIN2 TMM logCPM，在各模型分母内以 SD 表示。模型含 dataset 固定效应；报告 HC3 区间、t 自由度、杠杆、秩和分母。High 全集的 h=1 保留描述点估计。另列每数据集有多个记录患者的预设敏感范围。细胞、基因与重采样不增加独立 n。

## 复现入口

从项目根运行本目录的 `run_all`。只重绘时运行 `replot`，只刷新虚拟综合报告时运行 `refresh_virtual`。冻结输入文件及 SHA256 位于 `config/resource_manifest.tsv`，GEO 下载记录含 HTTP 状态，逐 GSM 映射见各队列表和扰动目录。每图 PDF、300 dpi PNG 及绘图 TSV 见 `FIGURE_INDEX.tsv`。散点的线为逐图 OLS 均值和 95% 置信带；热图固定 −1 至 1。使用 jins-routine-plotting 的公共主题与 Cairo / ragg 输出，统计图采用对应数据定制。虚拟接口仅在生产者 manifest 标记 complete 后导入声明的三份表，保留原量纲。
