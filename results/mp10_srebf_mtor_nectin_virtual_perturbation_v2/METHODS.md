# MP10–SREBF–mTOR–NECTIN 虚拟扰动 v2

## 输入和配对设计

直接使用本机 metabolism 项目快照及只读 v1 结果。核对 v1 的 METHODS、RESULTS、建模/准备/面板/汇总代码、scTenifold 运行配置及 manifest。旧 design.json 中的传播或 ridge 设想未用于描述 scTenifold。唯一主配置是 `config/main.json`。本轮只对五个主源分别置零，另有固定匹配对照、旧匹配对照基因。所有结果来自同一观察性数据。公开项目可用于核对工作流。[参考仓库](https://github.com/Kimsx05/transfer_for_neibor)

沿用 v1 主条码、原序、MP10 分数、标本内 H05/H10/H50/L50 标签、来源及上皮簇。8,000 个细胞覆盖 38 位记录患者、39 个标本和 5 个来源，T4463209a/b 共用一个患者键。原始 pre-harmonization H5AD 共 90,011 个 feature。对这 8,000 个细胞逐行核查非零支持、非负整数计数、检出基因数和完整总 UMI，逐值与 v1 稀疏缓存比对，总计 80,436,985 UMI。未重新抽主样本或计算 MP10。高组彼此嵌套。

HGNC 2026-11-18 快照用于标准化。优先匹配现行符号，继而核对 Ensembl ID、旧名和别名。多候选不猜测，用 RAW~ 前缀保留标识及候选记录；同一标准基因的原始行求和一次，重新验证每个细胞总 UMI。gene_mapping 保存逐原始 feature 支持，coverage 保存标准名层支持。至少 10 个细胞、3 位患者检出且原始计数方差为正，子集重拟合重新检查。原资格、入模、有效出边、拟合内有向路径及模拟合法性分别可查。

按 HGNC 细胞质核糖体大小亚基蛋白家族、线粒体核糖体蛋白家族（含 MRPL、MRPS 对应结构成员）及线粒体染色体精确去除，显式计入小亚基蛋白 RACK1。保留核编码 OXPHOS 及 RPS6KB/KA 激酶。完整版本化名单含未入模基因，in_P_full 指示本轮实际删去的行。核糖体与线粒体表达可能有生物学意义。[HGNC](https://www.genenames.org/) [UniProt RACK1](https://www.uniprot.org/uniprotkb/P63244/entry)

从 v1 面板、合格主源、预设读出、冻结程序、旧对照出发。对每个主源在原始未扰动表达中按 log1p 均值与检出率的标准化距离补充至多 15 个对照候选。两个 SREBF 的候选为 Lambert v1.01 人类序列特异 TF 清单；三个非 TF 的候选为 HGNC 蛋白编码基因，排除该清单及较宽的旧 TF 目录。主假说节点、预设程序、Ribo/Mito 不参与匹配。本轮在读取 KO 结果前将初步 motif 注释候选换成该有定义的 TF 目录，重新拟合两张 WT 并冻结名单。P_full 包含 2,233 个基因，P_clean 从中准确剔除 39 个，留下 2,194 个，无 HVG 回填。[Human Transcription Factors](https://humantfs.ccbr.utoronto.ca/download.php)

在两张 WT 均有出边的候选中，以 clean 主种子 WT 的 log1p 原始均值、检出率和 log1p 出度各自标准化，取欧氏距离最近的 5 个，按相同名单核查 P_full 距离和各变量差异。候选空间和匹配质量见 control_matching。旧名单独立标识，ANXA1、H2AFZ 等按旧匹配对照解释。五对照比较为描述性，无经验 P 值。

## scTenifoldKnk

固定 CRAN 源码 scTenifoldKnk 2.0.0、scTenifoldNet 1.4.3。使用 v1 的本地版本化源码包及只读 R 库。各模型用全基因组原始总 UMI 为分母：稀疏面板 counts 下另加一行面板外 UMI 总数并核对列和，随后只选真实模型基因做 PCR 或相关计算。该载体行不进网络/排名。v1 实际建网和方向计算使用面板分母。v1 到 v2 尚有映射、面板及有向门控差异，纯过滤效应看 v2 的配对比较。未额外拟合旧分母模型。

每张 WT 使用 10 个 PCR 子网，每次 200 个细胞，3 个 PC，边保留分位 q=0.9。张量 CP 秩 3、最大 1,000 次、阈值 1e−5、边保留三位小数。strictDirection lambda=0，自环归零，转置使行表示出边。heatT=10 的官方 hkManifoldAlignment 距离与 t=5 的 WT 共表达方向分别计算。方向主口径为原生未校正，另列完整文库大小回归及深度回归后移除 PC1 的诊断。用标准化向量 Gram 矩阵计算同一相关核，主种子子块与官方 knockoutDirection 结果的最大差为 2.22e−16。heat 核使用对称化；有向可达性另列。无有向路径的原生统计仍保存，生物学路径解释见 endpoint_path_coverage 与 program_path_sensitivity。

原生 .drStatistics 的 FC 为距离平方除以非 KO 行距离平方均值，P 为原包理论距离检验，噪声地板沿用包约定，BH 校正总体包含源自身。Z、FC、P、q 均按原生数值保存。基因按非源行距离降序排名，同距取均秩；分位为 (rank−1)/(N−1)，0 代表最靠前。方向评分不代表表达变化百分比。主 clean WT 另对五个主源和每源第一匹配对照删除该源出边行，官方 manifoldAlignment（d=2）及 dRegulation 对全基因计算。热核与经典分支的原始距离单位不作横向比较，两者共享的方向不重复计证。[scTenifoldKnk 官方代码](https://github.com/cailab-tamu/scTenifoldKnk)

全样本种子 123、456、789。20 次患者级无放回抽取 ceil(0.8×38)=31 位患者，抽样种子为 12000 加重复编号，每次从 counts 重新建 WT。对五个来源逐一留出、对至少 500 个细胞且 4 位患者的原簇各自重建 WT。每张 WT 只拟合一次，供所有合格主源及对照共享。稳定性范围取观察到的最小和最大值。患者及 High/Low 条件读出固定 clean WT，仅更换各层表达均值及离散度，至少 10 个细胞且源与终点变异可估，汇总以患者等权，同一比较只用共同合格患者。

## CellOracle

分析使用独立的 Python 3.10.20 venv，CellOracle 0.18.0，实际依赖见 config/oracle_requirements.txt。复用只读项目 Python 包，修复 velocyto 编译、gimmemotifs、genomepy、ipywidgets 依赖；正式计算不依赖 Jupyter 界面。官方人类 hg38 启动子先验为 20210630，gimmemotifs v5，fpr2，threshold 10。parquet 文件与官方 Git blob 5fb1b69eef4831fc8a9270b058b46923c5e07988 一致，SHA-256 及大小见输入指纹。39,315 条启动子记录、1,094 个原始调控因子列。双轴规范化、同靶启动子取并集、去自环并限制进入分析的基因后共有 17,619 条候选边，涉及 59 个源与 1,424 个靶。该先验源于参考启动子 motif，本队列缺少配对 ATAC。[官方先验](https://github.com/morris-lab/CellOracle/tree/master/celloracle/data/promoter_base_GRN) [官方教程](https://morris-lab.github.io/CellOracle.documentation/) [CellOracle 原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9946838/)

400 个 C0 细胞小试通过后，按原簇独立做预处理、邻居插补和 GRN 拟合。各簇从同批条码取规范名 clean 特征，再施行当地表达资格。簇至少 80 个细胞、4 位患者；每例的 PCA、插补及拟合均在该簇内，按来源留一时从头重建。将每个细胞按完整文库每万归一并 log1p，保留 raw_count；通过 import_anndata_as_normalized_count，30 维 PCA，k=20 普通最近邻插补，单线程。get_links 用 bagging_ridge（20 次 bagging，alpha=10），filter_links 按 P≤0.001、绝对系数排序每簇最多 2,000 条，实际边界条件以安装版本 filter_links 为准；get_cluster_specific_TFdict_from_Links 之后以 alpha=10 和 use_cluster_specific_TFdict=True 拟合模拟系数。逐簇断言非零系数边均在筛选边及先验中。该图是预测调控图。

正式主种子 123，另一种子 456 改变输入细胞排列及采样对应条码，包内使用其固有随机状态。另有 5 组来源逐一留出。对 SREBF1/2 与预先匹配、具有有效出边且符合官方检查的 TF 置零，主传播 3 步，主种子对两源另做 2/4 步及官方 randomized-GRN（123），并用空扰动校验数值归零。未绕过合法性检查。模拟采用原包 clip_delta_X=False，记录超出插补 WT 值域数、零边界数及官方警告。原簇过小、表达不足、缺先验源、无拟合出边、给定步数无路径及官方检查未通过均单独留痕。随机化只作为有限诊断；mTOR 三个非 TF 节点在该方法标为 NE。

imputed_count 是插补后的模型输入；simulated_count 是置零传播后的同尺度模型输出；delta_X 为两者之差，并逐值断言相等。主模型的两源关键成员另存三层数值和条码。可达簇内先求每患者每基因的 signed Δ、|Δ| 和除以簇内 WT 标准差的 signed Δ；跨簇按实际可评估细胞数汇集为每位患者的均值，再取患者等权中位数。只有源自身或在指定步数内可达的行进入患者汇总，患者/细胞覆盖随基因而变。原生不可达零值在簇级表留存。非源基因以患者等权的绝对变化降序排列。CellOracle 无原生终点 P/q；这些列为空。未用个别细胞作独立重复，未计算置信区间。

## 程序和解释

冻结 MP10 Top20（含 ACAT2）、CHOL_SYNTH、既有 SCENIC/CollecTRI/DoRothEA SREBF 靶集合。另加 MSigDB 2026.1.Hs 的 HALLMARK_MTORC1_SIGNALING。Full 使用原集合；M20 去掉与 MP10 Top20 重叠的成员，MP10 自身 M20 为空。保留原名、映射名、原始/入模/可评估数和逐成员响应。程序汇总至少需 3 个成员并覆盖原集合一半。分位中位数及方向比例为描述性读出，成员相关且部分预先保留，没有计算程序显著性。五主源比较限于共同非源基因。跨方法同时给各自原单位、原模型内排名和共同可评估集内排名。细胞状态、基因可达性及覆盖可改变比较范围。

RPTOR、MTOR、RICTOR 分别代表 C1 组分、共有催化核心及 C2 比较；转录结果仅用于评估相容性。mTOR 磷酸化、SREBP 切割及核定位、表面受体蛋白和 T 细胞功能均未测量。[Düvel 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC2946786/) [Peterson 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC3336367/) [Li 等](https://pmc.ncbi.nlm.nih.gov/articles/PMC4615269/)

## 运行及文件

在本目录执行 `./run_all`；只重绘执行 `./replot`。已保留本地输入依赖和版本化先验，预处理指纹变动或损坏会阻止复用；模型检查点核对代码、输入及输出哈希。每个阶段使用同一 systemd user slice，10 个可用逻辑 CPU 动态选取，内存硬限 200,000,000,000 字节。软阈值 180 GB 暂停启动新重任务，190 GB 停止本任务 scope；并行池内每个作业的 BLAS、OpenMP、R 计算线程均为 1。scTenifold 的上限为 6 个工作进程，CellOracle 为 8，两个重阶段顺序执行。资源记录包括作业时间、cgroup CPU、内存峰值和磁盘采样；早期环境安装及轻量交付编辑不在建模资源计量窗口内。所有大型输入仍在原位置。R 作图复用项目 jins-routine-plotting 的公共主题与导出器，排名、方向和重拟合图使用自定义几何图层，每图有 PDF/SVG/PNG 及数据。详见 FIGURE_INDEX.tsv、manifest.json 与 qc/resources.json。


### 运行中修复与审计

drop_InhouseData 的 C1 有 998 个合格输入基因，CellOracle 0.18.0 在按自身规则筛选后只有 995 个 CV 分数，却尝试读取索引 997。仅在隔离包中将该排序统计的索引封顶于最后一项；原表达、变异筛选、TF 资格、外推合法性与拟合步骤均保留。原本合法的索引得到相同分数，补丁和前后源码哈希见 qc/oracle_qc_compat.json。首次错误和重试日志保留。部分小簇收到包对基因数或 TF 数偏少的提示。

scTenifold 启动时修正过 data.table 对源行逻辑筛选的语法。两张 WT 完整性、最终成功的退出状态及所有源输出哈希均经复核；qc/checkpoint_reconciliation.tsv 记录启动时旧指纹与实际完成导出的脚本指纹。资源抽样日志中途有一次间断，qc/resources.json 提供最大采样间隔；内核的作业限额和峰值计数独立于轮询。源数据和 v1 对照文件仍按原始 SHA-256 核验。

本节点的 user slice 未暴露 cpuset.cpus.effective。CPU 绑定通过每次启动的 taskset 继承，另有 cgroup 的 cpu.max 总配额。验收在同一 scope 内读取实际亲和性和 cpu.max。
