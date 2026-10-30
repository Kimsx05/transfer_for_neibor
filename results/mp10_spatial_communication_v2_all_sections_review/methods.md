# 方法与解释边界

本轮在只读 v1 的已执行定义、隔离环境和缓存上扩展全部 63 张普通 Visium。v1 两片输入、参数与聚合口径一致的四个方法缓存以 CACHE_REFERENCE.json 指向原路径；所有新表图在 v2。母清单 selected 仅为上一轮试点选择，不是本轮准入条件。不复制大型 RNA 对象或作大文件全量 hash，不重新训练 MP10、cell2location、CNV 或区域模型。

## 输入、组成与状态

采用 hard-QC 后完整组织内位置、非负实测 counts；按每个 spot 全部实测基因的总量归一化至 10,000，再 log1p 一次。不只在 LR 基因上算 library size，不用 Harmony/scVI latent、integrated scale、负残差或 ALRA 代替实测 RNA。

accepted cell2location v2_30type 的全部 30 原始参考类型 posterior mean 以全类型总和归一化为 Q。没有 q05 替换、dominant-label 硬分类、按展示类型重归一化或逐位置丰度截断。保留 T_NK_unresolved 的不可辨识质量，Mural 包括 Pericyte/SMC；未建模 Tip/Cycling EC 明确 NE。

冻结 EpiMP10_Top20_AUCell，Epi anchor 为 Epi posterior-mean Q≥0.25；各片合格 anchors 内严格高于/低于 MP10 中位数为 high/low，等于中位数、不可评分、非 anchor 为 other。Q_EpiH+Q_EpiL+Q_EpiOther=原 Q_Epi；TME 共享同一份 Q。没有 malignant-like/CNV 过滤。EpiH/L 是高/低关联 spot 内的上皮成分，不是在同一混合 spot 分辨出两种上皮细胞。

类型支持：至少 10 个 spot 的原始 posterior mean abundance≥0.01 且全片 abundance 总和>1。Epi 状态用分配到状态的原始 Epi mean 逐一检查。这个条件只决定可评估性，不截断 Q。H/L 任一状态不足时仍推断完整位置网络，状态为 network_only，H/L 为 NE，不画空的“有效对照”圈图。单一类型不足不删除该片其他类型。19/20 MP10 仍运行并单列补充评分覆盖层，不与完整20/20主层合并。

## 原生网络与空间距离

每片每算法仅拟合一个共同背景位置网络，不分 H/L 两套网络，不按样本总体 MP10 分组，不跨片拼接坐标。真实普通 Visium 阵列 nominal pitch=100μm；Secreted Signaling/ECM-Receptor=200μm，Cell-Cell Contact=105μm，边界容差1e-6μm。QC 留下的最近邻间距不改变平台距离尺度。Contact 是数据库分类与空间范围，不声称 spot 分辨率确认细胞膜接触。

COMMOT 0.0.3，commit d117445bc07eaa19109c7609e97d9e35b26e99ca；human bundled CellChatDB 三类蛋白 LR，多亚基min规则，原表达检测比例0.05，cot_nitermax=10000、eps_p=.1、rho=10、四组件权重各.25、不smooth。每个完整原生组件接收相同完整 LR 集，不按不竞争小批次拆分 LR。官方API不直接暴露停止诊断；保留有限、非负和原实现等价性检查，不把“函数返回”称为收敛证明。

SpatialCellChat 0.1.0，commit be88f300464cf970b36a6165ea7d2e41a47a7511；原生MERINGUE feature_adjP=.05、variable_both=FALSE、use_AGAN=TRUE、Kh=.5、Hill n=1、scale_distance=.0101；native filterProbability 的位置抽样观察概率分位QC（nboot100、seed666、quantile.95，minlinks10、minsource10、mintarget10）。没有 group P 二次筛选；位置过滤不是患者 H/L 差异检验。实际逐片参数与 LR 清单以 parameters.tsv 和各原生缓存为准。

全量扩展包含一项稀疏矩阵尺寸兼容修复：SpatialCellChat 原生 computeCellDistance 的 sparseMatrix 未显式指定维度，末尾孤立位置可使矩阵缩短。v2 仅在进程内函数副本中指定 dims=c(NC,NC)，保留全部位置、既定距离以及原生对角规则；未改安装库或 v1。CanDis_2026_Sample2 的 10061×10061 恢复为 10062×10062，原有子矩阵误差为0；两张 v1 试点逐值误差为0，额外概率等价检查误差为0。MERINGUE/LR 小缓存经精确匹配604特征、1841 LR后恢复；旧失败尝试保留审计。验证表见 logs/scc_compat_qa/。


## 主聚合及分母

原生推断保存含同位置项的网络。主汇总去掉位置对角线，仅改变跨位置归属汇总，不是重新拟合“禁止同位置”的网络。major/fine 分别从同一位置网络聚合，不能相加 fine 均值冒充 major。

COMMOT：T_cross=T_total 去掉对角线；F=QᵀT_crossQ；D=outer(colSums(Q),colSums(Q))−QᵀQ；G=F/D。主指标 abundance_weighted_mean_cross_spot_strength；F 是 total_flow。D 包含空间距离约束导致的零配对，不是只统计可达配对。D=0 为 NE。

SpatialCellChat：逐 LR 去对角位置概率 P_l；F_l=QᵀP_lQ；K_l=I(Q>0)ᵀI(P_l>0)I(Q>0)；B_l=F_l/K_l；主指标 total_mean_strength=Σ_l B_l（按官方 Visium Avg 的零支持处理），并保存 ΣF_l 和各LR支持分母。主指标不是 ΣF_l/ΣK_l。类型总支持符合条件、已计算而无正支持时为0；缺类型/状态/可定义分母为NE。两方法原始数值单位不同，不能直接相减或平均。

## 对照、份额和图

每片先 H/L 再汇总；保存 mean、flow、H−L、仅当两者均正的log2(H/L)、零类别、支持、覆盖和身份。全TME份额分母取同片、同方法、同level、同方向 H/L共同可评估目标；fine谱系份额只来自相应fine和共同可评估目标。分母标签全部保存；零总分母的份额为NE，无任意伪计数。

每方法每H/L有效片：major全部8群及fine全部8谱系各画High/Low/Both三张，主边用全部实际纳入LR总mean，仅显示Epi→目标和目标→Epi，方向分弧；TME↔TME、EpiH↔EpiL、自环、EpiOther展示边隐藏但原生背景不改。三图共享目标位置、顺序、颜色、等节点尺寸和同一edge max映射，不按top-N隐藏正边。不同片、谱系和方法分别标尺度；0无边、NE明示。透明矢量PDF、300dpi PNG及少量白底预览。每片另有所有major/fine的效应与双份额比较图。

## 患者、数据集与统计

核验患者包括SOURCE_VERIFIED_PATIENT及用户接受的来源优先冲突映射，保留证据标记。同一核验患者、同组织条件和治疗状态先以切片有限log2(H/L)中位数汇总；未核验记录单列，不自动视为独立患者。tumor/normal/unknown和不同治疗状态绝不合并；未知治疗仍标unknown。多片单元方向若同时有增强和减弱记mixed；零类别和NE按全部片保留。spots、LR、两算法都不增加生物学重复数。

dataset分别报告核验患者/记录、可评估量、增强/减弱/相等/NE方向比例及零类别。dataset主效应为患者或记录有限log2(H/L)中位数；跨dataset描述给每个有有限效应的dataset同等权重，取其median-logratio的算术均值，并逐一leave-one-dataset-out。身份、组织、治疗和评分覆盖层分开。H−L仅逐片/辅助展示，跨片原始G受组织几何和完整配对分母影响，不能直接解释成batch或生物学幅度。

探索性检验预先写于EXPLORATORY_TEST_SPEC.md：同dataset/组织/治疗/complete20、每对照至少8名核验患者有有限比值才做双侧Wilcoxon和2000次患者bootstrap中位数区间；按方法×dataset×组织×治疗×方向×major/fine分族BH。不足量明确NE，不放宽条件或伪计数。主要结论仍为描述，不以显著性筛选展示类型。

## 共同 LR 与诊断

从已完成native逐LR缓存，按method×dataset×组织×治疗×评分覆盖层取实际精确多亚基复合物LR交集，逐片事后重聚合；保存交集、原清单、保留数/比例和分层片数。单片层交集等于原集，不能作为跨片稳健证据。COMMOT仍处原完整LR竞争背景，不是共同集重拟合。SCC同精确复合物可能有多个原生ID，保留其原计数而不擅自去重。方法间数据库/覆盖差异不强行等同。

诊断包括全63输入的spot UMI/nFeature、Epi/目标Q、状态支持、基因/LR覆盖、组织规模、105/200μm可达配对及Q加权可达质量占完整D比例。结果只是质量、组成和几何的描述性关联，不是批次校正，也不保证片内H/L消除局部质量混杂。保持19/20、未核验身份、unknown组织、均值与flow方向不同、dataset方向相反、覆盖不足结果。

## 恢复、资源和可追溯性

scripts/resume.sh以锁保护单调度器，按切片×算法独立缓存/日志，先2任务测量、按余量升至3；按用户后续要求，在保留原worker PID和checkpoint的情况下提升为最多5个并行任务。新SCC worker为3线程，已有6线程worker完成前保持原设置；模型输入与科学参数未变。本会话预算24CPU/128GiB不变，模型加后处理共享22CPU affinity、预留2核审查，任务预算21CPU；按实际进程树CPU/PSS和服务器余量动态准入，项目总限48CPU/200GB。execution_plan、parallel_upgrade_verified和resource_admission_latest日志记录实际执行变更。其后按用户资源要求进一步保守准入：模型内存预留合计≤80GiB，观测PSS>64GiB不启动新任务，session硬预算仍128GiB。磁盘按本账户quota（soft500GiB/hard505GiB）而非文件系统df空闲量判断，每30秒审查，给运行任务按max(0.5,spots×0.0002)GiB预留输出并另留8GiB安全余量；不足时任务保持queued并写admission_waits，不删除既有结果。此策略保留所有模型参数和完整位置背景，可能等待配额/持久存储调整，不会以降采样规避资源限制。COMMOT组件算入预算，BLAS/OpenMP禁止嵌套扩张。失败记录具体原因，降低并发/复用checkpoint合理重试后继续其他片。不抽样、不切ROI、不删LR缩规模。logs/task_status.tsv是运行状态唯一写入源；summarize/report只读它。仅所有126任务终态且成功任务后处理完成才允许最终报告。完成后依次common_lr.py、summarize.py、report.py；轻量review.zip不包含巨型原生网络或RNA副本，完整缓存路径、图索引和恢复脚本保留。
