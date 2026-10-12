# Neighborhood ring comparison v2 — 方法与范围

本轮用户于 2026-10-12 明确选择 q05 环境＋c2l posterior mean Epi gate。环境使用 neighborhood_v1 已冻结的 subtype_q05（accepted cell2location v2_30type 同一模型的 1000 posterior draws），按细胞类型名称对齐30列。中心 gate 使用较新 K12 shared 的 posterior mean Epithelial / 全30类 >=0.25。Cycling 不额外并入。q05 归一化仅称相对组成，不是细胞数、存在概率或恶性比例。旧 RCTD gate 差异保留在 coverage63.tsv。

主机沙箱初始化报 bwrap 网络命名空间错误；所有读写通过已审核的命令执行。历史 venv 解释器指向已不存在的 /python-env，复用已有 Python3.12 与既有 site-packages；仅新目录安装 texttable 1.7.0 以补齐 igraph 的文本输出依赖。未升级原环境。每计算进程1 CPU，单细胞来源审计与聚类可同时各占1 CPU；不使用GPU。顺序读取大对象，未重跑 c2l/CopyKAT/DESTVI。

## 完整组织与输入隔离

原历史表仅含265002个有组成 spots，已用 surviving Seurat QS 对63张目标图像重新核验 tissue。所有目标 tissue==1 位置进入图，额外25676节点缺组成时仍保留拓扑；不填为观测零，不补远端邻居。早期基于QC-only坐标的试算已归档为 INVALID。Seurat晶格使用 row/col；六角偏移(0,±2),(±1,±1)。R2为二步可达二值集合去中心去R1；同一节点不因多路径加权。切片分别建图，禁止跨片边。几何6/12与有效观测数分别保存。完整环、缺观测和缺口测试实际执行。

所有切片265002个q05向量的内容hash重验，30列一致。图含290678个组织节点；中心必须有既有有效组成/表达QC，Epi gate通过且有效R1>=4、R2>=8。共同anchors只由这些条件决定；待审区域身份不剔除anchors。q1/q2是在完整组织图固定邻居内按spot等权平均组成，再sqrt。X12为两个sqrt向量拼接/√2。中心组成、坐标、CNV、MP10及来源标签不进特征。

## 抽样、拟合、比较

共同20,000参考anchors：dataset等额→患者层等额→切片等额，每层不足向有容量单位再分配，无放回，seed123。缺patient_uid为明确unresolved层，不能算成独立患者。不能把用户假定的独立样本/配对称为已核实患者；来源已知患者数按source_patient_identity_known另报。

KM在原始30/60维Hellinger中拟合，K=2/4/6/8/10/12/16，50次初始化，最多300迭代，无标准化。Leiden参考集拟合中心化、不缩放、不白化PCA，保留累计95%最少PCs；精确Euclidean 30NN包含查询自身构建邻域集合，稀疏共享邻域交集转Jaccard，<=1/15边删去，无向无自环。使用既有igraph0.11.9/leidenalg0.10.2的RBConfigurationVertexPartition，10 iterations，seed123；resolution=.1/.2/.4/.6/.8/1。这是显式记录的Seurat语义对应实现，不声称与其默认approximate Annoy逐位相同。

共同固定2000参考UID在各自原始Hellinger空间评估silhouette，不生成全量稠密距离。选参：数值有效、多群、收敛；与最优silhouette相差<=.02时取较小K。保留tiny群、来源集中和profile冗余，不按dataset混合程度删除候选。5次层级平衡80%参考行无放回子采样，seeds124–128；重拟合该分支PCA/聚类，但不重建rings。ARI是拟合稳定性，不是生物重复。

KM原生应用为最近质心，Leiden为冻结PCA中15NN等权投票。另存六方案common-transfer：KM在Hellinger，Leiden在PCA。参考查询排除自身做QC；参考原生标签不被转移结果覆盖。投票平局按对应投票邻居平均距离，再按稳定整数cluster ID。votes<.6是工程低一致性标志；nearest distance >参考leave-self-out Q99标远端，不删除。投票比例不是后验概率。每输入共享一套参考集UMAP（30邻居、min_dist=.3、seed123），只展示不选参；全量空间图使用原生应用标签。主一致性仅共同参考。Leiden K在2–20内才允许补一个同K KM；R12同K参考ARI<.8才触发一个PCA→KM桥接。范围外明确不运行。

## 区域候选和复核

CNV score=mean(native final CNA profile²)，不是mean(abs(log-ratio))，不套0.005/0.3。保留原prediction、baseline未外部验证、mixture非收敛、profile可用性。aneuploid+Epi支持生成tumor-supported候选；diploid不单独证明正常。Epi不足不能证明非肿瘤；所有aneuploid含低Epi节点单独保留肿瘤证据mask。没有审定H&E字段，region_final留空。按同证据连通组件保留每个小岛及复核UID，不平滑跨洞。ROI review按明确审阅字段回填，可审阅后保留unknown，空白不能自动接受。

候选边缘仅是CNV证据到未定邻区的几何边，不命名病理肿瘤–基质界面。当前depth无效/NA。复核后：已知tumor内邻接已知外侧组织为inner boundary；外侧已知组织邻接tumor为outer，各自0 hops，分别在本侧BFS。最近未知区/组织缺口距离不大于已知界面距离时保守截断；不把组织边缘设0。只有明确stroma复核标签才使用tumor–stroma。更新区域身份不重拟合共同特征/标签。

## MP10和来源验证

聚类脚本不读取空间MP10列。原始Top20与Direct AUCell路径/genes存在性单独核验：全部原始表达层features降序平均秩，K=ceil(.05*n_features)，sum(max(0,K-rank))/sum(K-r), r=1..min(可用signature数,K)，不z-score、不换基因。不同可用基因/表达universe不视为统一尺度。正式M06等待真实区域复核与六工作标签冻结；给出全部、tumor-supported、normal、reviewed unknown、完整环及高Epi的患者级描述/效应。高Epi>=.7是本轮盲态预定的“至少70%总上皮组成”敏感性，并非文献恶性cutoff；覆盖在input_contract中记录，不随MP10改变。

实现的M06以切片内效应为点估计，按真实患者等权汇总切片、在dataset及基因可用性层内bootstrap患者；无spot独立标准误/P值，不做跨dataset分数混合。调整限中心Epi与log1p UMI；深度调整及不调整使用相同可靠深度子集。少于5位真实患者仅描述。区间是探索性pointwise，不宣称多重检验控制后的显著性。六方案同数据，不是六次独立验证。

M07必做并保留：当前有限单细胞来源审计在每类型最多300个细胞、内部source sample平衡，使用原始Top20和原始AUCell公式，全参考RNA universe。它与空间feature universe不同，不能横比绝对分数。一般pseudo-spots先混合counts再binomial thinning至目标7005UMI上限，再打分，保存细胞成员。按MP10本身划低/高仅检查混合后评分保真性，不是独立生物学验证。既有DestVI v8为19基因G0-UCell（缺HMGCS1），只作辅助，不可冒充本轮20基因AUCell来源验证。正式归因仍需冻结环境标签对应的composition-only、实际深度/基因universe匹配、恶性/正常细胞独立身份与高Epi真实关联复核。MP10不一定是上皮专属，中心混合效应未被排除。

## 文献做法—本项目采用—自定部分

| 来源与核验 | 文献做法 | 本项目采用与边界 | 自定部分 |
|---|---|---|---|
| [BLCA L1](https://www.nature.com/articles/s41467-026-75399-6)，Methods网页受站点限制；本轮未独立核验具体阈值 | 任务文件描述病理非肿瘤参考/ESTIMATE最高100以及mean绝对log-ratio；这些细节在未核对原文前仅作待核验条目 | 只复用已运行CopyKAT，不重新选择参考，不把本文score等同现有MSQ | Epi gate和候选ROI规则 |
| [TNBC](https://www.nature.com/articles/s41586-026-10469-9)，2026-10-12核验Methods | CNA cluster平均squared aberration>.005、平均相关>.3并结合已知事件；样本平均profile作相关参照 | 不把其阈值套现有任意CNVscore；相关参照不是CopyKAT正常基线 | 本轮不增加阈值分类 |
| [BLCA前沿](https://www.nature.com/articles/s41698-026-01417-x)，核验Methods与Results | SpaCET最大组分作初始标签；另用STdeconvolve；malignant>.7为后续敏感性 | 保留自动候选与形态复核层级，不把total Epi当malignant；GU复核细节未完整复核不声称已复现 | 本项目不新增全队列SpaCET |
| [Morph](https://www.nature.com/articles/s41586-024-08087-4)，Methods站点访问限制 | 任务文件所述>=50%恶性形态像素面积规则待原文细节核对 | 形态面积不能等同Epi或恶性细胞比例；不套删除小巢规则 | 保留所有连通小岛 |
| [SpaCET官方](https://data2intelligence.github.io/SpaCET/reference/SpaCET.identify.interface.html)，v1.4.0；[源码](https://github.com/data2intelligence/SpaCET/blob/eed9666a0d76507d46e70e5a168da2504086db0e/R/interaction.R) | malignant cutoff .5；六角邻接tumor的外侧非恶性spot标interface | 已核验源码；total Epi不是malignant，外侧不自动生物学stroma | 本项目双侧边界、unknown与truncation |
| [CellCharter](https://cellcharter.readthedocs.io/en/latest/generated/cellcharter.gr.aggregate_neighbors.html)，核验接口 | n_layers列表指定exact hops；整数包含layer0 | 借鉴语义，直接稀疏图实现并验证6/12 | 不声称复现其完整框架 |
| [Hellinger](https://link.springer.com/article/10.1007/s004420100716) | 组成生态数据的sqrt变换 | 先环均值再sqrt | 30类型、独立环、等权60维拼接 |
| [Seurat FindNeighbors](https://satijalab.org/seurat/reference/findneighbors) / [FindClusters](https://satijalab.org/seurat/reference/findclusters)，核验官方接口 | SNN共享邻域Jaccard、prune，Leiden algorithm4；指定实际降维/graph | neighborhood PCA建图，UMAP只显示；语义对应Python实现 | 95%PC、30NN、有限候选和比较规则 |
| [Harwood](https://link.springer.com/article/10.1038/s41467-024-52167-y) / [Polonsky](https://link.springer.com/article/10.1038/s41467-024-51186-z) | 邻域组成/表示接图聚类的背景先例，未用于推导本轮数值参数 | 本轮未独立核验其Methods；不声称平台/环/参数原样复现 | 整套六方案是项目设计 |
| [cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html)，核验官方教程 | posterior abundance summaries与q05 | 同一冻结q05归一化；大类仅组成求和，不称总量q05 | q05环境＋mean中心gate经用户确认 |
