# MP10 假说后续：实际结果

本轮完成 M0–M5 全部可运行分支。新增模型共1473项，924 completed、549 not_evaluable；NE不是阴性结果。主Epi10项、桥接18项和NC9项均可评价。详细状态见module_status.tsv。

## 1. 覆盖与身份

| profile_type | n_specimens50 | n_specimens30 | n_datasets50 |
| --- | --- | --- | --- |
| Arterial EC | 4 | 5 | 1 |
| Capillary EC | 12 | 14 | 4 |
| Epi_all | 39 | 39 | 5 |
| H05 | 24 | 30 | 5 |
| H10 | 33 | 36 | 5 |
| H50 | 39 | 39 | 5 |
| L50 | 39 | 39 | 5 |
| Tfh/CXCL13-like | 23 | 25 | 5 |
| Treg | 32 | 33 | 5 |
| Treg_CY | 28 | 31 | 5 |
| Treg_NC | 30 | 32 | 5 |

均为记录患者单位；其余身份仍derived/provisional。a/b从不同时入独立样本推断。40标本用于描述，主选择39个代表。原Epi亚群标签缺15224/116448(13.1%)，单列，不新造亚型。

## 2. 上皮支持超出了多少评分基因？

新Epi_all中mTORC1集合200成员与Top20重叠10，M20可用184/190；CollecTRI SREBF1/2 M20分别覆盖73/117和42/60targets。四分支M20覆盖为FA 7/7、甾醇17/18、摄取运输31/44、prenylation10/10；均达到评分门槛。详细交集与驱动基因见gene_overlap.tsv、regulatory_gene_membership.tsv.gz及Epi_mTORC1_driver_genes.tsv.gz。

MP10→SREBF2(M20)：β=0.402, 95% CI [0.155, 0.649], P=0.00227, q=0.0113, n=39。MP10→甾醇(M20)：β=0.535, 95% CI [0.286, 0.785], P=0.000118, q=0.00118, n=39。甾醇在所有LODO、Q95、a/b、RNA质量和原cluster0比例调整后仍同向；SREBF2所有LODO同向，但DoRothEA和独有靶集证据更弱。

SREBF1(M20)：β=0.276, 95% CI [-0.042, 0.594], P=0.0869, q=0.217, n=39；mTORC1(M20)：β=0.264, 95% CI [-0.079, 0.607], P=0.128, q=0.255, n=39。Full更强，去共享基因后降低了独立上游支持。配对层面mTORC1 M20的H50/H10/H05 NES约2.60/2.61/2.74，仍有排序富集，H05已完成全通路数值修复。**配对变化与跨标本关联是不同层级，不能互相替代。**

SREBF1/2→NECTIN2：β=-0.063, 95% CI [-0.357, 0.231], P=0.665, q=0.831, n=39；β=0.027, 95% CI [-0.290, 0.345], P=0.863, q=0.863, n=39。未建立直接转录桥。旧NECTIN2上皮配对上调保留，但MP10 Q50→Epi_all NECTIN2本轮q约0.465，不能从配对上调推出跨标本线性关系。既有SREBF1 SCENIC去Top20后减弱，SREBF2 Epi motif-pruned regulon仍NE；未重建SCENIC。

## 3. Treg：组成还是状态内部？

13204 Treg中S/G2M共6227，G1共6977，无无效分类。S/G2M约47.2%，但MKI67/TOP2A多数细胞仍零，分类仅是operational Phase。

MP10→All-Treg OXPHOS：β=-0.239, 95% CI [-0.617, 0.138], P=0.204, q=0.733, n=32。同30标本All-Treg：β=-0.266, 95% CI [-0.702, 0.170], P=0.22, q=0.944, n=30。同30标本NC-Treg：β=-0.271, 95% CI [-0.705, 0.164], P=0.211, q=0.894, n=30。点估计相近，不能凭某一设置显著与否宣布组成解释。NC加入线粒体比例后β约−0.271；去InhouseData后All/NC β均约−0.069，提示数据集依赖与不确定性。

计数模型的排序富集则保留：All/共同样本All/NC/CY OXPHOS NES约−2.73/−2.73/−2.59/−2.67；对应G2M约2.39/2.17/2.11/1.48。NC OXPHOS在全部LODO仍负向，去InhouseData后NES约−1.63、q约0.00246，强度降低。核编码多个呼吸链环节多负向，线粒体编码亚基多为正向；单基因全基因FDR支持不足。

因此**结果不宜归为仅有cycling组成改变，也不足以确立跨标本代谢功能变化**。GSEA强调许多基因的协调排序，HC3分数模型衡量有限样本的分数关联；二者相互补充而非独立重复。不得用NES差值计算组成解释率，不能把OXPHOS RNA下降写为呼吸下降或抑制功能增强。

## 4. 接触与脂质桥接

18项主桥接和9项NC均无本族q<0.05。NECTIN2→All-Treg OXPHOS：β=-0.169, 95% CI [-0.694, 0.356], P=0.514, q=0.923, n=32；NECTIN2→cycling的方向为负且不确定，未得到本轮新增机制支持。

mTORC1(M20)→cycling的log OR β约0.183，OR约1.20，主q=0.0797；甾醇(M20) β约0.193，OR约1.21，q=0.102。所有LODO均正，但不能用个别敏感性q替换主q。四脂质分支对OXPHOS未提供主FDR支持。Treg次要家族45项中3个程序覆盖不足，剩余程序也无主FDR阳性。

已有NECTIN2–TIGIT/PVRIG/CD226 receptor表达、Treg总丰度、正向NicheNet低target和反向NicheNet排名仅作背景；本轮未重新计算NicheNet，也不声称TIGIT主导或免疫抑制已经确立。

## 5. 内皮与Tfh

保留旧Tfh和Capillary EC状态、丰度、TF及受体结果。Capillary EC主12/4datasets，其中PRJNA662018贡献7，移除该dataset原模型不可估计；旧MYC V1/TNF正向排序值得继续，ANGIOGENESIS主负向不显著、Q95负向，不能称血管生成增强。

新HC3主模型受singleton dataset(h=1)或冻结signature覆盖约束，Tfh/Capillary主分数比较NE；这不是否定旧GSEA。Arterial EC主4/1dataset、df不足保持NE。跨细胞类型仅描述，未作不满足可比背景条件的交互，也不声称Treg特异。

## 6. 假说更新与下一步

当前优先级：先巩固上皮SREBF/甾醇转录状态，并用p-S6/p-4EBP1、成熟核内SREBP和甾醇通量检验上游；Treg需独立cycling、呼吸与抑制功能读出。接触与代谢两支均等待有区分力的证据，不据本轮筛出唯一机制。HD首先检验同切片局部NECTIN2接触距离与Treg状态，并与较大范围代谢环境关联区分；毛细内皮应扩充dataset覆盖。空间/代谢物转移均pending。

## 7. 交付与限制

先看F06_bridge_OXPHOS、F05_All_NC_CY_decomposition与F05_gene_rank_decomposition（必须并列）、F03_Epi_continuous、F07_LODO_OXPHOS以及decision_table.tsv。每图有直接表及PDF/300dpi PNG，图例拆分。模型、分层count/GSEA、旧证据和完整NE均保留。

当前新增目录约91.6 MiB，加代码约0.1 MiB；最终精确值见delivery_inventory.json。峰值RSS约41.76 GiB，发生于原对象读取；其后使用小型聚合。读取/聚合过程有两次类型/索引修复重试，日志完整保留，无旧结果覆盖。

没有运行湿实验、新HD、全量Seurat/聚类/CNV/注释、CellChat、SCENIC或NicheNet；没有推送远程仓库。独立验证、因果方向和最终活性分界仍未确定。

输入版本限定：历史v1 pseudobulk实际成员相对分位数表报告数，H50多8、L50少8、H10多2、H05多1个边界细胞（涉及11个不同细胞）；源评分/分位点导出精度与解析造成边界差异，逐标本每状态最多1个，未改变>=30/>=50资格。历史成员可由当时原输入文件的比较规则逐一复现。本轮保留冻结membership和原DE，不把二次导出的分位数表重新比较结果替换进去；新主关联用真实Epi_all和Treg profile，不以这些边界H/L细胞选择标本。详细见frozen_state_count_audit.tsv和quantile_reexport_roundoff.tsv。

G2M与Phase分类共享信息，NC中G2M正向属于原发现的拆解，不是独立的cycling验证；也可能包含分类界限内连续状态变化。
