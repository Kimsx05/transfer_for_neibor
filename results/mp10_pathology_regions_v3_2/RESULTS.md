# v3.2 实际结果

完成 **63/63 张、7 个数据集、265,002 个原始 spot** 的代码运行、点表、重新生成的 ROI/边段及 PNG/PDF 空间图。入口：[全队列图形索引](index.html)。旧 v3/v3.1 保留，原始分子数值不变，没有开展新的 MP10/TME 分析。

**关键限制：本轮没有获得足够支持的宏观病灶外缘。** 新版本不再把“未进入癌细胞富集 ROI”解释为已在病灶外，因此原版大部分 outside 改为范围未定。最终 63 张均不宜直接使用旧 margin 作“宏观肿瘤外缘”分析。这里的0张不是说生物学上没有边界，也不是没有计算成功，而是当前证据和明确保守规则不足以认可这些候选分界。已输出候选内部范围与具体未定范围；**不能把本版说成已完整确定63张的宏观病灶轮廓**。

## 位置与单一区域

| 项目 | 数量 |
|---|---:|
| 计算支持的 inside 候选 | 134,267 spot |
| 有正向非病灶背景支持的 outside | 5,644 spot |
| 宏观范围未定 | 125,091 spot（47.20%） |
| 原 outside → inside | 629 spot，其中 BC18 为503 |
| 原 outside → 范围未定 | 123,296 spot |
| 独立保存的新 ROI 连通组 | 442（含小候选，不等于442个独立肿瘤灶） |
| 单一区域组织块 | 75个全 inside、17个全 outside；多数是微小断片 |

**整片支持 outside 的是 GSE285715_N1**，1,417点、6个观测组织块；依据为已有癌旁来源与尿路上皮/平滑肌背景裁决，保留其原 SpaCET 高值冲突，未将其强制清零。没有整片被本轮强行改成 inside。

较大的全 inside 组织块（以下20点仅为报告筛选，不参与位置规则）：CanDis Sample7:T2（202点）、Sample6:T3（39点）、Sample12:T3（37点），以及 GSE319536_BC11:T4/T8/T12（51/86/36点）。每块由自身观测连通的局部支持及内部区室构成，外缘未观察到，不称 deep core。全 outside 较大块另有 B11:T4（56点）、BC04:T10（98点），有平滑肌/非病灶背景支持。

完整列表：[组织块状态](tables/tissue_component_summary.tsv)、[单一区域组织块](tables/single_region_components.tsv)。不能将微小断片数当作病例数或肿瘤灶数。

## BC18 与 T3

本地复算与前述线索一致：

| 原v3.1分区 | EPI中位数 inside/outside | SpaCET malignant中位数 inside/outside |
|---|---:|---:|
| BC18 | 1.1655 / 1.2276 | 0.7669 / 0.3216 |
| T3 | 2.1229 / 2.0377 | 0.5334 / 0.2333 |

这些是描述性检查，不是独立分割准确率。EPI已参与旧分数；CNV、SpaCET和marker亦多来自同一表达数据。

- **BC18**：作者 S1 区分 Tumor Cell 与 Tumor Stroma，中央/下部部分旧 outside 与二者定性重叠，支持“原分界切入病灶相关组织”的疑虑。按已记录的同片作者背景和局部域交错关系，503点扩展进入候选内部；最终2,003点 inside、844点范围未定、0点已支持 outside。没有把 Others 或整张切片全部吞入。作者图不提供可直接匹配的 barcode 标签，因此未将这503点声称为病理精确确认；其余范围仍未定，不计算宏观外缘分带。[对照图](by_dataset/GSE319536/GSE319536_BC18/01_location_comparison_HE.png) · [作者S1同片对照](by_dataset/GSE319536/GSE319536_BC18/05_author_S1_comparison.png)
- **T3**：最终1,896点 inside 候选、1,445点范围未定、0点已支持 outside。广泛上皮程序和 bladder cancer 来源不证明全片都在病灶内；原红蓝分界主要伴随恶性组分估计变化，不能直接作为宏观外缘。没有恢复 Tumor-419a↔T3 对照关系，未转移论文示例标签。[对照与H&E](by_dataset/GSE285715/GSE285715_T3/01_location_comparison_HE.png) · [无外缘状态图](by_dataset/GSE285715/GSE285715_T3/02_bands_distance_reliability.png)

GSE319536 的16张现有切片全部提供作者 S1 同片定性对照。BC01、BC13 的作者图可见 Normal Stroma 等不同区室，但没有逐点标签，不能据低分辨率图给未确定的 BANKSY 域伪精确定位；这些仍是本版的具体限制。S9 相对距离分区不等于本项目宏观 ROI 或200μm分带。其余来源和映射限制见[文献对照简表](tables/reference_comparison.tsv)。

## 小 ROI 与边段

共检查原版少于10点的284个 ROI；**57个缺乏局部支撑、依赖同域借分的小 ROI，共87点，撤回为范围未定**，分布在26张切片。没有按大小批量删除。

重点例子：B15为11个ROI/22点，B11为2个/2点，T4471438为1个/1点，NPJ_sample1为1个/1点；BC01为3个/6点。B16原来的5点和1点小 ROI 保留：没有足够依据仅凭体积删除；其6个“仅借分跨阈值”点没有形成应撤回的独立无支撑病灶，已有主体观测连接提供位置背景。B16另外1个原 outside 点因完整直接邻域支持归入内部；不将这些位置继承理解为恶性身份确认。

所有新位置重建 ROI/边段/参考ID。小ROI周围原margin已不再当作宏观外缘margin；同时发生了全局外缘证据重判，**不能把邻近分带的全部变化都归因于87个点**。逐ROI表包含局部借分、原barcode、修改依据及一圈邻域新旧分带计数：[小ROI处理表](tables/small_roi_decisions.tsv)。

最终检查还发现 T4465935 的4条候选外缘依赖一个内部非病灶支持区通过 QC 缺测通道被判作通向外部。原始全组织坐标不支持这一连接，已将其重新归为内部区室接触，撤回相关53个曾可判断的距离带标签。这是缺测拓扑检查的修正，未改变位置或分子值。[该片图形](by_dataset/Inhousedata__UNVERIFIED/Inhousedata_T4465935/index.html)

## 距离、图形和交付

最终无认可宏观外缘，所以主 signed distance 和100/200/300μm距离带均不产生数值 margin/core。inside 仍红色显示并注明“外缘未观察到”；范围未定另示，不以整片灰色替代已支持的内部位置。原 v3.1 距离和分带保留为对照列。48个有图像/原始组织支持的技术缺位用于几何模型检查，无新增分子观测。

每张4套核心 PNG/PDF：旧/新位置与H&E、分带/距离/可靠性、EPI多marker/SpaCET/CNA、域与小区域诊断；GSE319536另有16套作者图对照。另有7个dataset总览与逐片浏览页。KRT18有5片可用、58片 unavailable，未填0；所有基因采用共同队列色标。

核心文件：

- [逐片状态汇总](tables/section_summary.tsv)
- [全队列标签 Parquet](tables/all_spot_labels_v3_2.parquet) / [TSV.gz](tables/all_spot_labels_v3_2.tsv.gz)
- 各片 `spot_regions_v3_2.parquet`：完整原分子字段、旧/新位置、分数、ROI、边段参考、距离/状态。
- [分子证据分布](tables/evidence_distributions.tsv)、[旧/新边界两侧描述](tables/boundary_evidence.tsv)。无新MP10/TME比较。
- [保留的身份映射](tables/preserved_identity_crosswalk.tsv)：B5/B18仍同为GSE238145::P1；T1/T2/T3仍为三个未核验切片单位，未推断患者或配对。
- [方法与复现命令](METHODS.md)、[实际主脚本](scripts/run_v3_2.py)、[配置](config.json)。

全队列主运行两片并发耗时约6分52秒（不含输入核查、试跑、图形复核及汇总）；单个子进程记录峰值RSS约710MiB，远低于本任务资源限额。最终校验覆盖原始数值/缺失、barcode/坐标、同块几何、无外缘时禁止分带、身份映射和HTML链接。这里交付的是可复现的保守位置修订与明确未定范围，不把未解决的宏观外缘问题写成验证通过。
