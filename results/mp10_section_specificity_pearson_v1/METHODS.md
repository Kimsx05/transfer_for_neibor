# 逐切片 Pearson / 样本差异 (执行前定义)

## 输入与冻结对象
直接读取 full_v1、subtypes_v2、p95_front_programs_v1 的63片manifest及完整barcode。主MP10为EpiMP10_Top20_AUCell；query、tumor_epi、MP10_P95、pEMT_high_front、v3.4组织块/精确signed distance、R12 accepted_label均原样保留。不是重新按纯度、CNV或niche定义Epi。正常来源只分析既有Epi与可用候选距离，不定义肿瘤P95。距离坐标是上游nominal Visium 100μm点间距的x_um/y_um，非图像像素。

29个非Epi候选采用同spot全部30项posterior mean之和作分母。全片原有效组织内P95且比例>0，保留并列。7个已兼容亚型优先直接复用subtypes_v2目标及距离；其余同组织块精确近邻补齐。无同块目标距离NA，自身为参考点保留0。d_candidate、d_MP10和d_front分列。

## 85个程序与评分
原80套成员原样继承，新增Budden 4主+1去MP10重叠敏感性，共61主+24敏感性。PNPLA2不加入TG24。胆固醇27与Top20重叠9基因。CCN2/CTGF等按表达层实际可用名映射，保留原名与实际使用符号。OXPHOS既有集合与Top20无重叠，因此其原始版本就是无共享成员版本，不另制造第86个程序。

所有85项AddModuleScore逐片在原有效tissue spots统一运行。没有兼容SCT缓存，故从对应完整稀疏RNA补做Seurat5.5.1 SCTransform(vst.flavor=v2,ncells=min(5000,nspots),seed=20261104,min_cells=5,conserve.memory=TRUE)。glmGamPoi未安装，使用SCT原生R拟合。return.only.var.genes只影响scale.data；实际评分显式用assay=SCT,slot=data，即log1p corrected counts，全SCT可用基因池。验证data与log1p corrected counts一致。AddModuleScore nbin24,ctrl100,seed20261104；不分P95/距离/候选重跑，负值保留。新列AMS_SCT__开头。完整SCT表达不重复持久化，保存评分、成员和参数；复用评分缓存可重画。

旧80项UCell直接复用，新增5项沿用完整RNA universe、UCell2.14.0,maxRank1500,ties average,chunk50。比较二/三AMS为主，UCell为评分敏感性，两方法逐程序取同一有效spot交集、相同分组/距离范围。不能拿旧β/partial r称作新Pearson。既有重叠邻域模型继续保留其2000μm空间块等定义。

## 比较与不确定性
比较一：P95与全部nonP95距离中位数、均值、零比例；空间推断的效应明确为均值差，中位数差描述性。原base匹配对逐对复用，匹配距离差描述性，未再匹配候选距离或程序。

比较二：原始逐spot Pearson r及独立观测假设下双侧P，不作组成调整或范围择优。n≥2且双方有变异即保留r，n=2时P=1；n<10标very_sparse，10–29标sparse。全范围r与共同支持范围r分别保存。负r是总体负线性关系，不证明每一段单调。

比较三：在同一候选参考集合中计算r_P95、r_nonP95、delta_r。full_range仅描述；common_support为两组逐程序有效距离[min,max]交集，另固定common_support_0_1000供跨片可比描述。范围不重叠不回填。最大观测距离断档及>200μm标记保留，不插值。相关差异不等于μm斜率差。

空间不确定性以组织块×500μm方格为簇，seed20261104+manifest序号，1999次块bootstrap；每次重抽块保留multinomial重复权重，预聚合n,sumx,sumy,sumxx,sumyy,sumxy计算相关。共同一次抽样权重同时用于两组及两评分，delta区间来自逐次r差，非分别bootstrap。每组≥8个贡献块、≥1600有效重复才标可估计区间；未满足仍保留r/P/n。正式组间空间P用Fisher-z差及同一块内两组相关影响量差的cluster-robust方差、有限簇修正和t检验；不采用独立spot Fisher-z检验。full_range不作正式空间差异P。

补充partial correlation控制Epi fraction与log1p(nCount_RNA)。在精确d_front有值的同一子集上，另存加入前沿前/后的partial r和n，缺前沿不影响主Pearson。比较四在全部冻结Epi上计算连续MP10–候选距离Pearson和对应partial r。

普通Pearson P和q没有消除空间自相关。检验族按score_method、主/敏感性、population、range_type及来源分别定义；片内候选×程序全范围BH，以及相同族跨全部切片BH并列。比较一29候选、比较四29候选、比较三全部候选×程序空间P、仅前沿带九niche均独立成族。保留旧context q为旧版字段，另算包含新增程序的context q，不能冒充本轮Pearson全扫描q。

## 身份与图示
以section先展示；已知重复患者/analysis unit内在Fisher-z等权合并，绝不用spot数作患者权重。缺身份保留逐片描述，不创建患者。正式异质性主要在dataset内部，使用统一0–1000μm共同支持结果及空间不确定性；即使差异得到支持也称样本/切片异质性，有限重复切片不足以命名患者特有表型。完整单位及患者覆盖保存。

散点使用原始AMS(允许负值)和完整原始距离，统一两组坐标；主图标原Pearson r/P/n及delta的适用范围。固定核心程序面板与按统一覆盖/效应规则选的探索性散点分开，全部85×29组合可由spot源数据重画。沿用jins-routine-plotting规范，空间/矩阵/散点使用已有Matplotlib样式，白底矢量PDF、300dpi PNG、Arial约12pt、线宽≥1pt、图例三版。

参考：[Budden原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC13218938/)、[作者脚本](https://zenodo.org/records/18790615)、[Seurat AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)。本轮是文献对应的SCT AddModuleScore+原始Pearson适配，并非完整复现原论文研究设计或其全部默认环境。

## 跨样本汇总的实际实现
同一analysis_unit或核实patient的重复切片先在Fisher-z上等权汇总。区间为已有切片空间bootstrap z标准误的条件性合成(sqrt(sum(SE²))/片数)，假定不同切片空间误差独立；这不替代患者间重复。dataset内Q异质性检验仅在≥3个有可用空间SE的真实单位中估计，精度来自空间bootstrap而不是n_spot−3。Q检验的加权零假设均值与表中单位等权描述均值分开；单位对的Fisher-z差区间也是条件空间不确定性，不作为临床患者特有表型验证。未知unit仍保留逐片及dataset描述。

新增五项的旧上下文补充也包括前沿位置调整的P95环程序差、pEMT-high front支持、原base/front匹配参考集合的互斥距离曲线和图环。旧80项β/区间保留原义；含距离二次项的模型使用联合总体P进行上下文校正。85项上下文汇总继续用既有1999次单位及dataset层级bootstrap和dataset均值t检验，未用新Pearson P替换。

原始Pearson热图不以空间区间可用性筛选。|r|≤0.1只用于方向分布的描述性“接近零”分类，不是等效性检验或预设生物学无效区间。散点展示在固定核心面板中优先n≥30、各组≥8空间块，按|delta r|和相反方向挑选；覆盖不够时回退固定Treg/S1P与Matrix_CAF/TNF-EMT例子，仍保留真实r及稀疏标记。完整范围散点另附逐dataset补充分册，共同范围图不冒充完整范围。

单细胞表达来源复用上一轮同一19,053细胞样本(每样本×原注释类型最多60个)，完整RNA宇宙UCell，只用于来源支持；新增程序及核心解释可查看逐样本/亚型结果。没有构造未核实的空间–单细胞配对。

执行采用受限BLAS/OpenMP线程。初始6个逐片R评分任务，后续在资源限额内对后18片增加6个任务以缩短等待；统计4进程、上下文3进程、绘图2进程。已完成评分优先复用，最大同时12个R评分进程；按片对象读取，未复制全量Seurat或模型。SCT conserve.memory的中间日志“not return corrected UMI”来自初始残差类型none步骤；Seurat后续明确执行correct_counts，最终SCT data与log1p(SCT corrected counts)进行了数值核对，未用原counts代替。

空间区间以已估计的候选核心、冻结MP10-P95成员和评分为条件，不在每次bootstrap中重新训练SCT/c2l或重选P95。完整范围普通Pearson图与共同范围正式比较分册均可按源数据重画。MP10空间色标沿用上一轮统一上限0.1222405941465743，精确前沿显示范围±1000 μm；热图的显示截断均不修改数值表。原上下文q列明确加legacy80标记；新增85程序上下文q另命名。原MP10/P95概率和R12概率等非基因集端点使用独立端点族，不与新候选×程序Pearson扫描混合。

协变量完整点数5–9时partial r保留但明确标记sparse/descriptive。尤其加入前沿后剩余自由度极少时，接近±1不是稳健性证据；主Pearson可计算值不因该限制删除。
