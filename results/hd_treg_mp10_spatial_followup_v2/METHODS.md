# 实际方法与复跑边界

本轮为v1结果后的探索性复核。原counts、QC、评分、参考与结果只读；配置于2026-11-17 14:18 CST冻结，之后只进行了记录在common_support_reconstruction.yaml中的几何去重修复。没有基于MP10关联调signature、半径、ROI或高分门槛。

## 输入身份与版本

`tables/input_manifest.tsv`列出五片真实H5、坐标节点、18085输入基因、数据集和未知患者关系。Cervilla native8输入来自上一轮native8目录；cell计数来自segmentation_assignment_v2/Cervilla/cells.h5ad，旧可信上皮来源由hd_mp10_treg_endo_precommunication_v1/tables/epithelial_region_membership.tsv.gz限定。源2 µm归属映射用于恢复cell表达footprint，转换为row8=floor(row2/4)、col8=floor(col2/4)，不移动原点。

远程仅核对指定基线提交的三个小文件（旧共同中心脚本、配置、主表），通过GitHub只读接口；未下载新原始数据、未上传数据。精细参考使用subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad及v1固定成员。混合参考为17987个唯一symbol，完整参考背景秩；HD为18085原基因背景。差异继承v1，不将参考/HD基因背景混称相同。

## A：配对设计与统计

共同中心T32必须在两个200 µm来源下都有可用独立上皮结局及共同协变量，无共享原始T32 bins。50 µm不放宽资格。协变量固定顺序：T-lineage、T侧nCount、输入圆内可观察比例、原T区域上皮/成纤维背景平衡、bin上皮深度、cell上皮深度。两侧深度同时进入两套模型。面积无法可比，均不使用。竞争状态加入Tfh及Activated_T。

平均并列秩→总体标准差标准化（ddof0），在固定共同集合上一次完成；bootstrap不重新计算秩。按固定顺序计算每列对已纳入设计残差方差，>0.1才保留，相当于预定VIF<10原则。共同模型使用完全相同的X、相同Treg向量、相同中心顺序；未按P或方向删变量。缺失采用两结局及共同协变量联合complete-case，缺面积不填0。

区块由原component及floor(x/400),floor(y/400)定义，另作800 µm。999次有放回抽取区块；两模型使用完全相同索引，直接保存每次βcell、βbin、差值。CI为配对百分位区间。唯一主要新增P针对Ydiff=zrank(cell)-zrank(bin)，在不含Treg的受限共同模型下拟合零假设，用HC2残差、每区块一个Rademacher符号进行999次wild bootstrap，P=(1+|null|>=|observed|次数)/1000。它是依赖区块内相关、区块间近似独立与对称误差的片内探索性零模型，非患者检验、非普通bootstrap尾部P。百分位CI与该零假设P不互为严格反演。

800 µm、竞争状态、原热点排除和新增B原完整模型P进入统一辅助BH家族；来源各自β仅CI，不用两个P之差判断来源差。LOBO不重新选择变量/重做秩。原热点定义是v1两个候选256 µm矩形并集；本轮512 µm显示框不参与排除。

共同支持重建将整个可信上皮cell限定在两来源实际8 µm成员交集及单个原始64 µm网格内，删除与T32共享的来源。触碰的cell footprints按8邻接合并，避免两个历史cell区域在同一8 µm位置重复使用bin counts。整cell counts从原h5ad求和；bin端从相同footprint的原H5求和；两端不按面积分摊。170中心及其585对评分单元均保存原始成员。各T中心内部来源成员唯一；不同重叠邻域可引用同一评分来源，不被解释为独立重复。与原共同中心相比测量单位/中心发生变化，因此结果仅作有边界的敏感性。评分继承18085背景、平均秩、K905及20基因归一化，nFeature>=905、MP10检出>=5、上皮marker>=2；两来源同时满足才进入均值。

## B：混合校准

11基因、UCell maxRank1300、平均并列秩、T32资格及旧联合规则均不改变。旧门槛0.1995538920。无MP10、TIGIT、NECTIN2、CD96、NECTIN1、PVR、CD226用于定义Treg。文献signature只引用v1，未新选基因。

先按v1 RNG顺序重放原70%上皮2700混合，核对选定原始计数和全背景秩。扩展21600混合来自固定18,289训练/2,157留出cell、46/12样本。每表达profile按类型取5cell合并原counts再归一化为RNA概率，参考抽样不代表组织比例。背景Epi40/70/90%，及40%Epi+30%Fib或Myeloid；其余T RNA。竞争普通CD4、Tfh/CXCL13、耗竭CD8；Treg占T RNA为0/0.1/1/3/1。总UMI500/2000/5000，60profile重复。训练/留出来源完全分离，派生混合不跨分区。保留预设概率比例和用独立随机流条件归因的实现RNA来源量；后者只补充真值记录，不改变抽样表达、评分或门槛。

资格、原UMI、检出宽度、核心基因、秩窗口、旧/候选门槛及完整联合规则分别保存，并同时给顺序首次失败与单规则失败。UCell0例子逐项核查检测marker秩是否>=1300；不能把秩窗口截断等同无RNA。

唯一候选门槛是在训练无Treg混合中，各背景×竞争×深度满足T资格/多基因/核心的子集95%上取整分位的最大值。如此固定全局门槛是保守操作规则，不等于真实HD FDR。训练分层及留出完整规则性能全部输出。接受标准事先冻结：所有有支持负条件FPR<=5%、条件Wilson上界<=15%、整体阳性检出至少提高5个百分点且每条件不降低超过2个百分点。候选未通过，保留旧规则。没有第二次阈值搜索。500–5000 UMI之外不外推校准。

活化CD4补充检查使用训练常规CD4 profiles的Activated_T均值q75固定分层，在留出重复评价，不新造细胞类型。Wilson区间反映有限模拟的条件变异；父参考profiles重复使用，不能解释为患者泛化CI。高激活留出样本少，明确有限性。

## C：绘图和缺失语义

原24中心、候选/对照身份及统计半径完全不变；显示边长512 µm，质心累计圆50/200 µm。原始RNA计数来自v1已验证bin_evidence；全基因矩阵不重跑。T32与Epi64/旧cell区域成员按ID对齐，来源footprint和共享排除明确显示。每ROI的来源数量、深度、覆盖及邻居距离另存表，可直接复核。

白底Arial，PNG300dpi，PDF文字/坐标/图例矢量，H&E与密集点层栅格。提供完整、主体、图例PDF；统计图无图例者仅完整版本。灰色supplied背景不是低分，紫色叉为资格不足，色条0是合格计算零，白色无supplied，旧cell缺项明确未测量。截顶用色条上端三角及cap文本；比例尺64 µm。UC_SCC仅hires，MPBC无H&E与原掩膜，不能读单细胞形态或组织面积。

## 复用与依赖

软件版本与实际日志见tables/software_versions.tsv和logs。数值线程1、最多4核亲和性；不修改共享环境。所有新区域结果和来源映射为Parquet/TSV，没有新增全量表达对象。R加载器经Python桥读取轻量区域表，不要求再保存Seurat。旧counts通过manifest路径访问。

方法接口沿用已核验v1实现；本轮不切换到AUCell原生默认值，也不重新跑UCell参数选择。参考官方定义可见[UCell参数](https://bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html)，具体实现以冻结代码和数值核对为准。可复跑入口支持检查点；--force仅重算本轮输出，旧结果保持只读。replot.sh从表及已有图像重绘，不运行模型。
