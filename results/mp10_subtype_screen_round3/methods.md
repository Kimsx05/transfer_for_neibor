# Round3 方法：全 TME 聚合通讯

## 冻结范围与输入

仅使用 v2 的 `comparison_type=epi_state_within_specimen, condition=shared`。五个dataset和ALL，各有major/fine独立推断的一个对象，共12个。保留v2全部上皮、固定20基因评分、标本内相对中位数状态、标本资格与随机抽样结果。没有重新评分、重划状态、重注释、运行CNV、重推CellChat或重跑丰度模型。v1/v2保持只读，新结果都在round3。

本轮实际覆盖40个肿瘤标本、39个记录患者ID；所有ID沿用derived/provisional来源，不是已核实的生物学身份。两个epi状态来自相同标本，逐标本等量抽样；同一个TME细胞没有因两种状态复制。每标本epi每状态上限100、TME每节点上限100，seed123；CellChat1.6.1、human DB、triMean、raw.use=TRUE、population.size=FALSE、nboot100、min.cells30。详细来源和原对象参数保存在object_index.tsv。major来自自己的大群输入，fine来自自己的亚群输入；两分辨率绝对总量不互比。

主要读取普通list型qs中的prob、pval、weight、count、节点/贡献和参数。只为NECTIN2–TIGIT补充读取未抽样expression_fine.qs。没有复制旧网络、表达矩阵或RNA母对象；仅保存新生成的小型聚合矩阵。

## 聚合定义与缺失

在每个对象的完整已计算LR维度上，在每个source×target方向内按真实LR ID各计一次：

- W[i,j] = Σ p[i,j,l] × I(P[i,j,l]<0.05)，主要total interaction strength。
- C[i,j] = Σ I(p[i,j,l]>0且P[i,j,l]<0.05)，支持LR记录数量。
- R[i,j] = Σ p[i,j,l]，不按内部P过滤的辅助原始和。

全局C是跨伙伴细胞群对累计的支持LR记录数，不是全局去重LR种类或独立机制数。

不取两状态共同显著交集。仅一状态支持的LR保留其贡献。核对net$weight/count和v2 node_edges，按本轮明确定义生成W/C/R；当前核对完全一致或只有浮点精度误差，见aggregation_parity.tsv。

计算节点对没有支持LR时W=C=0；节点未达到v2最小细胞数或不在对象中时，扩展表中保留NA/NE。未计算/被筛出的单LR不是测得的生物学零。总量描述该对象已有LR维度，不代表所有可能分子通讯。正式数据从不使用fillna(0)补不存在节点；圈图为固定布局需要的零占位仅抑制绘边，同时保留[NE]标签，表中仍是NA。

## 节点、覆盖与排名

两个epi节点单独置于前方；其余按B、Plasma、T/NK、Myeloid、Fibroblast、Endothelial、Mast分组，组内固定排序。major固定9个展示节点；fine固定47个展示节点，实际可计算数随对象而异。unresolved、dataset_bias及既有其他标签全部保留，包括名称为excluded_myeloid的既有保留节点，不能把名字再次当作新的排除规则。

node_coverage.tsv报告原始数、抽样数、贡献标本/记录ID及可计算状态。每个目标另统计同标本low epi、high epi、该目标各≥30的配对覆盖，仅用于解释，不据此从全网络删除群体。计数只描述当前纳入和抽样覆盖，节点大小不编码生物丰度。

incoming_all为W列和，outgoing_all为W行和，self为对角线；other为all减self。C和R同样分解。求和只跨实际可计算伙伴，NE节点不作为零伙伴参与。目标本身NE则全部总量NA。incoming/outgoing各自排名，只在同对象、同分辨率可计算TME节点中进行，排除两个epi；并列采用min rank，不做任意加权总分。绘图可完整展示epi，但不赋TME排名。含自环与去自环结果同时保留；若相加incoming与outgoing会把自环计两次，本轮不以这种和作为独立事件量。

## 两类epi与全部TME

分别提取epi-low/high→每个TME及每个TME→epi-low/high的W、C、R和high−low差。两个方向独立报告；不取倍数、不加伪计数。low/high强度、差值、份额各自在当前TME集合内排名。

U_out为某epi状态到全部可计算TME的W和，U_in为全部可计算TME到该epi状态的W和；两个状态使用同一TME集合。明确排除epi↔epi和epi自环。份额为对应边W/U，分母0时NA。份额是组成描述，不取代绝对强度，不证明细胞特异偏好。

跨队列只汇总对象内方向、可计算性与匹配覆盖；不平均不同对象概率当作患者效应。五队列作为主要重复信息，ALL另列为pooled概览。W与R方向是否一致也保存。正负差与零差采用实际浮点差的符号，没有根据结果挑选阈值；原值全部可查。

## Treg分解

只从fine对象的单一Treg节点提取全网络incoming/outgoing、去自环、自环、排名及覆盖；没有构造high组Treg和low组Treg两套全局总量。major只提供T/NK背景。

全部伙伴分为epi-low、epi-high、其他TME、Treg自身，四部分之和核对全局总量。对Treg两方向，逐LR贡献严格为p_high I(P_high<0.05)−p_low I(P_low<0.05)；完整贡献和核对W差，原始概率差和另核对R差。解释图可显示最大正/负条目，但同时给“其余正贡献”和“其余负贡献”合计；不能用这些条目反过来筛选主网络。

## NECTIN2–TIGIT补充

按真实DB interaction ID提取同对象Epi→Treg概率、内部P及支持展示值；两边均支持就同时画出。数据库的NECTIN2_TIGIT和TIGIT_NECTIN2作为相关分子对记录，主W/C按原DB标准维度保留，解释时不作为两次独立证据或双向功能反馈。

未抽样expression_fine.qs按同标本两epi状态各≥30、Treg≥30提取NECTIN2均值/检测率及共享Treg TIGIT均值/检测率。两个状态的≥10%检测支持分别记录，不要求两边同时过阈值。基因不在缓存/目标缺失为NA；队列有该基因非零观测证据时保留节点内真实零。该小补充没有重建患者CellChat网络。

先计算标本high−low表达差，重复记录ID先平均标本差，再报告正/负/零方向和leave-one-ID均值范围；无差异显著性P。乘积只叫辅助表达乘积。CD8 Exhausted及其他实际TIGIT阳性T/NK接收群仅给简洁对照表。单轴表达复核不代表Treg全部总通讯完成患者验证。

## 图形与解释

主circle均为全LR聚合后的细胞群网络：A全网络含自环、B双epi→全部TME、C全部TME→双epi，每种分别W/C，全部节点固定位置和颜色、等大圆节点、真实定向箭头。top=1保留所有正边；同对象同指标A/B/C共用边宽最大值。不同队列独立缩放并注明，不能凭边粗横向比较绝对大小。fine保留大画布矢量circle，并以完整source×target热图作为精读对应；灰色是NE，计算零位于色标0端。含自环/去自环IO、全节点角色图、全部TME配对比较、差值与份额、方向矩阵均另给。

正式PDF与300dpi PNG背景透明；少量白底阅读预览放logs。角色散点全部节点显示，密集图采用数字键与完整名称表，Treg高亮。完成后实际查看ALL与Inhouse的major/fine圈图、Treg面板及方向矩阵。

内部P只用于CellChat支持过滤，不是患者状态差P。不同LR、细胞或pooled对象不是独立患者重复；不从现有pooled网络反推患者总强度。总W混合多种联系，不称“总免疫抑制强度”。没有新空间、niche、Treg注释、预后或中介分析，也没有上传GitHub。
