# 独立只读实现审阅

审阅范围：用户本轮附件 B/C、scripts/commot_loader.py、commot_run.py、spatial_postprocess.py、scc_run.R，以及实际vendor官方COMMOT/SpatialCellChat源码。未运行通讯推断，未修改被审脚本。本记录为首轮审阅，后续修订是否关闭应由主进程记录。

## 需要处理的两项

1. **距离边界浮点漏边（实质数值问题）**。直接由导出coords计算median dNN再取2倍，数值可为199.9999999999991 μm。对GSE171351_A1全部位置距离矩阵，几何距离200 μm的4112条有4000条被 `d>2*dNN` 判出；Inhouse对应21790条中14266条。即使固定200，坐标运算微误差仍使部分边略大于200。COMMOT原生cutoff和SCC computeCellDistance均使用直接 `<=`。建议使用round(dNN,8)名义值并在阈值加1e-6 μm纯数值容差，保持预声明生物尺度不变；同步主网、宽范围敏感性与参数记录。如果已完成旧推断，需重推断或显式标注旧边界结果，不能仅改文案。已即时告知主进程与SCC worker。
2. **空间Epi图注含义**。spatial_postprocess首次审阅版本调用render_triplet未传epi_caption，默认“MP10-high/low epithelial states”省略了“关联spot中的上皮成分”。renderer支持epi_caption，空间调用需明确Epi_H/L定义及共享TME，并保持Visium混合表达和c2l依赖说明。仅涉及重绘，无需推断重跑。已通知主进程。

## 已核查符合要求

- 每片完整背景分别推断，同一片H/L共享网络、Q与LR集合，不拼切片坐标、不按H/L独立求解。
- measured counts经全基因library1e4归一化log1p；不使用Q、MP10或Harmony作LR表达。现成log1p文件虽未复用，但从counts分别只做一次归一化，未双重log。
- COMMOT显式human，三类蛋白LR，复合物在filter与inference两阶段均heteromeric=True/min，检查所有亚基测量与5%位置检测。真实LR去重后同一批675等LR联合COT求解，非按通路分开求和。距离字典覆盖ligand×receptor索引域。
- COMMOT位置总网络从官方total-total读取，并与全部单LR之和核对。先原生求解，保存同spot贡献，聚合副本去对角。F=Q' T_cross Q，D=outer(sumQ,sumQ)-Q'Q，G=F/D；D含全部跨位置零值，不是仅正边平均。major/fine分别从同一T_cross聚合，非fine均值相加。
- SCC官方独立包与commit记录；无升级旧1.6.1环境。MERINGUE、位置quantile筛选、min.links/min.cells.sr分层记录；未伪装患者差异P。`filterProbability`源码是采样源位置附近概率的95%分位过滤，当前文案已准确区分。
- SCC保留完整组分Q，不使用dominant硬标签。仅一个whole_slice初始化ident；min.cells=NULL，避免按dominant计数删除Treg。全片所有正式Q列已实际通过原c2l支持规则；postprocess再按type_support屏蔽不可评估者。
- SCC缺测复杂体亚基由官方identifyOverExpressedInteractions中all-subunits条件排除。传LR.use使原生函数按类别排序，避免contact信号索引错位。
- SCC两处probability存储在位置过滤前同步，drop0清除显式零；最终每LR去同位置对角后调用官方Avg helper，且逐LR与F/K核对。K是Q>0支持的正P位置对数；无group P二次筛选。真实零保留0，类型不可评估才NE。
- 输出边source/target方向与矩阵行/列一致：R列主序rep规则正确，Pythonpivot与row/column读取正确；没有把反向加总为有向单值。
- circles仅Epi↔TME，排除EpiOther、H↔L、TME↔TME和群体自环；三图共尺度、固定坐标和target顺序、等大节点、线性权重、没有top-N隐藏。
- 高低份额使用同一可计算目标集合，allTME排除全部上皮列；谱系内来自fine自身矩阵。原Mural独立保留，T_NK_unresolved保留质量；单节点谱系份额1需报告其结构含义。
- 所有独立Q列都保留，包含从未占dominant的GSE A1 Treg；三个原低可辨识列没有被命名为独立生物结论，Tip/Cycling EC缺列明确NE。
- 现有身份字段保留：A1来源核验患者、Inhouse匿名记录单位，不应以两片宣称两个已核验患者。脚本未生成spot/LR伪重复差异P。
- SCC sum-per-LR mean与COMMOT G单位不同，脚本按method分开记录；后续报告只可比较内部方向/排序，不能绝对值相减或合平均。

## 后续报告注意

Treg/多数EC和CAF在已有c2l审查中为INTERPRET_WITH_CAUTION，结果需保留。SCC无置换group P，COMMOT未运行group permutation；“总通讯”不等同免疫抑制。两算法共用组成和相关数据库不是独立实验验证。名义阵列μm坐标并非显微镜μm/pixel恢复。原模型同spot对角剔除只是聚合选择，不宣称禁止同spot求解，也不声称删除真实自分泌。
