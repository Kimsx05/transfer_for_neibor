# 空间通讯方法与复现

输入详见[input_review/README.md](input_review/README.md)，涵盖counts、Q、坐标、冻结MP10与身份来源。先运行scripts/prepare_inputs.py（实际脚本名见scripts），再依次运行commot_run.py和scc_batch.sh，最后spatial_postprocess.py及spatial_finalize.py。复现脚本run_all.sh固定顺序和环境。

## COMMOT

官方源码0.0.3，commit d117445bc07eaa19109c7609e97d9e35b26e99ca；源文件未改。项目隔离runtime补POT/plotly/narwhals，复用现有科学Python库，仅按命名空间导入官方推断/database模块，避免无关tradeSeq/pysal/karateclub导入。每片完整背景，human CellChatDB、三类蛋白LR。多亚基在过滤和推断均heteromeric=True、min、下划线分隔。每亚基必须实测且在≥5% spots表达。精确ligand/receptor对去重保留一次，所有数据库记录另存。

所有LR一次联合求解，默认四种竞争COT权重各0.25，eps=0.1、rho=10、最多10000迭代、无平滑。距离字典覆盖当前源码会访问的ligand×receptor笛卡尔积；非真实LR不会形成竞争边。官方toy例与小矩阵F/D显式求和均通过。Inhouse为运行效率采用项目并发调度封装：四个完整原生COT组件分别接收相同全LR输入，并按原0.25权重与相同运算次序合并；没有拆分LR集合或修改组件算法。多亚基、多对多、零表达和边界小矩阵的逐LR与总矩阵parity最大误差均0，见commot/parallel_scheduling_parity.tsv。主/子进程整组资源另由日志抽样记录。API未暴露收敛诊断，不把有限/非负检查当作收敛证明。参数和实测资源详见逐片COMPLETED.json。

原T_total与逐LR网络保存；Tcross只为汇总将对角项置零，原求解不变。F=QᵀTcrossQ；D=outer(sumQ,sumQ)−QᵀQ；G=F/D，D=0为NE。D包括跨位置实际零边，不是有边条件均值。Q不是类型特异实测表达，纯软归属无法确定受体究竟来自同spot哪个细胞。

## SpatialCellChat

独立包0.1.0，commit be88f300464cf970b36a6165ea7d2e41a47a7511，隔离Rlib；旧CellChat1.6.1对象未改变。完整设置见[算法说明](spatialcellchat/METHOD_NOTES.md)和parameters.tsv。实测log1pCP10k，无ALRA补值；原生MERINGUE空间特征筛选和原生位置级概率/QC。没有group permutation P，主指标不做group P二次筛选。

每LR跨位置P经drop0清显式零，F=QᵀPQ，K=I(Q>0)ᵀ I(P>0) I(Q>0)。官方computeAvgCommunProb_LR_Visium_Avg给B=F/K，无链接贡献0，不可评估类型NE。总mean_strength=sum_l B，另存sum_l F和逐LR K。这不是每个真实细胞平均。major/fine分别从同一P聚合；固定LR矩阵不会因展示哪个谱系而变化。

## 距离、状态和呈现

使用真实Visium阵列坐标，名义dNN100μm。Secreted/ECM200μm，Contact105μm；1e−6μm仅数值边界容差。首片COMMOT Secreted300而ECM保持200；SpatialCellChat受原生接口限制Secreted+ECM同时300，Contact保持105。两算法逐片共背景/共LR/共状态阈值，绝不拼接切片坐标。

每个major和所有可辨识fine谱系均High/Low/Both独立圈图。显示掩膜只作用于矩阵副本，保留双向箭头，隐藏TME↔TME、Epi↔Epi、群体自环。等大节点；三图目标固定位置、颜色和共同线性映射，所有正边显示，无top-N。不同谱系、片和算法线粗不可直接比较。NE与真实0严格区分；Tip/Cycling EC不存在的空间组成列以NE注明。

两种份额分母分别为同对象共同可计算TME集及同fine谱系集；low/high使用同一集合。排名min-tie，单一子类份额1只是结构事实。所有表同时保留低/高、差、flow、分母、覆盖和身份状态。spot/LR和内部QC不当患者；这两片仅1个来源核验患者+1个未核验记录单位，不进行患者差异检验。

## 官方来源

[SpatialCellChat代码与教程](https://github.com/jinworks/SpatialCellChat/tree/be88f300464cf970b36a6165ea7d2e41a47a7511)、[官方Visium聚合实现](https://github.com/jinworks/SpatialCellChat/blob/be88f300464cf970b36a6165ea7d2e41a47a7511/R/modeling.R)、[COMMOT空间API](https://commot.readthedocs.io/en/latest/api/commot.tl.spatial_communication.html)、[COMMOT数据库API](https://commot.readthedocs.io/en/latest/api/commot.pp.ligand_receptor_database.html)、[COMMOT原论文](https://www.nature.com/articles/s41592-022-01728-4)。
