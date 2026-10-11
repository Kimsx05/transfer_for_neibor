# K32 lower80 v7：实际方法

本轮按用户附件授权，完成一个新 NB、精简 DEV 和两张真实切片探索。旧 v5/v6 保留只读；此前范围阻塞不是生物学失败。

## 成员、阈值与来源

逐 cell_id 一对一核对 v5 actual_cell_split_manifest 与 v6 member_label_manifest。K32=冻结FULL_Z_r1.00的C18。旧冻结上皮123699、K32 3282不变；18182范围外上皮仍独立排除，非低状态/阴性。合法Cycling来源按原cell_id保留。TRAIN/DEV/TEST 25/8/8个保守来源组保持，患者未独立核验；P07T/P07N在TRAIN。

TRAIN57008个冻结上皮的float64原始ucell10按NumPy quantile(method=linear,p=.8)，等价R type7：0.22016683621566624。K32优先保留全部；非K32且<=t为Epi_MP10_lower80；>t为withheld_MP10high_nonK32，非reference第三列。其余29类原标签不变。新TRAIN126791细胞（2395 K32、45544 lower80），31类，无均衡抽样。106个低于或等于t的K32保持正参考身份。来源支持与全部逐细胞标签见00_manifest。

## Counts、G0与NB

复用已审计真实整数RNA稀疏counts，冻结G0全部11230基因及顺序，SHA256顺序14316088d684099c8695b1427f6691e7e8945a96986d17d3a034f76ad3260954。旧大缓存继承已核验hash并记录本次stat；新输入重新hash。MT处理和G0过滤沿用旧规则，不扩展基因。Inhouse MT比例仍不可可靠评价。原始计数不由归一化矩阵逆推。

cell2location0.1.5/scvi1.3.3/torch2.6环境，单线程数据进程。RegressionModel labels_key=label_v7，batch_key=sample，无目标生物学nuisance。250初始、125延长、最大1000，batch2500、lr.002、train_size1；相邻50epoch loss均值相对变化<.01。signature=posterior per_cluster_mu_fg.T × 本次实际TRAIN细胞加权batch detection_mean_y_e均值。所有31类重估，无历史signature拼接。1000个参数后验用于signature，1000个固定TRAIN细胞作重建QC。有效但未稳定参考标PROVISIONAL_REFERENCE；下游开始后reference hash冻结。

TRAIN表达效应：原始计数均值/检出率，log2FC pseudocount .1，无细胞级p值。来源方向要求两侧>=30。空间展示基因从冻结MP10 Top20、两张实际可测、继承脂质通路注释、TRAIN K32检出>=.1集合中，按TRAIN原始均值log2FC排序，固定ELOVL5/INSIG1/LDLR；不按新空间图选择。MP10差异部分为lower80构建预期。

## DEV设计

436唯一bag、15实际sample、8保守来源组。与预计14的差异如实按精确来源清单保留。每个20个不重复、同一真实sample的DEV细胞；两重复，全部支持来源保留。A：上皮4/12，K32占比0/.25/.5/.75/1，余为lower80；B：对应withheld高分非K32剂量，exactK32=0；C：上皮0/4/12/20全部lower80。共享0剂量去重只映射一次。每sample/repeat固定池排序、背景细胞和嵌套替换，A/B共同可建背景成对比较。不可构建格保存，不借TRAIN/TEST、不跨sample。细胞复用次数保存。

p=min(1,7005/fullRNA bag总UMI)。每细胞G0整数counts二项稀释，非G0 RNA remainder单独二项稀释，数学上保持fullRNA稀释总量分布。所有20个细胞包括withheld均参与。无放大/归一化回7005；达不到目标记录。稀释后RNA零不改变真值。真值分别为冻结K32数、ucell10>t的上皮数、真实上皮UCell均值、全部上皮/总细胞数；31列truth和可以小于20（withheld另列）。N_epi0比例/上皮UCell均值NA。

## 映射、续训与后验

DEV与两张切片各独立完整31类cell2location映射。DEV输入obs仅sample，不输入真值/剂量/场景/p/UCell。真实切片GSE171351_B1 1017spots、GSE238145_B17 2399spots，全原QCspots。两个历史接受model_config均N_cells20/detection_alpha20，lr.002；不按新结果调先验、不筛上皮、不删vCAF。主参考固定G0同序。

初始3000，续训2000，最少5000、最大30000。实际完成epoch由callback记录；每600秒保存model/history/Pyro optimizer get_state，续训由set_state恢复真实Pyro优化器（不是Lightning dummy optimizer），保持模型参数，随机种子固定；每块软上限2100秒。正式256抽样丰度稳定性只比较相隔至少2000的里程碑，K32/lower80/两者和 abs(new-old)/(old+.1)中位数<.05且q95<.2，同时相邻500 loss均值变化<.005。未稳定不称通过。用户后续允许资源足够时并行，最多2映射进程；第二进程启动前须空闲显存>8500MiB且利用率<95%，仅管理本轮PID。

最后1000 w_sf draws优先；中间256draw初版。逐draw先计算总上皮与I_state=K32/总上皮，再汇总mean/median/q05/q95；draw分母<.1无效，>5%无效draw则spot比例汇总NA。E[ratio]和ratio-of-means分别保存。所有数量保留；不把总预测归一至20。技术重建QC使用100draw技术参数的均值代入，非联合预测区间。

## 空间分数和图

固定Top20来自冻结07_EpiMP_top50_genes.csv rank<=20，完整原矩阵两张均20/20；G0仅19/20，二者分别记录。B1复用精确spot/gene-set一致的既有Top20_A_UCell/AUCell；B17在18030完整实测基因上补算，沿用既有run_benchmark.py公式：rankdata(-counts,average)，UCell maxRank1500，AUCell cutoff=ceil(.05*nGenes)。未改Top20、未用G0代替排名宇宙。EPCAM/KRT8及可测KRT18和展示基因用log1p(CP10k)作表达图，缺失不造零。

坐标来自原obsm spatial，方向一致，无嵌入可靠组织图则坐标图。PDF/300dpi PNG透明底，MP10 plasma、真零黑、NA灰。两张 abundance pooled q99共同色标，裁切数单列且完整值保留；I_state0–1。原始MP10是混合转录组分数，与同程序基因一致仅表达支持，非独立验证。不把TRAIN q80阈值套入spot。

## 评价和预算

主权重先每来源组内sample等权，再来源组等权；每指标给有效来源/样本/spot数。A/C unique回收K32/总上皮/比例，误差/偏差/校准/90%区间覆盖；恒定总数20不计算相关/斜率。B的exactK32误差与广义状态响应分开，低状态/纯非上皮分布和vCAF等29类完整分配表保留。无二元阈值调优，无逐spot显著性检验。除非另外实际完成同输入DEV_OLD，旧v6 MAE仅历史背景，不能算本轮标签改善量。

t0=2026-10-10T16:16:54Z（北京时间2026-10-11 00:16:54），训练软截止07:16:54，交付目标08:16:54。单调时钟采集点和墙钟关系存CONFIG。时间限制、数值状态、生物学证据分开。实际epoch/draw数见MODEL_INDEX及各RESULT_STATUS；不以请求值冒充完成值。

K32来自旧全发现集，DEV/TEST此前已查看；本轮为探索性迁移，不是全新独立验证。伪spot不能完整模拟平台差异。当前范围不扩展63张、G1、阈值扫描、NMF/CNV/预后/CellChat/niche。

本轮K32 TRAIN中P07T贡献1461/2395（61.0%）；v7只训练一个新reference，没有新增去P07T拟合，因此来源方向与DEV迁移不能单独排除P07T依赖。继承的旧轮次敏感性结果只能作历史背景。

90%区间覆盖另按K32阳性/零真值及内部比例(0<真实比例<1)分层。连续正丰度后验的q05通常大于0，比例区间通常不包含精确0/1，因此总体覆盖率受边界真值组成影响；不单凭包含大量零值的总覆盖率断言内部区间校准，且不把q05>0当检出。

续训可靠性附加核验：将B17的11000-epoch检查点复制后在CPU重新加载，独立256draws与原GPU256draws比较。K32/lower80/总上皮相对差异中位数0.742%/0.306%/0.314%，q95约2.23%–2.55%，远小于当时跨2000epochs的变化；不同设备随机流不要求逐位相等。安装源码确认on_load的一步初始化先于load_state_dict恢复保存权重，不计作保留参数的新增训练epoch。记录见00_manifest/reload_audit。

## 最终实际运行量

NB 250 epochs；DEV 9000 epochs、1000 draws、STABLE（NUMERICAL_STABLE）；GSE171351_B1 30000 epochs、1000 draws、PROVISIONAL_NOT_CONVERGED（EPOCH_CAP）；GSE238145_B17 29000 epochs、1000 draws、STABLE（NUMERICAL_STABLE）。

未追加旧S31同输入配对基线：附件要求默认三个mapping均已稳定后才启动可选基线，本轮该前提未满足。不是因DEV生物学效果不佳而跳过。旧v6数值仅为历史背景，不计算“标签修改带来的MAE改善”。第三张切片也未追加。
