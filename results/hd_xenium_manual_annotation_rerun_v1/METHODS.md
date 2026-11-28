# 计算方法

## 输入和 QC

本轮独立实现表达流程。HD 读取既有 bin2cell 聚合的稀疏整数 X、原对象 ID、来源和几何标记，并以 ID 对齐 GEX secondary 的既有影像审阅标记。逐对象复算 UMI、检出 feature 数，与源对象保存的两列逐一核对。HD 保留 UMI ≥10、检出 feature ≥10、基础几何正常且未标记 secondary 影像异常的对象。对象以重建单位计数。

Xenium 从各样本原始 10x H5 读取 barcode、feature ID、gene symbol 和整数矩阵，原细胞多边形、质心/面积按 cell ID 对齐。仅 Gene Expression 进入库大小，逐细胞总数与原始 transcript_counts 核对。各类控制 feature 和计数汇总见 input_manifest.tsv。纳入检出 feature ≥10、多边形至少四个顶点记录、面积严格介于同片原始对象的 0.5% 和 99.5% 分位之间的细胞。各过滤规则记录全部命中原因。旧身份与状态均不参与筛选。

每个平台先检查 feature ID 唯一，跨样本按真实 ID 对齐并核实 symbol。保留重复 symbol 对应的各独立 feature。marker_availability.tsv 区分未测与已测基因，各样本实际的 feature 集一致后才进行 pooled 统计。

## 表达空间

各平台以全部合格对象、全部 GEX 的库大小计算 CP10k、log1p；根据这些值的稀疏矩估计方差，按 variance/(0.1+mean) 选择 1,200 个训练特征。该项目特有的 dispersion 计算使用归一化对数表达，scVI 的似然输入使用所选列的原始计数。样本为 batch，随机种子 1729，负二项似然、每基因 dispersion，1 个隐藏层、64 个隐藏单元、latent 12 维。HD 固定抽取每片最多 35,000 个合格对象组成训练/验证框架，再投影全部合格对象；Xenium 使用全部合格对象。新拟合的模型、各 loss 轨迹及训练/验证成员均保存。训练最多 24 epoch，90% 训练、10% 验证，batch 512，每 epoch 检查验证值、early stopping patience 5。与 retry_v1 的差别包括 HD 先落实最终影像 QC 再选训练对象、Torch 线程调整至 2，两个平台都从新初始化训练。

在新 12 维 latent 上分块求精确欧氏近邻。UMAP 的 30 个输入位置含自身，随后构建默认局部连通性 1 的对称 fuzzy union 图，并去掉对角线。Leiden 在同一无向加权高维表达图上使用 RBConfigurationVertexPartition，三个固定分辨率均迭代至实现检测不到进一步改进，参数 n_iterations=-1。每个分辨率完成即写标签并按大小重编号。UMAP 使用 umap-learn 的 simplicial_set_embedding，在该图上以 PCA 初始化后优化二维坐标 200 epoch，min_dist=0.3、spread=1、负采样率 5、初始学习率 1、固定种子 1729。无空间邻接或既有标签输入。样本单图保留联合 UMAP 坐标与坐标范围。

## 全基因统计

每平台、每分辨率、每簇与本平台其余合格对象比较。对全部实测 feature 的 log1p(CP10k) 进行双侧 Mann–Whitney 渐近检验。零值合并为一组并列秩，正值按实际表达并列求平均秩；U 以簇内为第一组，方差使用全体 tie correction，无连续性校正。全零、常数或比较任一侧少于两个对象保留对应状态及 NA P/FDR。每簇以全部可检验 feature 在 log P 域做 BH。极小的 P/FDR 在浮点格式中可能为零，另给负对数值。预先不限制检测率、效应或显著性。每分辨率另用 SciPy 独立复核一组不同检出率基因的 U 和 P。

mean_norm_in/out 为 CP10k 的算术均值；log2FC 为两个算术均值各加 1e-9 后之比的 log2。检测定义为 raw counts >0，pct 与 delta_pct 均用 0–1 的比例尺度（差值可负）。各样本内使用相同公式，每侧至少 20 个对象，缺少背景保留原因和 NA。样本内实测全零的效应为零，另标记 measured_all_zero。样本方向一致数、比例及效应范围使用可评价样本。样本可以代表切片，患者关系未用于推断。

显示表取 FDR ≤0.05，正 marker 要求 log2FC ≥0.25 且 pct_in ≥0.10，负 marker 要求 log2FC ≤−0.25 且 pct_out ≥0.10。每簇每方向最多 20 项，按 FDR、P、feature ID 排序。图内选择每簇最多前三个正 marker 的并集供浏览，并绘制用户指定的大类和 Treg marker。dotplot 圆面积对应原始检测比例，颜色对应 log1p(mean CP10k)。热图按基因标准化。UMAP 与空间图由 Matplotlib 直接绘制，沿用项目模板的颜色习惯、簇标签和独立单图交付。正式人工字段保持空白。

## 资源与统计解释

HD 与 Xenium 依次工作，每次运行一个主要阶段。systemd cgroup 限制 TasksMax=10、MemoryMax=195,000,000,000 字节，禁用该任务 swap。可运行 CPU 亲和性设为 2，BLAS/OpenMP/Numba/Arrow 等限制为 1 个内部线程，Torch 计算线程 2、interop=1、loader workers=0。使用 CPU。监控进程每 0.4 秒记录整个任务 cgroup 的总线程、全部进程 RSS 之和及 cgroup 内存；内核 pids.peak 与 memory.peak 同时保存，包含短于采样间隔的峰值。logs/resources.tsv 从左到右为阶段、退出码、秒数、采样线程峰值、采样 RSS 合计峰值、阶段采样 cgroup 内存峰值、服务生命周期内核线程峰值和内核内存峰值（字节）。

这些检验使用已由表达聚类产生的分组，适于寻找人工复核的 marker；对象数量、样本组成和稀疏性可影响 P 值。样本方向检查是描述性证据，不能推广为患者级独立验证。

