# 两片候选联合域的重建与空间比较

## 输入与空间单位

两片分别处理。input_manifest.tsv 记录实际 2 µm filtered H5、positions、H&E、scalefactors、核掩膜及输入审计。使用已有 fullres 坐标乘本片 native_mpp 转为 µm，x 向右，y 向下。D1 使用原全图像素坐标，读取捕获区图像时先减 crop 原点 (23764,35851)。图像重采样的半像素偏移、两片各自的阵列仿射与旋转都来自原 transform/input_audit。2 µm 格反映 RNA 被捕获的位置。没有 supplied 格的图像范围缺少观测。

## 候选与核心

全片读取每个 supplied 格在 T 面板和 11 基因程序上的原始整数，并保留未归属格。历史来源包括两套参考的 T/NK、Cervilla 精细 T 标签、v4 核/归属/网格结果、T16/T32 原区域与原生 8 µm 注释。对两片完整旧 cell 矩阵重用了 v4 冻结多类别分类器，保留 T 类第一预测，包括低支持。所有重建细胞有 T 或程序计数或第一预测为 T 的均加入。不同证据单元具有不同前缀。同一对象的规则合并，区域与格通过实际成员关联一个或多个新细胞。原始 RNA 及历史预测数目不能相加解释为细胞数。原分类器的归一化、特征和排除基因均由冻结模型读入。

所有候选半径 200 µm 圆的并集定义分析域 D。细胞是否在 D 内由候选树到核中心的真实距离决定。ROI 编号用 8 µm 栅格对圆并集取八连通分量，仅用于分组与绘图。ROI_manifest.tsv 的坐标列是候选中心的包围框；完整域见 D_geometry.npz。计算整片 supplied 网格，任何候选的 200 µm 查询均可跨计算块，真实边缘仍限制可观察性。candidate_coverage.parquet 记录在最近 2 µm 格心处的圆形网格覆盖审计，中心量化误差最多约 1.42 µm。正式中心机会基于原 v1 同组件、supplied 且非明确组织外的 8 µm 格心覆盖，圆内格数乘 64 除以圆面积，最高取 1，要求至少 0.5。两种分辨率对应不同的覆盖账目。

## 核与 RNA 的归属

沿用两片均经过原图配准的 StarDist 2D_versatile_he 核，prob_thresh 0.1，nms_thresh 0.3。上游曾比较 Cellpose-SAM 并报告边缘小核遗漏。本轮没有再次推理或声称修复所有漏核。Proseg 3.2.0 源码提供 Visium HD 入口。本地编译器 rustc 1.75 低于要求的 1.88，隔离工具链下载未成功。真实安装日志在 logs/proseg*。这里实际使用核种子驱动的替代归属，没有生成 Proseg 后验。参考 [Proseg](https://github.com/dcjones/proseg) 与 [Cellpose](https://cellpose.readthedocs.io/en/latest/settings.html)。

在有观测且可用于归属的 2 µm 格中固定核内标签。在八连通网格上传播核标签，每一步的成本为物理长度乘以 (1 + 2 × 三通道平均绝对色差 / 255 + 0.15 × 两格 log1p 总 UMI 的绝对差)。总成本上限 10，平局取较小核 ID。每格最多归属一核，未分配计入账目。核外范围为推定边界。本方法只用总 UMI 和 H&E，不依赖基因身份。512 格计算块各含 8 格 (16 µm) 重叠，超过最大传播距离。每片一个历史 ROI 对照局部计算与整片拼接后的所有格。原先规划的外层 64 µm 缓冲通过整片计算覆盖，未额外加入分析域。核内固定与逐基因守恒均已执行核验。负对象 ID 在放大图索引中表示未归属捕获格，编码为 -(source_bin_index+1)。

## 身份与评分

实际预测来自 v4 的 30 类 L2 逻辑回归，C=0.1。上皮、CD4、Tfh/CD8、髓系、内皮、成纤维等共同竞争。使用模型冻结的共同 symbol，先相加同名 counts，再进行总量一万归一化和 log1p。判别特征排除 MP10 Top20 和预定通讯读数。原始 RNA、第一/第二预测、支持与分差分别保留。具体原始面板和身份条件见 frozen_config.yaml 与 reconstruct.py。模型分数仅用于参考支持。推定混合依据 T 谱系与上皮/胶原/髓系面板同时成组检出。Cervilla 的 ROI 2–4 来源保留，几何复核发现碎片或孤立点，排除正式身份。两片在这一步之前满足主 Treg RNA 与模型条件的对象已为零。

从新 counts 直接调用历史 selected_ranks 和 mp10。平均并列秩包括实测零，18,085 基因为背景。得分为 Top20 各基因 max(905-rank,0) 之和除以 rank=1…20 时的同一分子。按每片 8 个历史 Epi64 的完整原生 8 µm 成员重算，核对结果见 historical_score_check.tsv。新评分对象须在 D、身份合格且 nFeature≥50。保留合格计算零。H05/H10 为本片 D 内可评分上皮严格大于线性 q95/q90。以所有 ROI 为分母。

## 效应与读取

各片各自计算核中心 50 和 200 µm 累计欧氏圆。每个上皮仅计一次，每个中心同时记录全部合格上皮、可评分数、high 数与均值、覆盖和成员 ID。模型匹配限制为同组件、400 µm 区块、上皮机会和 T 侧 UMI 层。块是局部组织背景代理，没有独立病理区室标注。比较病例与匹配对照均值，先在块内聚合，再等权平均有匹配的块。至少 5 块才报告区块 bootstrap 区间。候选级分支单列。

按切片读取 cells_scored.parquet、neighborhoods.parquet、candidate_ledger.parquet 和 candidate_cell_links.parquet。counts.npz 的第 i 行对应 id=i+1；列与 genes.tsv 同序。assignments.parquet 按 source_bin_index 接回原始 H5。无分数的行保留 NA。运行入口 rerun.sh，表读取函数 code/load_result.py。数值线程默认 1，KD-tree 并行查询最多 4，进程限制在最多 16 个可用 CPU。未使用 GPU。详细软件版本见 tables/versions.tsv。

## 图像审阅记录

8 个联合 ROI 先在隐藏 MP10 的情况下记录焦点核、主要竞争归属及组织边缘，再查看评分图。补充放大拼版的查看持续到评分后，未用分数修改边界或身份。原图页和显示拼版的对应关系保存在 tables/detail_review.tsv。每个放大视野由 zoom_index.parquet 的对象或捕获格 ID、file 和 slot 定位。拼版层面的 AI 辅助核查与按冻结规则生成的逐对象身份处置分别保存。0.5 µm/像素为模型显示网格，原始图像精度见 input_manifest.tsv 的 um_per_fullres_pixel。MP10 表中 positive_genes_avg_rank_lt_905 为平均秩在有效窗口内的实测非零基因数。原始 2 µm 覆盖汇总按来源台账条目计算，同一细胞可有不同类型的筛查线索。
