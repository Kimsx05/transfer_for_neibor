# Cervilla 与 D1_Johnson 固定区域方法

## 取样与坐标

复用 v3 登记的两片 ROI，以半开、256 × 256 μm 的范围纳入格中心或核中心。两个 Cervilla 候选相距约 32 μm，属于同一热点。逐片合计对重叠 ROI 的同一对象去重。捕获阵列原点固定。格中心由 v1 阵列仿射变换和 native_mpp 获得，原始位置以 fullres 像素乘 native_mpp 换算，x 向右、y 向下。模型图采用既有各轴 scale、半像素公式及 D1 crop offset。H&E 使用对应核模型的配准图。没有 supplied 格的地方保留缺测。

## 计数与历史对照

读取带 20 μm 上下文的 supplied 2 μm 原始数据，并对选中核补取其全片 8 μm 归属的所有来源。每个核内和 4 μm 来源都包含在后者内。8/16 μm 使用阵列索引整除 4/8；上下文覆盖完整父格，不预先做 UMI 过滤。逐条核对 barcode 和 UMI，逐基因核对每个分支的纳入与排除账目。使用既有互斥 nucleus_label、v1_cell_label、cell_label。4 和 8 μm 的表达逐元素核对并读取历史 h5ad。逐核查看范围扩张时各面板的计数。FOXP3 捕获格表给出局部 12 μm 中同核和邻核的 T 支持，12 μm 仅为查看来源的小窗口。既有核分割保持固定，本轮没有独立的漏核或粘连真值。

native8、T16、T32 的旧标签和细胞原 SingleR 从历史文件直接读取。T32/T16 先聚合表达再评定。historical_native8_context 可用 parent8 接回成员表。Cervilla 的旧局部 RCTD 仅在原拟合可匹配 bin 上有结果。D1 旧细胞 BLCA 为大类分辨率。FOXP3_raw_screen 表示原始 FOXP3 非零，FOXP3_T_screen 另要求同一新单元内有 T 面板 RNA。v3 的 FOXP3-supported 和 strict 仍属于原 T32。

## 分类与测量资格

复用 v1 冻结的参考清单和 sample_holdout，样本为 orig.ident。参考含 30 个细类。按共同 symbol 对齐，HD 同名基因先合并原始 counts，再用共同基因文库总量归一至 10,000 并 log1p。仅在训练折内筛选每类 24 个相较全训练均值高的表达特征，再取并集，排除 Top20 及六个通讯读数。L2 多项逻辑回归，以训练样本分组三折、100 UMI 抽稀的 macro recall 比较 C=0.1、1。用同一训练折的 out-of-fold 结果确定支持阈值，见 method_config.yaml。保留全部原始 Treg 模型候选和 FOXP3 阳性，无可用计数的预测记为空。25 UMI 为描述界线。进一步证据要求预测 Treg、有原始 T 分子并达到训练阈值。该等级仍需留出精度与逐核来源共同解释。抽稀以二项采样获得目标预期 UMI，留出集同时保留原深度。混合测试复用 v2 的冻结父对象和既定比例，以 500 UMI 生成本模型的有限表达实现。训练与留出按样本分隔。所有指定的非 Treg 留出对象进入 FPR 分母。

## MP10 与来源排除

直接引用 v1 全片 Epi64 区域和成员，v3 合格区域、坐标、分数、H05/H10、阈值逐 ID 对照。背景 18,085 基因、冻结 Top20、平均并列秩、K=905 均为原计算配置。本轮不重算。图上对成员着色只表示相同区域的实际支持范围。邻域统计按独立 region_id 计数。50/200 μm 均为累计物理半径，同片同原输入组件，采用单位质心与全片 Epi64 区域质心之间的欧氏距离。所选对象任一来源 2 μm 格的父格属于 Epi64 时整对排除。组件或父格不明时独立机会记缺失。可观察支持字段为圆内已 supplied 的同组件 8 μm 格数，不能推断未测区表达。未进行全片富集推断。

## 方法依据

[Seurat HD bin 工作流](https://satijalab.org/seurat/articles/visiumhd_analysis_vignette)和[细胞分割工作流](https://satijalab.org/seurat/articles/visiumhd_analysis_cell_segmentations)用于区分固定格与核支持对象。[Visium HD 原始研究](https://doi.org/10.1038/s41588-025-02193-3)提供小 T 细胞的核级表达定位依据。[膀胱癌 T 细胞研究](https://doi.org/10.1016/j.cell.2020.05.017)支持对照相近 CD4 状态。本轮阈值与训练策略是项目内的受控比较定义。
