# 运行环境与恢复
起始核查记录（2026-11-27，Asia/Shanghai）：主机可见 256 个逻辑 CPU，约 1.8 TiB RAM、1.3 TiB 可用，项目文件系统约 11 TiB 可用。Quadro RTX 6000 显存 23,040 MiB、驱动 580.159.04，GPU 已有其他计算进程，占用约 340 MiB。所有正式阶段使用 CPU，没有申请或干预其他进程。
已有 bin2cell/CellTypist/scVI 及 CellCharter 组合环境可用。TACCO 等仅装在本轮 shared/deps。StarDist 2D_versatile_fluo 权重使用官方 release，SHA256 4ad678d0758eed6e55625f1b5ae30771e59adb79f1239e09b9772eac8846c3dd。TensorFlow 首次推理因 private threadpool 创建失败，改用同一 Keras 网络直接 inference，随后固定 ROI 与全片均完成。TACCO 两处 pandas 3 只读数组错误在本地工作拷贝修复，保留原始错误日志。
主要计算均固定在同一 CPU 亲和集合 2–5，BLAS、OpenMP、Numba 和数据加载线程限制见 code/run.sh。Torch scVI 最大 3 个算术线程，数据加载 worker 为 0。单阶段 cgroup 限 17 tasks、275 GB。后续共享临时用户 slice 设 20 tasks、285 GB；没有改动系统全局设置。logs/resources.tsv 每行记录命令、秒、最大 RSS（KiB）及状态。共享 slice memory.peak 包含文件缓存；初段使用单阶段 scope，无同口径全程合并峰值。
HD 首次全量 scVI 训练在 2026-11-27 01:00 左右退出 143，未写成模型。当时 slice 无 OOM 记录；未将该次尝试计为完成。最终用相同架构与固定随机种子每片均衡抽取 35,000 个合格原始 counts 对象训练，对全部合格对象编码。输入特征按两片全量表达选择。正式结果以 hd/scvi_model 和 hd/scvi_balanced_fit_ids.tsv.gz 为准。Xenium 五片均用全量基线合格 counts 训练。共同 niche 均按完整切片坐标构图，抽样只用于模型拟合。
Xenium 分类器在检查历史参考抽样框时剔除了只为状态验证/转移抽取的细胞；最终训练和留出数据均限 identity_ref，使用 ct_identity_* 模型。Xenium 表达嵌入和共同 niche 未使用这些参考细胞或身份标签。

