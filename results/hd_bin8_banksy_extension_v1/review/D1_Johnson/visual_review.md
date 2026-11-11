已实际打开本轮全片角蛋白图、conflict_focus 单基因对照图，并结合上一轮已复核的同坐标 H&E/域图。14个指定基因全部存在，没有以近似基因替换。

补图让主要冲突区域的上皮证据更明确：lambda0.8簇1/2的 KRT14、KRT6A、KRT17、KRT5 检出率分别约79–85%、77–82%、79–83%、82–84%，相对其他 bins 的平均线性标准化表达均为正效应。COL17A1也有支持。相反，EPCAM约3.6–4.4%、KRT7约0.46–0.81%、UPK2<0.03%；只看后一组会低估这里的上皮相关表达。

conflict_focus 中 KRT14/KRT5/KRT6A/KRT17 呈连续而非单个孤立点信号，COL1A1/COL1A2 也广泛检出。需要区分“检出”和“富集”：簇1 COL1A1检出84.75%，但相对其余 bins 的 log2(1+均值)效应约−0.001，几乎没有区域特异富集；簇2 COL1A1为89.35%、效应约0.533。不能把胶原阳性直接等同为成纤维细胞身份。

建议保留原区域标签，新增补充描述“上皮信号明显，同时伴胶原信号”。这更支持原参考冲突存在面板/参考响应差异，不能据此消除实际共检出的胶原，也不能区分 bin 内多细胞贡献与细胞自身 ECM 表达；不推断 EMT、双细胞比例、CAF或恶性。双方一致上皮候选区域作为比较，不视为纯度真值。

原来的两套类型候选均保留；lambda0.8组织域配合lambda0.2局部结构仍适合继续核验。是否可进行细胞类型特异归因仍未解决。此次没有重跑旧片模型、逐binSingleR或RCTD。

冻结模型核验提供了具体参考机制线索：六个扩展上皮相关基因在BLCA test.genes中全部可匹配，但实际选用的732个判别marker中仅KRT17在内，KRT14/KRT6A/KRT6B/KRT5/COL17A1均不在；HPCA则选用了全部六个。四个胶原/间质基因被两模型共同选用。这是有效判别特征选择差异，不是query缺基因或Ensembl未对齐。它可以解释部分角蛋白明确区域仍被BLCA偏向间质参考的现象，但未进行因果消融，不能量化各因素贡献或据此把HPCA全部判为真值。完整检查在`tables/reference_expanded_panel.tsv`；冻结BLCA参考的对应基因表达分布另存，不重新训练。

The frozen BLCA epithelial reference itself detects KRT14/KRT6A/KRT5/COL17A1 in only9.2%/13.0%/19.0%/10.3% of its2,807 epithelial cells, whereas D1 conflict clusters detect these much more broadly. This describes reference/query phenotype mismatch; per-cell and per-bin detection fractions are different units and cannot be treated as an accuracy or direct purity comparison. The same reference also has nonzero keratin detection in mesenchymal classes (e.g. KRT17:31.0% fibroblast,26.0% mural); this warrants caution about marker exclusivity, not relabelling the frozen reference here.
