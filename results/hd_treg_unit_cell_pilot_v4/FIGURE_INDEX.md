# 图册读取

ROI_atlas.pdf 按 roi_manifest 顺序排列，Cervilla 四页、D1 六页。各页使用同一片的原配准，左上为 H&E 与核边界，x 向右、y 向下，黑色比例尺为 50 μm。第二面板叠加相同核轮廓。所有坐标均为 μm。

从左到右，第一行给出 H&E、2 μm T 面板与 FOXP3、2 μm IL2RA 与 CTLA4、8 μm 的 T 面板与 FOXP3。蓝色为前一信号、红菱形为后一信号、紫色为同一捕获格或父格内共检出。浅灰格为有来源的位置，白色可查看区域超出供应范围或格间间隔，不能读成阴性。逐格计数、原有排除、未归属状态与同核/邻核关系可由 sources.parquet 与 FOXP3_capture_localization.tsv 查回。

接下来的五幅图表示五种单元。紫色为模型首选 Treg，绿色为达到参考支持等级，橙菱形为原始 FOXP3 阳性，橙圈表示双重筛查命中。浅灰为其余有计数的单元，深灰小方形为无有效计数的单元。低置信候选保留。原始单元计数见 unit_comparison.tsv，逐对象证据可查每片的 evidence.parquet。

T32 历史面板中，青色为合格 T32，橙菱形为原 FOXP3-supported，紫三角为原 strict。复用的旧细胞模型在 Cervilla 展示精细 Treg 首选，在 D1 展示 BLCA T/NK 首选。MP10 颜色标到原 Epi64 的完整成员格，灰色为有 native8 输入但未进入可评分上皮区域，白处不据此解释上皮状态。Epi64 是一个测量单元，成员着色不增加样本数。独立邻居机会用全片查询，可能位于页面外。

reference_confusion.png/PDF 展示 100 UMI 样本留出原始首选类别的逐行比例。完整四种检测量和支持阈值的分子分母见 confusion.tsv 和 reference_validation.tsv。
