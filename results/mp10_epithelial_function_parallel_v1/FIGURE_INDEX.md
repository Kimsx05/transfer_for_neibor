# 图件与绘图数据

每个编号均输出独立 PDF 与白底 300 dpi PNG。效应图和连续药敏散点在 jins-routine-plotting 中没有直接对应的模板，本轮使用该技能的 Arial 主题和 jin_export。区间为 95% CI。

| 图 | 内容与读法 | 绘图数据 |
|---|---|---|
| F01_project_modules | 24 个标本 H05/L50 配对。去 MP10 重叠的程序均值差，右侧代表 High 较高。去重叠后 MUFA 无剩余成员，见覆盖表。 | tables/project_statistics.tsv |
| F02_project_footprints | 同 24 对的 MAPK、PI3K 有符号转录足迹差，去重叠；单独标尺。 | tables/project_statistics.tsv |
| F03_erastin、F03_ML162、F03_1S_3R_RSL_3、F03_ML210 | 各药一图。横轴为原配套基线 MP10，纵轴为实测 AUC。点为合并同患者后的 BLCA 模型单位。 | plot_data/CTRP_units.tsv |
| F03_GPX4_2879_ | 基因依赖单图，27 个 BLCA 单位；纵轴越负表示基线依赖越强。 | plot_data/DepMap_units.tsv |
| F04_induction_rescue | 每条骨髓瘤细胞系 3 个 passage 配对。RSL3−UT 与 FR−RSL3；FR 指 Fer-1 预处理加 RSL3。 | tables/GSE182638_modules.tsv，plot_data/GSE182638_profiles.tsv |
| F05_CHMP4B_localization | SLO 后膜 CHMP4B 相对信号；两种肺癌模型，三种条件。点是源表值，短横线为均值。 | plot_data/HMGCS1_source.tsv |
| F06_sterol_pairs | 去除 MP10 成员后的甾醇程序。每根线连接同一标本的 Low 与 High。 | plot_data/F06_sterol_pairs.tsv |
| F07_model_associations | 四药和 GPX4 gene effect 并列。五行均为 BLCA 的 Spearman ρ 和 bootstrap CI；正的 gene effect 相关表示依赖减弱。 | plot_data/F07_model_associations.tsv |

文件位于 [figures](figures/)；统计量、真实单位及模型标识见 [tables](tables/) 和 [plot_data](plot_data/)。
