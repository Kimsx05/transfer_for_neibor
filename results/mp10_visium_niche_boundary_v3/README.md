# 普通 Visium · MP10 的逐片距离曲线与免疫空间结构 v3

**先看逐片散点，再看组间推断。** 63 张切片都在[逐样本索引](SAMPLES.md)中。122 个有观测界面的组合各有原始、中心化的近窗和全域 PDF，共 488 张。四个 ENN+ENS 掩膜为空，索引保留原因。每张图从左侧外部正距离读到右侧内部负距离；灰点是 non-E 混合组织的 signature 分数，阴影为 IQR。

60 张肿瘤片的仅使用免疫空间特征、对 MP10 设盲的 PAM 分为 C1/C2/C3，共 18/28/14 张。簇间的局部 T 信号、进入份额和异质性不同，但标签对网格与来源敏感，没有足够证据将三簇命名为“浸润好／浸润差／mixed”。两套边界的负侧原始曲线存在探索性的组间差异。相对各组自身水平，两种边界均未显示稳健的纯形状差异。

|目的|入口|
|---|---|
|查看每片的两套边界、完整分布与免疫地图|[逐切片索引](SAMPLES.md) · [正式图与图注](FIGURE_INDEX.md)|
|读科学结论、病例与局限|[RESULTS](RESULTS.md) · [METHODS](METHODS.md)|
|查逐片分型、稳定性与来源|[逐片标签](classification/section_classification.tsv) · [PAM medoid 和 margin](classification/PAM_sections.tsv) · [GMM 比较](classification/GMM_models.tsv) · [重采样](classification/bootstrap_stability.tsv) · [尺度变化](classification/sensitivity_sections.tsv.gz)|
|查曲线差异及支持|[四项主检验](tables/four_primary_tests.tsv) · [全部函数检验](tables/curve_tests.tsv) · [单位与来源分母](tables/curve_support.tsv.gz) · [拟合差与点位区间](tables/model_differences.tsv.gz) · [整体与组成比较](tables/auxiliary_tests.tsv.gz)|
|查看轻量数据及计算依据|[点表索引](cache/sections.tsv) · [字段说明](cache/FIELDS.md) · [连续窗口数据](cache/immune_windows.parquet) · [特征词典](audit/feature_dictionary.tsv) · [几何和组别配置](config.json) · [分型锁](classification_lock.json)|
|核对实际输入、运行和验收|[Git blob 与本地 SHA256](audit/actual_inputs_vs_commit.tsv) · [临床文献证据](audit/section_literature.tsv) · [验收](audit/acceptance.json) · [资源](audit/resources.json) · [轻量复算](audit/lightweight_replay.json) · [文件清单](file_manifest.tsv) · [模块完成与支持限制](module_status.tsv)|

## 复现

在本目录运行以下入口。项目 Python 环境由 `run.sh` 选择，可用 `VISIUM_PY` 指定相同依赖的解释器。R 组成图使用本机已安装的 jins-routine-plotting 模板。保留与 v2 的相邻目录布局和只读点表。

```bash
bash run_all.sh
bash replot.sh
bash regroup.sh new_500_T 500 T 20261127
```

第三个命令只使用免疫特征缓存、原分箱和轻量输入，为新版本重新拟合、锁定并比较。支持网格 `500`、`1000`、`x250`、`y250`、`xy250` 及 T/I 轴。输出写入 `versions/新名/`；主版原始逐片图可共享。默认切片分数、P95、high/low 和距离均沿用 v2。全流程无需读取 counts、Seurat 或 H5AD。实验随附逐文件清单，未制作重复的大对象归档。
