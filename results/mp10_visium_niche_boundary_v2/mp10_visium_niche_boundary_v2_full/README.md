# 普通 Visium 的 MP10 与上皮相关生态区边界 · v2

完成了 7 个来源全部 63 张切片，含 60 张肿瘤来源、3 张邻近正常。133,064 个 E 富集混合 spot 都有可用原分数。两类边界近侧的总体原始 MP10 均低于深部，调整上皮组成和测序深度后区间跨零。免疫格局呈连续变化，临时规则只给 19 张肿瘤片分配了方向明确的标签。组间比较的精度受来源重叠和可分类切片数限制。

|要查的问题|文件|
|---|---|
|六个问题的结果、效应与限制|[RESULTS](RESULTS.md) · [方法](METHODS.md)|
|单图及每张切片的地图|[正式图](FIGURE_INDEX.md) · [63 张审图](review/INDEX.md) · [全指标按来源查看](metric_explorer.html)|
|逐样本分型、选择依据、置信度、重复标本|[逐片名单和审图](classification/SAMPLES.md) · [标签与不确定原因](classification/section_classification.tsv) · [阈值](provisional_classification.json) · [免疫盲态模型](classification/model_selection.tsv) · [分类输入和标签校验](classification/classification_lock.json)|
|输入、临床信息与评分覆盖|[manifest](audit/input_manifest.tsv) · [逐片原文映射](audit/section_clinical.tsv) · [研究来源](audit/source_studies.tsv) · [MP10 审计](audit/score_coverage_audit.tsv) · [逐片真实交叉表](audit/score_cross_section.tsv)|
|观察到的边界、低覆盖、缺界面及缺深部|[逐片逐边界状态](module_status.tsv) · [所有切片和敏感性效应](tables/section_effects.tsv.gz) · [距离方法比较](tables/border_distance_agreement.tsv.gz) · [成分内效应](tables/component_contrasts.tsv.gz)|
|可重用数据|[逐点坐标、分数、组成和双距离](spots/) · [边段及两侧区域](boundary/) · [局部窗口、组织块和 E 分量](tables/window_component_metrics_with_state.parquet) · [分量外围的独占归属](tables/component_unique_periphery.tsv.gz)|
|连续指标、统计与分母|[指标值及分子分母](tables/section_metrics_with_states.tsv.gz) · [指标词典](audit/metric_dictionary_with_states.tsv) · [全部关联](tables/all_metric_associations.tsv.gz) · [主检验](tables/main_six_tests.tsv) · [共同支持](tables/paired_summary.tsv)|
|整体 MP10、暴露和组成|[所有切片](tables/all_section_levels.tsv) · [high/low 与并列](tables/unit_exposure_groups.tsv) · [区域和 niche](tables/region_niche.tsv.gz) · [组成比较](tables/level_composition_tests.tsv.gz) · [high/low 比较](tables/highlow_levels_composition.tsv.gz)|
|运行与校验|[参数](config.json) · [软件](audit/software.tsv) · [资源](audit/resources.json) · [文件大小及 SHA256](file_manifest.tsv) · [验收](audit/delivery_validation.json) · [完整包校验](audit/archive_validation.json)|

## 重绘、重分组与全量复算

在本目录执行。Python 默认使用项目已验证的 spatial 环境，其他机器可设置 `VISIUM_PY` 为已安装相同依赖的解释器。

```bash
bash replot.sh
cp provisional_classification.json /tmp/mp10_new_rules.json
# 编辑该 JSON 中待讨论的阈值，使用一个新版本名
bash regroup.sh revised_rules /tmp/mp10_new_rules.json
```

重绘只读轻量结果。重分组复用冻结免疫窗口、邻接、聚类坐标、片内效应和评分。输出单独写入 `versions/revised_rules/`。审图优先使用 manifest 的原 H&E，离开原服务器时原图可能缺失，已交付的 63 张审图保留 H&E。R 组成图使用项目已安装的 jins-routine-plotting 模板，其余图有 Python 入口。已用 mixed 份额 0.25 在独立目录实测重分组、重统计和重绘。63 张切片的结果与预置严格敏感性一致。[实测记录](audit/threshold_rerun_test.json)

在含原输入的项目主机，从本目录运行 `bash run_all.sh`。入口核对冻结轻量输入哈希，逐阶段运行，按切片复用已有基础结果。完整重建需要 manifest 中的 H5AD、c2l、v1 及上游脚本。在不同输入上重做时应使用新结果版本。详情见各脚本及 [原始运行日志](logs/)。本目录与完整压缩包均采用项目既有的本地交付方式。
