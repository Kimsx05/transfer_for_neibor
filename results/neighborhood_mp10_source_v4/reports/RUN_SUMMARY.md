# 本轮实际结果

固定九群；共同 anchors 129,632，正式接受 126,075。中央冻结 MP10 覆盖 61 张有共同 anchors 切片。完整20基因 matched Epi 129,412 anchors（含正式接受标签的125,855个）；另有 B5 220 个19基因 anchors 单列。

七数据集技术 pilot 全部通过，Epi 最大绝对误差 9.09e-13；冻结分数复核 61 切片，最大差 9.89e-17。

来源：9600 真实单细胞，360 蓝本×5=1800 composition-only；99 成对状态蓝本。

病理检索覆盖7数据集/63切片；可直接应用的区域标签0，133 ROI 仍需复核。


|niche|完整整体N|整体中位数|匹配Epi中位数|同中心整体效应|Epi效应|独立单位|来源级别|
|---|---:|---:|---:|---:|---:|---:|---|
|niche0|43871|0.0331|0.0336|+0.0002|-0.0009|44|UNRESOLVED|
|niche1|28250|0.0465|0.0459|+0.0062|+0.0036|32|UNRESOLVED|
|niche2|19089|0.0299|0.0340|-0.0050|-0.0027|50|COMPOSITION_CONFOUNDED_OR_ATTENUATED|
|niche3|14303|0.0361|0.0420|-0.0071|-0.0008|36|COMPOSITION_CONFOUNDED_OR_ATTENUATED|
|niche4|4992|0.0540|0.0528|-0.0028|+0.0002|15|UNRESOLVED|
|niche5|4102|0.0265|0.0360|-0.0057|+0.0026|32|COMPOSITION_CONFOUNDED_OR_ATTENUATED|
|niche6|5272|0.0350|0.0380|-0.0025|-0.0041|25|UNRESOLVED|
|niche7|3819|0.0397|0.0411|+0.0032|+0.0004|30|UNRESOLVED|
|niche8|2157|0.0335|0.0403|-0.0073|-0.0024|21|COMPOSITION_CONFOUNDED_OR_ATTENUATED|


整体列为原冻结背景，Epi列为模型背景，两者中位数不能直接相减。两效应列均来自同中心、同模型背景。


逐切片状态见 ../03_epi_attribution/model_and_gene_universe_manifest.tsv。无共同 anchors 的切片：GSE171351_D1, GSE238145_B13


入口：[双层图](../06_figures/ALL_DATASETS/matched_layers.png)、[来源证据](../05_source_validation/source_evidence_all9.tsv)、[病理复核](../01_pathology_search/user_review/review_index.html)、[限制](WHAT_IS_STILL_UNCERTAIN.md)。


续跑：bash code/resume_mp10_source_v4.sh（在本结果目录内；已完成阶段校验后跳过）。