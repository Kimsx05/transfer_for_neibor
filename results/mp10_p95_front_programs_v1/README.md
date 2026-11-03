# MP10 P95 / 前沿 / 程序 v1

本轮以 MP10 为中心，复用 full_v1 的63张切片与正式v3.4，补充文献程序、真实图邻域、P95参考集合及前沿位置匹配。正常来源单独保留。

- `overview.pdf`：五组主图和预先确定的代表切片。
- `summary.md`：研究问题、客观结果和限制。
- `methods.md`：固定定义、推断和分母。
- `full_sections.tsv`：逐片覆盖。
- `figure_index.tsv`、`atlas/`：63片空间叠加与数据集图册。
- `tables/summary_effects.tsv`：section/unit/verified三级结果、CI、FDR和方向。
- `tables/dataset_effects.tsv.gz`、`tables/unit_effects.tsv.gz`：可核查的重复单位效应。
- `tables/effects.tsv.gz`：全部逐片检验(包括不可估计原因及空间块数)。
- `tables/curve_summary.tsv.gz`、`ring_summary.tsv.gz`：实际距离/图环曲线及覆盖，组织与Epi分别保存。
- `tables/matching.tsv.gz`、`matching_balance.tsv`：逐对匹配及质量，不隐去不理想匹配。
- `tables/niche.tsv.gz`、`bands.tsv.gz`：九类niche和P95分母。
- `tables/gene_sets.tsv`、`signature_index.tsv`、`new_gene_set_coverage.tsv`、`literature_symbol_coverage.tsv`：完整成员、映射和检测。
- `tables/signature_detection_summary.tsv`：短签名稀疏性与单基因计数贡献(非UCell贡献分解)。
- `tables/sc_sample_scores.tsv`、`sc_gene_source.tsv`：已有单细胞注释内来源支持，不是空间配对验证。
- `spots/`：必要几何、冻结分数、P95/活跃前沿/niche和80项主/敏感性程序；`centers/`：R0–R10程序均分与邻居数。
- `scripts/run_all.sh`：从项目根目录执行；已有评分与完成切片可复用，不覆盖上游。
- `review.zip`：报告、主图、核心小表；不含counts或大型中心缓存。

`main`为57项(冻结55+Puram30+Flores10)；23项敏感性为公开Puram100和22项去MP10重叠。NA表示无覆盖或不可估计，不表示零功能。研究来源的肿瘤样本不等于逐spot病理恶性确认。
