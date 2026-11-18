# MP10 hypothesis follow-up v2

实际结果目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_hypothesis_followup_v2`。代码目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_hypothesis_followup_v2/scripts`；结果目录内 `scripts/` 为同版本轻量交付副本。

本轮以冻结 **H05 vs L50** 为上皮主比较，H10/H50 为敏感性。M0 后连续完成所有可运行模块；旧结果未修改，未发布远程仓库。

## 主要发现

1. 原 SREBF2 六个候选模块均在 ranking 成员映射门槛被排除，并非没有 TF 表达、GRN 边或 motif 支持。PRAG1 的可核实历史符号 SGK223 修正使 top50 从40/51到41/51，通过原门槛，恢复36成员 regulon。H05 活性增加，去 Top20 后仍有支持。
2. 正式集合支持 H05 胆固醇合成；直接与脂肪酰 CoA 程序比较，D=0.679 (95% CI 0.442–0.917)，M20 D=0.808 (0.537–1.079)。这是固定尺度的转录变化，不是通量比；Full Bloch 比较不支持更大变化。
3. 上皮 mTORC1、EC 的 Q50 关联和空间混合 spot 结果来自不同表达单位。EC 固定评分覆盖98/200不足，空间调整后的跨数据集支持减弱。
4. 上皮胆固醇外排程序未获明确支持；两项上皮→Treg 摄取主关联及五项 EC 探索关联均因覆盖不足为 NE。可估计的 cycling 主关联未获支持，不能据此证明脂质转移。

## 阅读顺序

- `RESULTS.md`：逐问题结果、效应与局限。
- `figures/F01_SREBF_module_mapping.pdf`、`F01_recovered_SREBF_activity.pdf`：技术原因与修正后活性。
- `figures/F02_mTORC1_singlecell_sources.pdf`、`F02_EC_filter_difference.pdf`、`F02_mTORC1_spatial_source.pdf`：三类来源及过滤。
- `figures/F03_H05_formal_lipid_effects.pdf`、`F04_sterol_FA_specimen_deltas.pdf`：正式过程与直接比较。
- `figures/F05_prespecified_routes_status.pdf`：预登记路线的 NE；次要状态见 F05 其余两图。
- `tables/model_results.tsv.gz`：1328条完整登记结果；`model_registry_complete.tsv`、`FDR_family_summary.tsv`、`decision_table.tsv` 配套。
- `METHODS.md`、`module_status.tsv`、`FIGURE_INDEX.md`、`qc/`：实际方法、模块状态与核验。

## 复现

从项目根目录运行：

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_hypothesis_followup_v2/scripts/run.py --resume
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_hypothesis_followup_v2/scripts/run.py --replot
```

`--resume` 核对输入、代码、冻结配置及输出哈希；兼容完成结果复用。缺失产物可按冻结管线恢复，发生不兼容变化时停止覆盖并要求新版本。`--recompute` 在同一输入/配置下重算受管理阶段；不会重跑 GRNBoost2、CellChat 或 NicheNet。`--replot` 仅依赖结果表和绘图配置，不读 Seurat、不重新统计。无原始数据的接收者可直接执行 `Rscript --vanilla scripts/11_replot.R`，设置 `MP10_RESULTS_DIR` 指向本结果目录，并安装 `logs/sessionInfo_scores.txt` 对应绘图依赖。

所有运行限制单进程计算线程1，R使用 `--vanilla`；顺序复现不超过总4线程。输入绝对路径、大小、时间、SHA256见 `tables/input_manifest.tsv`；原大对象仅引用，不复制。交付体积见 `tables/resource_usage.tsv`。未运行湿实验、HD、环境RNA校正、全量GRN重建、第二套完整seed456 motif网络或新排名库搜索。
