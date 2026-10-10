# Module 06 已完成：intra-spot 预测依赖

状态：**PASS_WITH_WARNING**。覆盖 63 张冻结切片、7 个数据集、265,002 个 tissue spots。仅运行 Module 06；未启动后续模块。

A：按批准映射汇总 9 个 major posterior mean abundances，逐切片训练；状态分布：{'PASS': 63}。
B：冻结 MP10 raw AUCell 为目标，技术 baseline = 原 Epi fraction + log1p(nCount)；完整模型增加全部 29 个非上皮 TME factors。同一 complete-case 集合与相同空间 folds。59/63 张切片完成五折空间评估；58 张 delta R² > 0；53 张同时达到 full R² > 0、delta R² > 0 的 importance 解释条件。

| dataset | n_sections | n_spatial_cv_complete | n_not_estimable | n_positive_increment | n_reportable | median_baseline_spatial_r2 | median_tme_spatial_r2 | median_delta_r2 | median_delta_rmse |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CanDis__UNVERIFIED | 13 | 13 | 0 | 13 | 11 | 0.01264 | 0.10521 | 0.06312 | 0.00128 |
| GSE171351 | 4 | 3 | 1 | 3 | 2 | -0.00927 | 0.04185 | 0.02014 | 0.00048 |
| GSE238145 | 10 | 8 | 2 | 8 | 8 | 0.00926 | 0.05856 | 0.0545 | 0.00077 |
| GSE285715 | 6 | 5 | 1 | 5 | 5 | -0.0283 | 0.11208 | 0.18311 | 0.00258 |
| GSE319536 | 16 | 16 | 0 | 16 | 15 | 0.03692 | 0.14492 | 0.07748 | 0.00142 |
| Inhousedata__UNVERIFIED | 7 | 7 | 0 | 6 | 5 | 0.04765 | 0.1715 | 0.08173 | 0.00107 |
| NPJ__UNVERIFIED | 7 | 7 | 0 | 7 | 7 | 0.11593 | 0.36832 | 0.21372 | 0.00283 |

以上是切片层面的描述性预测结果，不是患者独立重复的推断、显著性检验、通讯或因果证据。负提升保留，不计为技术失败。

## 主要交付

- `cellcell_importance.tsv`：逐切片 major 网络 raw 与 standardized importance。
- `mp10_baseline_vs_tme.tsv`：相同 spots/folds 的未截断空间预测表现与增量。
- `cellcell_dataset_median_standardized_importance.tsv`、`mp10_dataset_median_standardized_importance.tsv`：按 dataset 的中位数及有效切片数，同时保留全量与性能合格汇总。
- `abundance_mp_technical_associations.tsv`：与 importance 分开的方向关联表。
- `by_dataset/`：模型、完整预测、输入、性能、审计与 PDF + 600 dpi PNG。
- `review.html`、`review_bundle.zip`：紧凑审阅入口。
- `variable_dictionary.tsv`、`EVALUATION_METHODS.md`、`parameters.json`、`input_manifest.tsv`、`acceptance_validation.tsv`、`handoff.json`。

## 限制与来源

MISTyR 1.18.0（原文 1.10.0）与 distances 已在模块私有 R 库安装。官方集成输出会截断负 R²，本项目另存未截断 OOB/内部 meta-fold 指标；空间评估完全独立重训训练 fold 的森林与 OOB 校准器。A 仅内部验证；B 为切片内连续区域外推，无 buffer，不代表新患者验证。

标准化 importance 的负值表示低于相对平均贡献，不表示抑制。丰度与 counts 同变、变量相关性和混合 spots 均限制生物解释。patient_uid 和未核实 dataset 标签沿用冻结输入，不以 section 代替患者。

方法与解释详见 [EVALUATION_METHODS.md](EVALUATION_METHODS.md)。已在预测依赖/增量预测信息层面停止。

原生输出缺陷修正：MISTyR 1.18.0 的系数表为 5 列表头/6 列数据，原生 importance 误用了截距 P 值。已逐 target 用保存的 OOB 数据重建校准回归，验证真实 slope P 值并修正汇总权重；保留 native_standardized_importance、原生文件与完整缺陷审计。此修正不改变预测性能。
