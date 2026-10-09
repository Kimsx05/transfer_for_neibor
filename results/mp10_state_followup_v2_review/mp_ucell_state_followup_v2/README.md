# MP10状态后续分析 v2

已完成六项固定比较、原usage衔接、sample配对edgeR QL、marker/来源/QC复核及冻结坐标出图。

工作目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_followup_v2`

首先打开 [REPORT.html](REPORT.html)，再看 [candidate_decision.tsv](candidate_decision.tsv) 和 [METHODS.md](METHODS.md)。

|比较|合格配对sample|dataset|检验基因|
|---|---:|---:|---:|
|r04_C8_vs_rest|28|4|15902|
|r06_C12_vs_C16|10|3|7298|
|r06_C12_vs_rest|25|4|15967|
|r06_C16_vs_rest|12|3|15739|
|r04_C8_vs_C6|15|4|11479|
|r04_C8_vs_C3|21|4|9676|

C8原MP10 usage在20/28合格样本更高；C12/C16有表达差异，但C16在全部10对中nCount/nFeature更高，且69.5%来自HRA两个sample。优先保留完整C8作独立marker评估，C12/C16作并行状态对照，不能宣布已验证的新亚群。

展示候选：C8 EGR1/GADD45B/ELOVL5；C12 HIPK3/GSE1/ELOVL5；C16 NRARP/IDI1/LDLR。各自局限和样本反例见03_marker_review。MP10非共享IDI1/INSIG1/LDLR等支持脂质/甾醇相关解释，但不是整条通路一致激活。

原NMF/Seurat/UCell/labels/UMAP全部只读。FULL含Normal，未删细胞。患者映射未核实，不能称独立患者验证。GSE326225保留完整视图但六项比较无合格配对；无配对图的模块不伪造结果。未运行空间、CNV、预后或额外聚类。

输出目录：00_audit输入/资格，01_usage_bridge原usage，02_pseudobulk_markers完整DE/过滤/pseudobulk缓存，03_marker_review短名单/功能/来源/QC，04_figures_by_dataset透明PDF+300dpi PNG及实际绘图数据，code与logs。

复现入口：`bash code/run_all.sh`。代码固定绝对输入路径，配置一致的聚合/DE检查点可复用，原UCell和UMAP始终读取。审阅包为`mp10_state_followup_v2_review.zip`，排除原Seurat、完整表达/usage矩阵及大型pseudobulk缓存，包含完整DE压缩表、报告、代码、关键图及图索引。
