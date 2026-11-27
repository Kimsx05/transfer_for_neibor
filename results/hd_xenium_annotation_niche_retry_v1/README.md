# BLCA 高分辨率细胞注释与共同 niche
主入口为 [REPORT.md](REPORT.md)、[输入清单](tables/input_manifest.tsv) 与 [逐片可读对象](tables/analysis_objects.tsv)。参数、依赖和作者代码定位见 [配置](shared/analysis_config.json)、[版本](shared/versions.json) 和 [方法](METHODS.md)。

每片 analysis_index.parquet 的索引是 original_cell_id。多切片拼表必须同时使用 sample。空间坐标 x_um / y_um 可以直接用于同片作图。major_identity、fine_identity、review_level、display_label、conflict、eligible_expression、niche_id、graph_degree 分别给出大类、可信细类、四层状态、展示标签、争议原因、niche 资格、平台共用编号、有效邻居数。unassigned_QC 保留观察行。graph_degree=0 表示缺少邻域信息。human_review=pending。

HD 同目录 cells.h5ad 的 X 是本轮唯一完整细胞×基因原始聚合矩阵。var.feature_id 与 var.symbol 保留重名探针，唯一化的 var_names 用于索引；不应把同名探针无意丢弃。成员表 final_label=0 为未分配，正值对应 cells.h5ad 的 obs.label；old_id 只标记 primary 对应关系。geometry_flag 的最新审阅结果读 analysis_index，GEX 图像支持另见 secondary_image_support.parquet，模型成分见 tacco_weights.parquet。niche 或重新评分应首先明确是使用所有 source，还是 primary 及形态支持子集。

Xenium 直接读取 input_manifest 中的原始 H5，只选 Gene Expression；用原始 barcode 与同片 analysis_index 对齐。old_gate_audit.parquet 保存旧拒绝原因与分差，prediction.parquet 保存新分类器候选、竞争项与分数。历史代谢读数另在 v1/v2 中按原单位保留；本轮未重算。

每片 latent.npy 按原始表达行序；低表达/基线几何不合格者为 NaN。HD 后续 GEX 图像质控可能进一步缩小 niche 集合。niche_graph.npz 的行序读 niche_graph_ids.tsv.gz。平台级 scvi_model/、latent_genes.tsv、cellcharter_model.pkl 及 gmm_training_indices.npy 支持在同一兼容环境复查模型。两平台同号 niche 各自独立。

