# 点表和观测域

`sections.tsv` 以 section_id 唯一关联 dataset、来源、核实患者 ID、analysis_unit、证据、原 H&E 及每片 spots_path 和 boundary_path。按这些路径读取 v2 的唯一点表，再以 section_id+barcode 连接。一个 MP10 字段供两种边界共同使用。用于发现的独立元数据在 immune_metadata.tsv。机器可读的点表字段类型见 [point_schema](point_schema.tsv)。

|字段|含义|
|---|---|
|x/y，x_um/y_um，array_row/col|原始图像坐标、物理坐标 μm、六角阵列位置|
|block，all_E_component，E_without_I_enriched_component|组织块，以及各自掩膜的内部连通分量；外侧原字段 −1|
|region_main，region_minus20|v1 主标签、固定去 Top20 版本；此轮不重聚类|
|all_E，E_without_I_enriched|主标签四类 E 或仅 ENN/ENS|
|expression_QC，UMI，nGene|v2 center_qc 且 UMI>0、nGene≥150，以及测序深度与检出基因数|
|MP10，new_epi_domain，P95_new|原冻结分数，主 E 评分域，以及此域每片共享 P95（含并列）；真实零保留|
|legacy_query，legacy_inside，legacy_tumor_epi，P95_old|历史域与历史空间分位标签，与新域分开|
|E/T/I_ab、E/T/I_f、total_abundance|c2l posterior mean 指定谱系估计量或占全部 30 类和的比例。T 不含 NK|
|dist_all_E，dist_E_without_I_enriched|同块同侧连通分量到实际共享 hex 边的距离，μm，内负外正|
|strict_*，status_*|截断或缺测边更近时 strict 为 NA；status 保留原因|
|accepted_label|经过既有验证对齐的 SAMPLE_R12_LEIDEN 0–8，在原 anchor 分母上解释|
|analysis_unit_uid，patient_uid，unit|保留原映射证据，统计键优先核实的患者，其次既有分析单位。两者均缺失者只进入分型描述；发现 bootstrap 用逐片暂时代码并另作排除敏感性|

局部表的 grid、window 对应组织块内网格；E_n、E_area、E_total、ab_T/I、lag_n、lag_T/I、outer_n、outer_ab_T/I、outer_total 均可追溯分子分母。边段原表保存两侧条码、有限端点、组件和外侧类型。逐箱的 block / 当前掩膜内部组件贡献见 [几何支持](../tables/bin_geometry_support.tsv.gz)。

逐箱的两侧组件贡献另外按同一组织块、同侧六邻接重新编号。side_component 是每片每边界的局部图标识，不能跨片合并。每箱组件分母经点表重计核验。[含外侧分量的支持表](../tables/bin_side_component_support.tsv.gz)
