# v9 最终结果摘要

完成：99/99，含 inhouse 57/57，不含 42/42。2026-10-19 北京时间 07:34 主分析完成；08:47 图包与敏感性完成；08:49 最终核验完成（902 项通过）。

主分析：REFIT、SOURCE_P05、k10、ALL_POLICY_UNITS、两层、各组合全部 niche 的 O/D；planned family 2,876 条，其中 2,621 条可估计，255 条科学性不可估计均保留。

| 结果类别（q<0.05） | 关系数 | 涉及组合数 |
|---|---:|---:|
| 组合内 BH 阳性（含全局阳性） | 172 | 59 |
| 全局 BH 阳性 | 85 | 51 |
| 仅组合内 BH 阳性 | 87 | 35 |
| 全局 BY 支持 | 17 | 16 |
| 全局 BH 阳性且 inhouse 有实际支持 | 37 | 22 |



完整的三个类别组合名单、真实 dataset 成员及各类别关系数见 [final_combination_categories.tsv](final_combination_categories.tsv)；各类别对应计数列大于 0 即为该类别名单。权威完整关系见 [primary_relations.parquet](primary_relations.parquet)，全局命中见 [global_primary_hits.tsv](global_primary_hits.tsv)，仅组合内探索性命中见 [within_combo_exploratory_hits.tsv](within_combo_exploratory_hits.tsv)。

85 条全局 BH 阳性中，spot 层 80 条、epi 层 5 条；O 重叠减少 69 条、增加 1 条，D 距离增加 14 条、减少 1 条。含 inhouse 的全局阳性共 39 条/23 组合，其中 37 条/22 组合有可估计且同向的 inhouse 支持。实际支持不表示 inhouse 单独显著。

## 重点组合

| 组合 | 数据集范围 | 重点关系（spot 层） | 单位平均 O 效应 | 全局 BH q | 全局 BY q | inhouse 可估计单位数 |
|---|---|---|---:|---:|---:|---:|
| 1111111 | 全部七个数据集 | niche3 的重叠减少 | -7.06 个百分点 | 0.00392448 | 0.0335211 | 2 |
| 1011111 | 除 GSE171351 外的六个数据集 | niche3 的重叠减少 | -6.82 个百分点 | 0.00164566 | 0.0140565 | 5 |
| 1010111 | 除 GSE171351、GSE285715 外的五个数据集 | niche3 的重叠减少 | -7.25 个百分点 | 0.0031158 | 0.0266137 | 3 |



niche 编号属于各自重拟合，不能仅凭同编号认定跨组合身份完全相同；上述 niche3 的环组成均以上皮为主、Matrix_CAF 为主要非上皮成分之一。这是 MP10 热点空间重叠关系，不是代谢物浓度或主动排斥机制的证据。

全七数据集的 spot 层 niche2、niche3、niche8 均通过全局 BH，q 分别为 0.02936、0.003924、0.01006；仅 niche3 通过全局 BY。epi 层 niche2 的组合内 q=0.02128，但全局 q=0.06029，属于仅组合内探索性结果。

## 敏感性与限制

选中的 3 个组合共 8 条关系：逐一剔除 dataset 或 analysis unit 后均未出现平均效应方向逆转；组成匹配后的 REFIT/FROZEN_LABEL 方向一致 8/8。k6/k18、几何支持、既有技术分层和组件等权的单位平均效应也均同向，完整数值见 [review/sensitivity_summary.tsv](review/sensitivity_summary.tsv)。这些是方向和影响诊断，不等于剔除后依然显著或排除单个样本对 P 值的影响。

1011111 的 niche1 虽是全局 BH 阳性，但 inhouse 的平均效应与总体反向，未被计入实际支持。inhouse 七个样本患者身份尚未完整核实，沿用独立样本政策；不能称为七名已核实独立患者。99 组合共享数据，不能视为独立重复验证；BH/BY 不消除结果选择偏差。

全部组合统计、辅助 family、紧凑标签、切片/单位效应及复现记录均已保留。3 个重点组合共 164 张 PNG 及对应 PDF，另有全局 overview。无尚未完成的计算项；身份核实与独立验证仍是科学解释的限制。

资源：最终核验时实际分配约 1.13 GiB（预算 3 GiB），累计清理本轮可再生临时文件约 6.94 GiB。新增摘要为轻量衍生文件，最终字节数见 audit/delivery_summary_manifest.json。

完整执行报告与核验见 [RUN_SUMMARY.md](RUN_SUMMARY.md) 和 [audit/final_validation.json](audit/final_validation.json)。
