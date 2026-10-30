# 全部普通 Visium 输入审计

母清单63片全部准入，v1 selected仅保存为历史标记。冻结RNA counts、AUCell、accepted posterior mean、身份索引均只读。两个试点counts/Q/coordinates经逐元素等价验证后软链接v1；meta/support为本轮独立文件。每片READY.json为队列入口；失败记录FAILED.json并继续其他切片。

Q以全部30参考类型归一化；Epi严格中位数状态；支持为原始mean abundance≥0.01的spot至少10且总abundance>1，不截断Q。固定100μm平台节距，200/105μm通讯范围；QC缺失邻居不改变节距。所有全片通过QC位置保留。没有足够H/L仍network_only，19/20评分覆盖单独汇总。

组织条件来自已存GEO SOFT明确注释并记录tissue_treatment_evidence.tsv；未核验组织/治疗保留unknown，不由文件名前缀猜测。身份沿用正式revised_by_user_v2交叉表，患者未核验不视独立患者。

恢复：/usr/bin/python3 scripts/prepare_inputs.py，跳过已READY切片。无大文件hash，无母RNA对象复制。

图像仿射残差仅作诊断。权威阵列索引来自 surviving VisiumV1@coordinates$row/$col，经严格section+barcode连接；63片coverage/topology原审计均PASS。保留每片相对残差及warning，不以任意图像像素阈值删片。geometry_denominator_diagnostics.tsv为固定半径图下潜在Q配对质量，不是推断通讯。spot_quality_diagnostics.tsv提供Epi状态内UMI/nFeature及Q权重质量诊断。

## 用户样本表补充标注（2026-10-27）

用户授权采用《空转+ScRNA样本汇总.xlsx》：NPJ 7片按直接sample别名标tumor、NMIBC、pT1；CanDis 13片按全队列均为MIBC/NMIBC的证据标tumor，逐片亚型不指定；Inhouse 7片仅补MIBC临床背景，组织条件unknown（T与RMYY对应及LN受累未知）。最终53 tumor、3 normal、7 unknown。证据层级、Excel行号和限制见workbook_annotations.tsv；患者核验、治疗、Q、评分及网络均未改变。
