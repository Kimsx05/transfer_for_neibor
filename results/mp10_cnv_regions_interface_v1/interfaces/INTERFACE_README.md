# 只读位置接口

主键：dataset_id + section_id + position_id（保留完整barcode后缀）。265,002行，每个有效观测位置一行；外侧不按多个region重复计数。

- region_id：最近有效界面所属Epi-rich region。position_region_id：位置实际所属区域，外侧为Outside；无界面时相应键缺失或保留位置区域。
- region_label / region_evidence / region_confidence：组成支持区域，不能解释为病理确诊肿瘤。
- boundary_id：对应boundary_summary.tsv.gz与outline GeoJSON中的单条有效线段。
- distance_to_boundary：非负欧氏最短距离；signed_distance：内部负、外侧正；signed_distance_available=False时距离NA。
- distance_unit / coordinate_unit：本片原图像像素；coordinate_reference为原图像绝对路径。不同图像像素不是同一物理单位。
- distance_dNN：signed_distance / 本片中位最近邻距，适合相对分辨率展示，不是微米。
- distance_band：inner_core、inner_edge、outer_neighborhood、far_tissue或not_evaluable。
- evaluation_status=evaluable：可进入主分带描述；censored_near_support_edge距离可读但主分析排除；no_valid_interface_in_component不能当作0。
- distance_to_support_edge：距有效观测支持边距离，包含QC缺失与组织/截断不明边。
- patient_id仅为来源可核实ID；record_id包含匿名/假定分析单位，不能统称独立患者。

边界表中eligible_interface=True才是区域内外界面；missing_QC_observation_edge为完整tissue网格中有位置但本次QC观察缺失；tissue_or_capture_edge_unresolved无法区分组织外缘和捕获取材截断。support_topology另分outer_envelope与interior_hole。参照methods.md使用，不把缺失当作阴性。
