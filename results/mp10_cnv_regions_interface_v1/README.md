# MP10 CNV regions/interface v1

完成63张Visium切片；主区域为Epi-rich region，不是病理肿瘤标签。先读 report.md、methods.md 和 review.zip。

只读接口：interfaces/position_interface_distances.tsv.gz。唯一键dataset_id、section_id、position_id；主评价筛evaluation_status=evaluable。signed_distance内部负外部正，单位为本片图像像素；distance_dNN可用于分辨率相对展示。region_id是参照边界所属区域，position_region_id为该位置本身区域。

缺失不填0，患者与匿名记录分开。所有旧输入只读。
