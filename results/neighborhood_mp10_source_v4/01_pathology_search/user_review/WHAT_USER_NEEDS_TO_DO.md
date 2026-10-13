请打开 review_index.html，按表顺序复核尚无可信病理标签的 ROI。最先五项：CanDis Sample8、GSE319536 BC11、CanDis Sample1、CanDis Sample3、GSE285715 T1。

在 ROI_review_to_fill.tsv 填 reviewed_label、review_comment、reviewer、review_date。允许值 tumor_epithelium / non_tumor_epithelium / other / mixed_needs_split / uncertain。混合区域需 parent_roi_id、polygon_coordinates 与 coordinate_transform 或 spot_members_file。

GSE319536：请向作者索取 Fig1.R 引用的每切片 SeuratObj/*.seurat.umap.spatial.rds 或 V2503.figure1.metainfo.summary.txt，以及 Official_Annot/Level1_Annot/Level2_Annot 字典和坐标配准。CanDis 与 NPJ：需要本地 Sample 名到公开标本的原始交叉表。GSE238145 的 tissue 只是 B12/B13 标本拆分，不能当病理。GSE285715 的 N 是邻癌来源而非全切片病理正常。

未知病理不影响本轮评分。所有图仅为证据，未替代用户判定。
