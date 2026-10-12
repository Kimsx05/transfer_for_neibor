# 你现在需要做什么

本轮 niche 的计算不等待 ROI 判断。请先复核下列 5 个高优先级区域，再按需要处理其余 128 项。133 项是计算连通块，不等同于独立病灶；原 8257 个连通块无需逐一人工复核。

打开 [待填表](ROI_review_to_fill.tsv)，保留 roi_id，填写 reviewed_label、review_comment、reviewer、review_date。允许标签 tumor_epithelium、non_tumor_epithelium、other、mixed_needs_split、uncertain。不要填写不存在的病理依据。

若区域混合，不要求给整个 ROI 一个标签：在 review_comment 写明需拆分，另返回 polygon 坐标（所附图像像素坐标）或 spot_uid 成员表，保留 parent_roi_id；后续据此生成稳定 child UID 和新 ROI version。

## CanDis_2026_Sample8__uncertain_or_conflicting__0001
共同 anchors：10161；mean Epi=0.415；aneuploid=0，diploid=10459。

[证据图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample8__uncertain_or_conflicting__0001.png) · [原嵌入 H&E](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample8/HE_embedded.png)

明确问题：请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。
配准状态：IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION；请确认 spot 与组织边缘吻合。

## GSE319536_BC11__uncertain_or_conflicting__0001
共同 anchors：4669；mean Epi=0.536；aneuploid=0，diploid=4614。

[证据图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/GSE319536_BC11__uncertain_or_conflicting__0001.png) · [原嵌入 H&E](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/GSE319536_BC11/HE_embedded.png)

明确问题：请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。
配准状态：IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION；请确认 spot 与组织边缘吻合。

## CanDis_2026_Sample1__uncertain_or_conflicting__0001
共同 anchors：4291；mean Epi=0.729；aneuploid=0，diploid=4482。

[证据图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample1__uncertain_or_conflicting__0001.png) · [原嵌入 H&E](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample1/HE_embedded.png)

明确问题：请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。
配准状态：IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION；请确认 spot 与组织边缘吻合。

## CanDis_2026_Sample3__uncertain_or_conflicting__0001
共同 anchors：4269；mean Epi=0.726；aneuploid=0，diploid=4535。

[证据图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample3__uncertain_or_conflicting__0001.png) · [原嵌入 H&E](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/CanDis_2026_Sample3/HE_embedded.png)

明确问题：请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。
配准状态：IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION；请确认 spot 与组织边缘吻合。

## GSE285715_T1__uncertain_or_conflicting__0001
共同 anchors：3913；mean Epi=0.532；aneuploid=0，diploid=4165。

[证据图](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/GSE285715_T1__uncertain_or_conflicting__0001.png) · [原嵌入 H&E](/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_harmony_v3/05_user_review/evidence/GSE285715_T1/HE_embedded.png)

明确问题：请结合 H&E 判断：肿瘤上皮、非肿瘤上皮、其他组织，还是仍不确定？如 ROI 混合，请拆分并返回 spot 成员或 polygon；不要以 diploid 判正常。
配准状态：IN_BOUNDS_REQUIRES_VISUAL_CONFIRMATION；请确认 spot 与组织边缘吻合。

## 缺少什么

[逐切片图像清单](image_availability.tsv) 列明可用的嵌入 H&E、坐标和比例因子。已有图像属于 Seurat 保存的分辨率，不能代替原始全分辨率 H&E；尚缺病理标注。若无法判断，请针对对应 section_id 提供原始 H&E、配套 scalefactors_json.json / tissue_positions 文件及病理标注来源，不需要重做 CopyKAT。

## 如何返回

将填写后的 ROI_review_to_fill.tsv 保存到此目录或上传，必要时附 ROI 拆分 polygon/spot_uid 表。空白表示尚未复核，uncertain 表示已看证据仍不能判定。表内不预填人工姓名或日期。返回后只进入后续区域/界面确认；不会自动触发 MP10、解盲或来源验证。
