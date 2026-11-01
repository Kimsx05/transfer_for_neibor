# v3.4 正式候选区域的使用与复现

[全队列标签表](final_spot_labels.tsv.gz) 中每行对应一个原始spot。用 `section_id + barcode` 合并，`lesion_location` 表示计算候选位置，`position_support` 和 `reference_state` 单列支持来源及冲突。主位置使用 inside、outside、unlocalized、quality_excluded。坐标 `x/y` 对应原位图像像素，`x_um/y_um` 采用既有普通Visium标称100 μm点距。原样本和分析单位沿用项目既有映射。

```python
import pandas as pd
from pathlib import Path

root = Path("results/mp10_pathology_regions_v3_4_final")
key = ["section_id", "barcode"]
labels = pd.read_csv(root / "final_spot_labels.tsv.gz", sep="\t",
                     dtype={k: "string" for k in key}).set_index(key, verify_integrity=True)
# adata.obs 的这两列应保存原 section 和原 barcode，包括原后缀。
order = pd.MultiIndex.from_frame(adata.obs[key].astype("string"))
assert order.is_unique and order.isin(labels.index).all()
selected = labels.reindex(order)
for col in ["lesion_location", "lesion_roi_id", "region_role", "position_support",
            "signed_distance_um", "exact_distance_status", "distance_band_200um"]:
    adata.obs[col] = selected[col].to_numpy()
```

## 局部修订

冻结非GSE319536切片及BC05、BC17的位置。其余切片复用 v3.3 的同片H&E配准、可靠参考核心和原作者S1图。用原图图例的RGB颜色（与最近色块欧氏差<28）提取 Tumor Cell＋Tumor Stroma 及 Normal Stroma＋Urothelium，在缓冲前的图上分别取得四连通父区。只启用含匹配可靠核心的父区。保留中心原作者类别，Others、TLS、免疫、肌肉、坏死等继续作为组织属性。[原作者资料](https://github.com/VividKai/SpatialBC) · [补充图](https://ndownloader.figshare.com/files/69527807)

对每个父区，在已有配准误差半径对应的像素邻域内计算父区所占比例。一个父区的支持仅在其比例≥0.15、能沿观测spot六邻接走到同父区同侧核心的部分生效。相反核心构成阻挡。多父区重叠时，同侧取最强的父区支持。待优化范围包括这些部分及一圈真实观测邻居，其外沿用 v3.3。禁止跨缺测、裂隙或独立组织块传播。

在上述范围做一次二元最小割，所有切片使用同一组参数。内外侧采用对称的可靠核心约束。inside相对outside的单点偏好由三项相加：v3.3位置先验（内+1、外−1），六倍的内外父区支持差，以及弱分子项 `0.2 × clip(已有域局部分数−0.42, −1, 1)`。同一BANKSY域连通部分的相邻点异色代价1.35，跨域相邻点0.60。核心和冻结位置的约束代价100，容量乘1000取整。两侧各自接受同色先验或核心支撑；完全翻转、且脱离同侧原位置和核心的孤立结果回到先验。局部恶性分数仅参与过渡部分判断。参数没有按切片或MP10结果调整。[配置](config.json) · [分割代码](scripts/regions.py)

对被红区包围的完整蓝区另做对象级判断。须满足全部外邻位实测、外邻域全inside、周围包含由肿瘤参考支持的原outside，且本区过半观测中心属于同一个可从可靠核心到达的肿瘤父区，该父区也出现在包围环上。存在正常核心、正常父区支持或正常类别中心时不作这项归并。分别检查原图的封闭区及最小割后的局部封闭区。满足条件时将整个观测区室赋为 inside。该规则不插补空位，不处理无作者证据的普通孔洞。作者图参与计算校准，所有输出均保留候选解释。

## 几何及分带

实际位置发生变化的切片重新建立ROI连通组、hex面、边段、技术缺位判定、已知区域K及距离。其余切片复用v3.3几何。最终表和几何包采用一致的 `:F:` 标识，ROI ID只在所属切片内使用。真实采样单点可以独立存在。

可观测候选外缘要求观测到内外点直接接触，外侧观测分量能达到组织外侧。QC缺测通道不提供这一路径。内部区室接触、切缘、捕获截断、缺损和未知接触分别记录。对于仍保留outside位置、但没有观测外侧通路的区室，清除对其他区室外缘的距离引用，`region_role` 记录区室状态，分带为 `outside_no_external_route`。没有外缘的inside块保持位置，距离NA。

主距离以共享边为零点，内部负、外部正。inside只参考自己的ROI，outside只参考同一组织块中的可观测候选外缘。已知K包含所有当前可定位的观测hex和满足 v3.3 全组织网格、H&E及一致周边条件的技术缺位。技术缺位只用于几何，不产生分子观测。令 g 为到可观测候选边的距离，q 为到K外未知区域的距离。绝对距离下界为 min(g,q)，直线路径可见时上界为 g；沿用25 μm采样和58 μm同块支持半径。q>g+1e−6且路径可见时才保存精确有符号距离。

100、200、300 μm 三种分带共享同一套轮廓。上界≤w判相应侧margin，下界>w判 `candidate_inside_beyond_w` 或 `candidate_outside_beyond_w`，其余情况保留未定或无外缘状态。精确距离NA与固定宽度分带未定分别存储。`observed_outer_contact` 与 `internal_compartment_contact` 表示点是否直接接触相应边段。`signed_author_parent_id` 的正负区分作者图中的肿瘤/正常父区，0表示该点没有直接父区ID。

每片的 `geometry.json.gz` 包含微米坐标ROI GeoJSON、主边段、未知区域边段、几何缺位及独立的小巢接触面。每个表用 columns/data 存储，主边段带原始barcode。可以直接读取：

```python
import gzip, json
with gzip.open(root / "by_dataset/GSE319536/GSE319536_BC13/geometry.json.gz", "rt") as f:
    geometry = json.load(f)
tab = geometry["roi_edges"]
edges = pd.DataFrame(tab["data"], columns=tab["columns"])
```

## 运行与上传

在项目根目录运行 `python3.12 -B results/mp10_pathology_regions_v3_4_final/scripts/finalize.py` 和 `python3.12 -B results/mp10_pathology_regions_v3_4_final/scripts/package.py`。前者按切片两片并发复用已有缓存，后者生成主表、核对重点对象并生成浏览入口。输入依赖本项目完整 v3.3 本地数据及其已保存的参考图和几何帮助函数。最终表、几何和浏览页面可独立读取。

有目标仓库写权限时，将完整目录复制到该仓库的 results 下，使用 `git add -f -- results/mp10_pathology_regions_v3_4_final`，并核对 `git diff --cached --stat` 和 `git ls-files -- results/mp10_pathology_regions_v3_4_final`。强制纳入交付目录可避免全局忽略PNG或压缩表。当前是完整本地交付，尚未发布到GitHub。
