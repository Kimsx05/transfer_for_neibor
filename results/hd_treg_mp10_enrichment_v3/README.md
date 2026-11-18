# hd_treg_mp10_enrichment_v3

首先阅读REPORT.md。定义与统计见METHODS.md、frozen_config.yaml；所有计划行(含NA)保留。

```python
import sys
sys.path.insert(0, '/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_enrichment_v3/code')
from load_result import load_enrichment, attach_to_regions
result = load_enrichment('Cervilla')
T = result['T']             # 所有T32，含T_eligible/FOXP3_supported/strict
pairs_summary = result['neighbors']  # T_eligible × 50/200μm；U_r标记
# metadata = attach_to_regions(metadata, 'Cervilla')  # 唯一region_id，新增v3_前缀
```

```r
source('/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_enrichment_v3/code/load_result.R')
z <- load_enrichment('Cervilla')
```

R有arrow时直接读取；否则调用已冻结Python环境到临时TSV，读完删除。region_id不是bin/cell ID。不得将T32标签传播成成员bin的真实细胞身份。成员与原始配对仍引用`../hd_treg_mp10_spatial_all_v1/<section>/T32_members.parquet`、`Epi64_members.parquet`、`region_pairs.parquet`；原表达矩阵不复制。

`rerun.sh`逐步续跑，检查点校验源码/输入身份；`replot.sh`只从现有表重绘。一般不需要FORCE。手动修改定义应新建版本，不与本结果混用。未发布或上传。

`review.zip`为轻量审核包；完整图件PDF和紧凑区域表留在结果目录。旧原始证据图册通过ROI_interpretation.tsv引用，并选固定代表图纳入review。所有结果限于区域空间关联。
