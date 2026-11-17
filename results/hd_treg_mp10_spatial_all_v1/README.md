# hd_treg_mp10_spatial_all_v1

实际完成五片全量可执行分析。先读 [RESULTS.md](RESULTS.md) 和 [METHODS.md](METHODS.md)。主结果未检出跨片可重复正向关联；部分终点可评估性不足，详见逐行状态。

## 目录与依赖

结果根目录：`/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1`。
原8µm矩阵、空间节点、旧细胞结果与冻结参考均通过 `tables/input_manifest.tsv`、每片 `prepared.json` 和代码引用，没有复制整套counts/RDS/H5AD。
Python环境：`/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python`；R库路径见 `code/04_native_UCell.R` 引用的项目common.R，不修改共享库。

- 每片 `T32_regions.parquet`：所有32µm区域，分数、原始支持基因UMI/名单、资格和证据层。
- `T16_regions.parquet`：唯一测量尺度敏感性。
- `Epi64_regions.parquet`：上皮支持连通区域、MP10及合格分组。
- `*_members.parquet`：唯一原始bin→区域映射；不将区域分数复制成独立bin观察。
- `T32_neighborhoods.parquet`/`T16_neighborhoods.parquet`：50/200µm，两种来源规则、覆盖/资格/邻域均值。
- `region_pairs.parquet`：精确距离、共享来源标记；不包含全空间N×N图。
- `bin_evidence.parquet`：定位及原始目标panel表达，不是完整表达副本。
- `*_measure.parquet`：用于精确复跑的选定基因raw counts/完整背景秩缓存，非全基因矩阵。
- `tables/association_primary.tsv`、`association_all_results.tsv`：全部终点、效应、CI、P、BH及实际协变量。
- `tables/measurement_coverage.tsv`、`section_summary.tsv`：可评估性分母。
- `figures/`、`FIGURE_INDEX.md`：全图及索引；review.zip不带原矩阵和训练参考。

## 读取

```python
import sys
sys.path.insert(0, r"/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1/code")
from load_result import load_section
r = load_section("Cervilla", include_members=True, include_pairs=True)
t = r["T32"]
# 统计使用区域表；按bin ID挂回只附成员ID。
```

```r
root <- "/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_treg_mp10_spatial_all_v1"
source(file.path(root, "code/load_result.R"))
t <- load_treg_mp10(root, "Cervilla", "T32")
e <- load_treg_mp10(root, "Cervilla", "Epi64")
# existing <- attach_treg_mp10_membership(existing, root, "Cervilla")
```

## 复跑与重绘

`bash rerun.sh` 断点续跑，完成步骤跳过；输入身份与冻结配置保留检查。`bash redraw.sh` 仅从已保存表/原图像重绘，不重跑统计。完整强制复跑 `FORCE=1 bash rerun.sh`，按已冻结参数执行，不重新挑signature/ROI/最佳关联。上游依赖必须存在；不能用review包代替原输入。
本轮用户已授权立即计算，runner设置HD_ALLOW_DAYTIME=1；默认CPU总亲和性4、数值线程1。其他任务的夜间约束不由此改变。

## 字段与限制

`measurement_eligible`是本轮用途资格，原QC不改；`evidence`区分Treg-like、多基因一般T、信息不足。`MP10`/各UCell值对不合格区域也可保留计算诊断，但仅在资格True时用于正式终点；不可测量不等于Low。`tissue_status=unknown`仍为未知。
`effect`是标准化秩回归效应，`status`标明片内探索性或仅描述性；不是细胞比例、因果效应或多患者验证。
`support_area_um2`为成员网格面积，MPBC不能解释为病理组织面积。

检出基因名单和UMI直接来自raw区域计数。原始symbol/Ensembl完整映射在tables。未评价正式通讯，所有formal_communication_ready=False。
