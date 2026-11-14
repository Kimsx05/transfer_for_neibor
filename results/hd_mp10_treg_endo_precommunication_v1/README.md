# Cervilla precommunication delivery

完成内容与结论见 [RESULTS.md](RESULTS.md)，方法见 [METHODS.md](METHODS.md)，图索引见 [FIGURE_INDEX.md](FIGURE_INDEX.md)。本目录为增量交付，上游输入不变。Treg不可评估不是零；正式MP10关联单位是Epi-only区域，不是复制区域分数后的单细胞。

## 读取

在项目父目录运行R (R --vanilla；原项目包路径见common.R)：

```r
source("highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/code/common.R")
source(file.path(root, "code/load_result.R"))
x <- load_precommunication(root, include_all=TRUE)
# 原158144对象、RNA counts和QC保留；新增pc_字段。原分析集合：
x_analysis <- subset(x, cells=colnames(x)[x$pc_in_analysis])
r50 <- load_precommunication_regions(root, radius=50)
r200 <- load_precommunication_regions(root, radius=200)
# 正式区域关联输入，而非将区域MP10复制到每个细胞：
r50 <- r50[r50$association_eligible, ]
```

Python使用原隔离环境 `../hd_epi_niche_multimethod_v2/.venv/bin/python`：

```python
import importlib.util
from pathlib import Path
root = Path('/home/data/t070721/codex_workspace/Bladder Metabolism/highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1')
s = importlib.util.spec_from_file_location("pc", root / "code/12_loaders.py")
pc = importlib.util.module_from_spec(s); s.loader.exec_module(pc)
metadata = pc.load_precommunication(root)
a = pc.load_precommunication(root, include_expression=True)  # upstream h5ad, backed read-only
# a.file.close() when finished
```

加载器按唯一ID挂回pc_字段；RDS/H5AD计数来自 `../segmentation_assignment_v2/Cervilla`。不生成另一份全量对象。区域表在 `tables/region_result_r50.tsv.gz` / `r200`，完整邻域条件分母在 `Cervilla_cell_neighborhood_r*.parquet`。D1只有轻量诊断表和图，不存在本轮正式结果对象。

## 主要文件

- `input_manifest.tsv`、`frozen_config.yaml`、`models/analysis_eligibility_frozen.json`、`niche_frozen.json`: 真实输入与冻结定义。
- `tables/unit_selection_audited.tsv`, `unit_comparison.tsv` (原冻结选择与修正前表另存): 同物理ROI8/16μm/cells选择与限制。
- `tables/reference_membership.tsv.gz`, `reference_model_features.tsv`, `annotation_validation.tsv`: 冻结参考抽样/样本留出/实际判别基因。
- `tables/Cervilla_cell_formal_metadata.parquet`: 93,593个原QC分析对象，raw panel/候选/互斥支持层/评分资格。
- `tables/epithelial_region_membership.tsv.gz`, `Cervilla_epithelial_regions.parquet`: Epi-only区域成员、分数、测量/分组资格。
- `tables/association_results.tsv`, `niche_MP10_association.tsv`, `association_tissue_context.tsv`, `association_Endo_mixture_sensitivity.tsv`: 完整主结果和有限敏感性，包括NA。
- `tables/*rctd*`, `*ROI_contributions.parquet`: 原生权重、状态、对应关系；uncertain/reject不等于可靠比例。
- `tables/section_target_evaluability_final.tsv`, `completion_status.tsv`, `analysis_coverage.tsv`: 可评估性与执行/关联状态分别记录。
- `models/Cervilla_bin16_counts.npz`: 唯一必要16μm稀疏counts，行序来自bin16_measurement.parquet，基因序来自hd_genes.tsv；不需要重保存2μm或8μm矩阵。
- `models/fine_reference.rds`, `fine_trained.rds`: 本轮精细参考训练必要缓存，不打入review。
- `figures/`: PDF/PNG与图例文件；`tables/`为直接绘图数据，H&E/masks通过manifest路径复用。
- `review.zip`: 报告、代码、关键图与紧凑结果，不含全基因矩阵、上游RDS/H5AD或参考缓存。

## 执行与重绘

```bash
bash highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/rerun.sh
bash highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1/replot.sh
```

rerun默认复用身份匹配的完成缓存，只补缺失步骤并更新报告；`--from-inputs`按保存定义完整重算本目录产物，禁止借此修改上游。重计算等待北京时间23:00–07:00，4 CPU affinity、BLAS/OpenMP1、RCTD4workers。最终绘图/表汇总可直接读取已有紧凑结果。可视化不用重新拟合。

必要依赖的绝对路径/尺寸/mtime在manifest；Python和R软件版本在tables/Python_versions.tsv、logs/sessionInfo_delivery.txt。16μm矩阵以及参考训练缓存均不可被review包代替；真实上游路径必须可访问。
