"""Build human-readable delivery from actual final tables; no model selection or fitting."""
from common import *
from datetime import datetime
import zipfile,re,sys,shutil
v=pd.read_csv(R/'tables/validation_checks.tsv',sep='\t');rv=pd.read_csv(R/'tables/R_loader_validation.tsv',sep='\t');assert v.pass_check.all() and rv['pass'].all();assert (R/'figures/R_loader_readback_niche50.png').exists()
def md(d):
 def fmt(x):
  if pd.isna(x):return 'NA'
  if isinstance(x,(float,np.floating)):return f'{x:.4g}'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x)for x in row)+' |'for row in d.itertuples(index=False,name=None))
coverage=pd.read_csv(R/'tables/analysis_coverage.tsv',sep='\t');u=pd.read_csv(R/'tables/unit_selection_audited.tsv',sep='\t');a=pd.read_csv(R/'tables/association_results.tsv',sep='\t');b=pd.read_csv(R/'tables/niche_MP10_association.tsv',sep='\t');comp=pd.read_csv(R/'tables/niche_composition.tsv',sep='\t');val=pd.read_csv(R/'tables/regional_MP10_validation.tsv',sep='\t');d=pd.read_parquet(R/'tables/Cervilla_cell_formal_metadata.parquet');reg=pd.read_parquet(R/'tables/Cervilla_epithelial_regions.parquet');ref=pd.read_csv(R/'tables/reference_membership.tsv.gz',sep='\t');nt=int(ref.partition.eq('train').sum());treg=d[d.fine_prediction.eq('Treg')];endo=pd.read_csv(R/'tables/Endothelial_mixture_diagnostic.tsv.gz',sep='\t');ctx=pd.read_csv(R/'tables/association_tissue_context.tsv',sep='\t');sens=pd.read_csv(R/'tables/association_Endo_mixture_sensitivity.tsv',sep='\t')
# Input identity: small metadata/stat checks, no full-data hashing.
inp=pd.read_csv(R/'tables/input_manifest.tsv',sep='\t');paths=[('Cervilla_bin_to_cell',A/'segmentation_assignment_v2/Cervilla/bin_to_cell.parquet'),('Cervilla_4um_ROI_source',A/'segmentation_restart_v1/Cervilla/seurat.rds'),('Cervilla_cell_counts',A/'segmentation_assignment_v2/Cervilla/cells.h5ad'),('Cervilla_Seurat',A/'segmentation_assignment_v2/Cervilla/seurat.rds'),('Cervilla_nuclei',A/'segmentation_assignment_v2/Cervilla/inputs/nuclei.npz'),('Cervilla_assignment',A/'segmentation_assignment_v2/Cervilla/inputs/assignment_mask.npz'),('Cervilla_transform',A/'segmentation_assignment_v2/Cervilla/metadata/transform.json'),('D1_transform',A/'segmentation_assignment_v2/D1_Johnson/metadata/transform.json'),('Cervilla_bin8_spatial_nodes',A/'hd_epi_niche_dual_methods_v1/Cervilla/nodes.parquet'),('D1_bin8_spatial_nodes',A/'hd_epi_niche_dual_methods_v1/D1_Johnson/nodes.parquet'),('frozen_Top20',P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'),('historical_MP10_code',P/'neighborhood_mp10_calibration_v5/code/common.py'),('MP10_groups_code',P/'mp10_immune_endothelial_mechanisms_v1/scripts/01_pseudobulk.R'),('type_dictionary',P/'mp10_sc_c2l_final_annotation_v1/tables/final_label_dictionary.tsv'),('Cervilla_existing_BANKSY',A/'hd_bin8_banksy_v1/Cervilla/metadata.tsv.gz')]
for section in ['Cervilla','D1_Johnson']:
 assets=json.loads((A/'segmentation_assignment_v2'/section/'metadata/assets.json').read_text());paths.append((section+'_HE',H/assets['image_project_relative']))
inp=pd.concat([inp,pd.DataFrame([dict(role=role,path=str(path),size_bytes=path.stat().st_size)for role,path in paths])],ignore_index=True).drop_duplicates('role');inp['mtime_ns']=[Path(x).stat().st_mtime_ns for x in inp.path];inp.to_csv(R/'input_manifest.tsv',sep='\t',index=False);put(inp,'input_manifest.tsv')
# Compact per-analysis-object review table; detailed raw panels remain in server tables.
cols=['cell_id','x_um','y_um','nCount_RNA','nFeature_RNA','fine_prediction','fine_pruned','type_supported','annotation_confidence','annotation_reason','epi_anchor','epi_mixed_flag','mp10_epithelial_eligible','mp10_single_score_eligible','MP10_detected','mp10_historical_auc']
for c in ['fine_delta','fine_max_score']:
 if c in d:cols.append(c)
compact=d[cols].copy();compact['section_id']='Cervilla';compact.to_parquet(R/'tables/analysis_cells_compact.parquet',index=False)
# Scope/status rows are categorical and separate from association results.
status=pd.DataFrame([['Cervilla','fixed_ROI_bin8_RCTD','completed','partially_evaluable','not_an_association'],['Cervilla','fixed_ROI_bin16_RCTD','completed','partially_evaluable','not_an_association'],['Cervilla','full_existing_cells_fine_reference','completed','partially_evaluable','not_an_association'],['Cervilla','Treg_identity_and_endpoints','evaluated_not_identifiable','not_evaluable','not_tested'],['Cervilla','Endothelial_and_MP10','completed','partially_evaluable','density_no_clear_support; distance_exploratory_positive'],['Cervilla','niche50_and200','completed','evaluable_composition; partial_MP10_coverage','K4_no_BHq_below_0.05; K2_200_sensitivity_support'],['D1_Johnson','light_target_registration_check','completed','partially_evaluable','not_tested']],columns=['section','branch','execution_status','evaluability','association_status']);put(status,'completion_status.tsv')
result=f'''# Cervilla MP10–Treg / 内皮空间验证结果

完成日期：{datetime.now().isoformat(timespec='seconds')} (Asia/Shanghai)。这是单张切片内的空间支持分析，不是患者重复验证。全部授权可执行分支已运行；Treg 相关推断因身份信息不足明确记为不可评估。未运行通讯、配体受体、NicheNet、TF 或恶性判定。

## 完成状态与正式单元

{md(status)}

本轮冻结的操作性主线为 **既有 cell 的保守谱系层用于位置与邻域 + 上皮限定的 64 μm 网格内连通小区域用于 MP10**。并非每个细胞都有可靠 MP10，也不意味着 cells 普遍优于 bins。全片直接细分参考预测覆盖93,593个原分析对象、60,670,091原始UMI；原158,144几何对象及QC/counts未改写。8个固定ROI的三路比较如下：

{md(u[['unit','n','fitted','reliable_fit','Treg_candidates','Treg_multigene','Endo_multigene','Epi_score_eligible','mixed_Epi_fraction','selected']])}

8 μm 保留原始定位优势，但相同ROI仅6个Epi高贡献bin通过个体MP10测量资格；16 μm虽提高检测量，按原生selected-model规则，Epi高贡献bin中约30%仍有Epi/胶原共检出，不能用模型Epi比例恢复纯上皮RNA。cells约30.5%同类共检出也不被当作纯上皮评分。**修正后的比较不支持cells比16 μm混合更少，bin/cell没有明确通用赢家。** 保留已冻结cells操作性主线，是为核支持位置与Epi-only连通区域求和，不能宣称身份或评分准确率更高；主结果因此保持探索性。8 μm原始信号应配合使用，不按Treg数量、分配UMI或显著性重选。16 μm从全部297,675 supplied8 bins聚合成75,293格，64,528,669 UMI逐基因及全局守恒；74,730格通过本轮宽松表达QC。没有先删除低UMI8 μm格。

## 目标可信程度

- 冻结30类精细参考训练实际{nt:,}个细胞；按 `orig.ident` 样本留出，不宣称患者独立。TIGIT不参与分类/RCTD判别。Treg内部参考精度在预计100/300/1000 UMI约0.61/0.80/0.89，且CD4/Tfh仍有混淆；这不是HD准确率。
- 上皮参考1800个细胞中，KRT5/KRT14/KRT6A/KRT6B/KRT17/COL17A1分别在352/165/239/64/1144/161个细胞检出，说明参考存在这些表达模式，但不是独立病理亚型验证。当前SingleR实际DE特征包含KRT7/KRT19，而上述基底/鳞状面板、EPCAM/UPK2未入选；它们仍进入原始RNA复核面板，不能把“基因存在”等同“分类器实际使用”。这限制对基底/鳞状上皮覆盖的强断言。Treg的FOXP3/IL2RA/CTLA4/IKZF2实际进入判别，LRRC32未入选，TIGIT明确排除。
- Cervilla有{len(treg)}个原始Treg预测候选，其中{int((treg.T_detected>=1).sum())}个具有至少1个原始T-panel基因，{int(((treg.raw_FOXP3>0)|(treg.raw_IL2RA>0)).sum())}个有FOXP3或IL2RA；没有足够的成组Treg/T-lineage联合证据。候选组平均/检出图也不支持把全部候选当作Treg。8/16 μm RCTD、既有4 μm保守归属ROI均未提供可信补救。**Treg、Treg/CD4、Treg/T、Treg距离及反向Treg/Tconv/Tfh比较全部NA**，不是零Treg，也不是不存在联系。
- 204个Endothelial支持候选有参考和至少2个原始内皮marker；PECAM1/VWF等群体证据明显。108个同时有至少2个Epi或Fib marker，因此不解释为完整、纯净血管细胞图谱。固定排除这类共检出后剩96个位置用于一次敏感性分析。内皮亚型未作为正式终点。
- 38,027个支持上皮对象作为niche中心，30,009个未标记混合的上皮对象可参与Epi-only区域构建；仅1,394个单对象达到预定MP10支持条件。8,008个候选连通区域中1,209个可评分。所有20基因均在完整18,085基因输入中，未使用HVG/BANKSY表达或TME表达评分。
- 主MP10严格复用历史**平均并列秩**公式，K=905；与原SciPy实现1800个检查实例最大误差0。原生AUCell随机并列/排除零秩仅为敏感性。合格区域半深度排序ρ={val.loc[val.comparison.eq('MP10_half'),'spearman'].iloc[0]:.3f}，预计3000UMI排序ρ={val.loc[val.comparison.eq('MP10_depth3000'),'spearman'].iloc[0]:.3f}；检测量调整后仍只能做区域层面解释。部分合格区域AUC=0表示已观察基因未进入高秩窗口，并不等于未测量；不合格区域不分Low。

## 50 / 200 μm 空间关系

{md(coverage)}

中心自身被排除，邻居池是全部原QC cells，包含周围Epi和全部可支持TME，Unknown从组成分子和条件分母同时排除。50/200为累计真实物理圆，不是环带、kNN或图跳数。类型贡献是核支持对象计数，不是真实无漏检细胞数。区域MP10没有复制成独立细胞观测。

主要条件结果（效应为标准化秩斜率；400 μm空间区块，499次）：

{md(a[a.model.eq('conditional')&a.family.eq('continuous_primary')][['radius_um','predictor','n','n_blocks','effect','low','high','P','BH_q']])}

**内皮支持密度未检出明确关联；距离有小幅正向关系，即较高MP10偏向离已识别内皮位置更远，而非内皮富集。** 最近距离本身不依赖邻域半径，两尺度共享916个可观察区域，不能当成两次独立复制。更大800 μm区块的距离检验P=0.048/0.040、敏感性家族q=0.192；排除混合Endo后距离q约0.099。加入已有BANKSY区域、边缘及可注释支持的条件分析方向仍相同，但这是同一片内的相关证据。因此总体结论为**探索性距离信号，未达到跨质量设定均稳定的强支持**。

{md(sens[['radius_um','predictor','n_target_sensitivity','n','effect','low','high','P','BH_q']])}

距离超过最近观察边界时删失；缺测、未拟合与零值分开。现有组织区域条件分析只将区域质心附近<=8 μm的原bin域作为空间背景协变量，不转移细胞身份。未据其结果改写niche或MP10定义。

## niche 组成与 MP10

两个尺度均按冻结规则选K=4，并保留K=2。50 μm K4的种子最低ARI0.986、区块ARI0.972；200 μm分别1.000、0.952。空间连续性不是准确率。组成描述：

{md(comp[['radius_um','niche','n','composition_description','mean_known_fraction']])}

这些niche主要解释上皮优势、上皮–间质梯度及较大尺度少量髓系/浆细胞成分；没有形成可验证的Treg niche。图像和原始Epi/胶原信号支持广义组织差异，但共享marker不是独立真值。两尺度编号不对应，交叉表在 `niche_scale_crosstab.tsv`。

K4的MP10主要家族没有q<0.05（最小q={b.loc[b.family.eq('niche_primary'),'BH_q'].min():.3f}）。50 μm niche1较高、niche4较低和200 μm上皮优势niche2较高仅为方向性趋势。50 μm间质优势niche3只有1个可评分区域，正式P撤为NA；不能把该niche判为MP10低。K2的200 μm粗结构敏感性出现上皮优势区域较高MP10(q=0.016)，但不据此替换冻结K4主结果。H50/H10/H05与L50的数量、比例及描述性优势比完整保留，不选择最显著分组。

**解释范围**：稳定组成差异成立；细粒度niche–MP10关联未检出明确多重校正支持。评分资格主要偏向含足够可信上皮的区域，不能推广到全部38,027中心或上皮稀少niche。

## D1 轻量检查

原364,805 QC8 μm bins及裁剪原点保留。原始T-lineage多基因支持14,443 bins，内皮多基因支持2,993，Treg-panel与T-lineage联合原始线索232，均不等于已确认亚型。既有 `conflict_focus` 和一致对照ROI的H&E/原始Epi/T/Endo图已实际查看；可支持粗组织配准，不能凭图定Treg。3000个固定抽样Epi候选全部20个MP10基因在输入内，单bin检出成员中位数1、126个检出>=5，提示仍有单bin评分稀疏问题。没有对D1做新全片精细分类、RCTD、niche或关联。D1评估结论为**部分可评估，值得后续目标验证但本轮不扩展**。

## 空间复核与后续适用性

已查看固定T-lineage、vascular、dense、far-UMI及补充原始panel图。`vascular_focus` 的命名仅来自预先选中的内皮panel热点，不是病理确认的血管ROI；H&E显示核支持结构，但现有扩展边界仍可能跨越邻近细胞。D1冲突ROI有组织化Epi表达，T/Endo信号较稀疏。无输入覆盖范围不作为目标阴性区域。

可进入下一轮的是：可追溯的上皮限定区域评分、广义上皮–间质组成、保守内皮支持位置及其不确定性。**当前不具备Treg特异通讯推断条件**；内皮相关结果也需独立细胞/血管身份或蛋白验证、更多切片及深度稳健性后再解释机制。未检出明确支持不是证明无生物学联系，本轮不输出通讯增强结论。

## 完成、验证与限制

Python {len(v)}项针对性核验、R原对象实际重读均通过：ID、整数counts、坐标、原QC、16 μm逐基因守恒、Epi-only区域求和、邻域分母/中心排除及两尺度独立字段。重读对象生成 `R_loader_readback_niche50.png`。原counts及正式上游文件只读，没有保存新全量Seurat/H5AD。

过程修复透明记录：RCTD工作进程命名空间问题以进程内FORK适配解决；原生输出的位置权重列按实际第一/第二类型解码；最终又核对get_beta_doublet，对singlet采用第一类=1，只有certain doublet保留两类权重，未定贡献为NA。该修正撤回“cells混合低于16 μm”的初始理由；修正前汇总保留，原冻结cell主线/MP10/niche和关联未变，单元选择置信度明确为provisional；历史MP10平均秩在关联前纠正；全片邻域训练索引修复为位置索引。最终统计检查发现极小niche组不能借总样本量产生P，增加每比较侧>=20区域且>=8区块要求，50 μm niche3撤去P。后续只增加任务要求的组织背景和一次内皮混合敏感性，未改动主单元、评分、中心、K或主要效应。

没有未启动却标为完成的计算分支。未进行有效检验的目标是Treg及其相对比例/距离/反向检验；D1正式关联依范围未做；所有通讯分析依停止条件未做。主要剩余限制为参考跨平台适用性、Treg稀疏、细胞归属混合、内皮不完全检测、Epi评分覆盖偏倚、单片空间自相关及边界删失。
'''
(R/'RESULTS.md').write_text(result);(R/'REPORT.md').write_text(result)
readme=f'''# Cervilla precommunication delivery

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
root = Path({str(R)!r})
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
'''
(R/'README.md').write_text(readme)
# Supplemental settings and statistical-validity repair, separate from pre-association freeze.
(R/'execution_addendum.json').write_text(json.dumps(dict(primary_unit_and_niche_unchanged=True,primary_MP10='historical_average_ties',post_fit_validity_repair='binary contrasts require >=20 regions and >=8 blocks per side; withdrawn P for 50um niche3 with n=1',bounded_sensitivity=['800um blocks for both main endpoints','known-fraction adjustment','existing BANKSY context and observed edge','exclude Endo with >=2 Epi or Fib markers'],Treg='not_evaluable; no fake zeros',time=datetime.now().isoformat()),indent=2))
# Count compact outputs; resources are observed, not invented whole-job memory peaks.
metrics=[]
for f in (R/'logs').glob('*.log'):
 txt=f.read_text(errors='replace');m=re.search(r'Maximum resident set size \(kbytes\): (\d+)',txt);t=re.search(r'User time \(seconds\): ([\d.]+)',txt)
 if m:metrics.append(dict(step=f.stem,recorded_process_max_RSS_KiB=int(m.group(1)),recorded_user_seconds=float(t.group(1))if t else np.nan))
put(pd.DataFrame(metrics),'recorded_process_resources.tsv')
# Preview/image review notes: exact observations, not pathology truth.
(R/'VISUAL_REVIEW.md').write_text('''# Visual review\n\nActually opened: raw T-focus/vascular/dense/far-UMI previews; registered mask/8μm/16μm/cell vascular ROI; full regional MP10 and r50 niche; D1 conflict H&E and raw Epi/T/Endo panels; target-group marker detection plot; final registered overview. Nuclei and inferred extensions are explicitly different boundaries. The vascular ROI label denotes a selected raw-panel hotspot, not confirmed vessel pathology. Raw T/Treg signals are sparse; supported EC group has PECAM1/VWF but mixed Epi/collagen evidence persists. D1 Epi pattern is organized and registration is plausible at this level; no subtype identity is inferred visually. Missing supplied areas remain unmeasured.\n''')
# Lightweight archive. No references, raw/full-expression matrices, model caches or giant per-bin raw-panel tables.
keyfig=['Cervilla_core_registered_overview','MP10_epithelial_regions','target_group_marker_support','ROI_bin_cell_T_lineage_focus','ROI_bin_cell_vascular_focus','ROI_targets_MP10_vascular_focus','ROI_model_and_score_far_UMI_hotspot','niche_r50','niche_r200','niche_composition_r50','niche_composition_r200','niche_stability','MP10_historical_depth_sensitivity','Endo_primary_and_sensitivity','niche_MP10_r50','niche_MP10_r200','D1_light_conflict_focus','D1_light_reference_agree_control','R_loader_readback_niche50']
allowed_tables=['analysis_cells_compact.parquet','Cervilla_cell_neighborhood_r50.parquet','Cervilla_cell_neighborhood_r200.parquet','region_result_r50.tsv.gz','region_result_r200.tsv.gz','epithelial_region_membership.tsv.gz','Endothelial_mixture_diagnostic.tsv.gz']
archive=R/'review.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6)as z:
 for f in list(R.glob('*.md'))+list(R.glob('*.yaml'))+list(R.glob('*.json'))+[R/'input_manifest.tsv']+list(R.glob('*.sh')):
  z.write(f,f.relative_to(R))
 for f in (R/'code').glob('*'):
  if f.is_file() and f.suffix in ['.py','.R','.sh']:z.write(f,f.relative_to(R))
 for f in (R/'tables').iterdir():
  if f.name in allowed_tables or (f.suffix=='.tsv'and f.stat().st_size<3_000_000):z.write(f,f.relative_to(R))
 for name in keyfig:
  f=R/'figures'/f'{name}.png'
  if f.exists():z.write(f,f.relative_to(R))
size=sum(f.stat().st_size for f in R.rglob('*')if f.is_file());resources=dict(completed=datetime.now().isoformat(),new_directory_bytes=size,review_zip_bytes=archive.stat().st_size,default_CPU_affinity=4,BLAS_OMP_threads=1,RCTD_workers=4,recorded_max_single_process_RSS_KiB=max([x['recorded_process_max_RSS_KiB']for x in metrics],default=0),observed_RCTD_parent_plus_worker_RSS_note='~16GB sum of RSS observed; FORK shared pages included; not a measured unique-memory peak',temporary_large_exports_removed=not any((R/'temporary').glob('*.mtx')),night_window='23:00–07:00 Asia/Shanghai',raw_inputs_modified=False)
(R/'RESOURCE_USAGE.json').write_text(json.dumps(resources,indent=2));(R/'RESULTS.md').write_text(result+f'\n新增目录占用(含review)：{size/2**30:.3f} GiB；review.zip：{archive.stat().st_size/2**20:.1f} MiB。实际资源记录见RESOURCE_USAGE.json和logs；默认4 CPU，未升级共享环境。\n');(R/'REPORT.md').write_text((R/'RESULTS.md').read_text());
# Refresh changed reports and resource record without duplicate zip members.
tmp=R/'review.tmp.zip'
with zipfile.ZipFile(archive)as old,zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED,compresslevel=6)as new:
 for info in old.infolist():
  if info.filename not in ['RESULTS.md','REPORT.md','RESOURCE_USAGE.json']:new.writestr(info,old.read(info.filename))
 for name in ['RESULTS.md','REPORT.md','RESOURCE_USAGE.json']:new.write(R/name,name)
tmp.replace(archive);print('DELIVERY_FINALIZED',json.dumps(resources))
