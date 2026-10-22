"""Create concise audit summaries and the review report after independent validation."""
from pathlib import Path
import json,hashlib,datetime,pandas as pd,numpy as np,anndata as ad,shutil
r=Path(__file__).resolve().parents[1];assert not (r/'metadata/FROZEN.json').exists();project=r.parents[2]
a=ad.read_h5ad(r/'cell_annotation.h5ad',backed='r');o=a.obs;q=json.loads((r/'metadata/reconstruction_qc_summary.json').read_text());s=json.loads((r/'metadata/annotation_summary.json').read_text())
def table(df):
 def fmt(v):return str(v).replace('|','/')
 return '| '+' | '.join(map(fmt,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(fmt,row))+' |' for row in df.itertuples(index=False,name=None))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
# Disposition by final lineage, preserving independent low-depth and morphology flags.
rows=[]
for field in ['major_label','subtype_label']:
 for label,g in o.groupby(field,observed=True):rows.append(dict(grouping=field,label=label,n=len(g),low_information=int(g.low_information.sum()),qc_flagged=int(g.qc_flags.ne('none').sum()),morphology_uncertain=int(g.reconstruction_flag.ne('none_confirmed').sum()),inference_only=int(g.model_role.eq('inference_only').sum()),mixed_or_ambiguous=int(g.annotation_status.eq('mixed_or_ambiguous').sum())))
pd.DataFrame(rows).to_csv(r/'metadata/flags_by_final_identity.tsv',sep='\t',index=False)
st=pd.read_csv(r/'metadata/subtype_support_summary.tsv',sep='\t');st['supported_fraction_of_retained']=st.subtype_supported/len(o);st['supported_low_information_fraction']=np.divide(st.low_information_supported,st.subtype_supported,out=np.full(len(st),np.nan),where=st.subtype_supported>0);st['candidate_unresolved_fraction']=np.divide(st.candidate_unresolved_n,st.hierarchical_candidates,out=np.full(len(st),np.nan),where=st.hierarchical_candidates>0);st.to_csv(r/'metadata/subtype_support_summary.tsv',sep='\t',index=False)
# Hash immutable external reference decisions, coordinate inputs and task attachments.
ref=json.loads((r/'metadata/reference_manifest.json').read_text());sources=[Path(ref['object_path']),project/'subtype/v2/final_annotation_pre_reference/frozen_annotation_barcode_mapping_v2.csv']
for name in ['B_round2_final_annotation_frozen_v2.md','TNK_round2_final_annotation_and_parameters_frozen.md','Myeloid_round2_final_annotation_frozen.md','Fibroblast_round2_final_annotation_frozen.md','Endothelial_round2_final_annotation_frozen.md','B_cell_annotation_final_framework.md','Myeloid_annotation_final_frozen_framework.md','Endothelial_annotation_final_frozen_framework.md','fibro_CAF_annotation_final_plan.md']:
 sources.append(project/'subtype/v2'/name)
sources.extend([project/'subtype/v2/reference_training_v2/10_final_report/reference_v2_final_30type_freeze_report.md',Path('/home/data/t070721/.codex/attachments/6fdf867c-7225-4ebd-9266-2c75fd24c7b2/VisiumHD_permissive_QC_subtype_freeze_prompt.txt')])
for row in json.loads((r/'metadata/input_manifest.json').read_text()):
 p=Path(row['matrix']).parent;sources.extend([Path(row['image']),p/'spatial/tissue_positions.parquet',p/'spatial/scalefactors_json.json'])
provenance=[]
for p in sources:
 assert p.exists(),str(p)
 provenance.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
(r/'metadata/source_provenance_checksums.json').write_text(json.dumps(provenance,indent=2))
# Only regenerable files produced by this run are removed, never original data/models.
cleanup=[]
for p in [r/'.pip_cache',r/'inputs/reference_holdout_temporary.h5ad']:
 if p.exists():
  n=sum(f.stat().st_size for f in p.rglob('*') if f.is_file()) if p.is_dir() else p.stat().st_size
  cleanup.append(dict(path=str(p),bytes=n,reason='Task-only download cache or reproducible small held-out reference matrix; verified models, selections, seeds and predictions retained'))
  if p.is_dir():shutil.rmtree(p)
  else:p.unlink()
(r/'metadata/cleanup_manifest.json').write_text(json.dumps(cleanup,indent=2))
# Module guards prevent accidental rewriting once the frozen marker exists.
for p in (r/'code').glob('*.py'):
 if p.name=='marker_rules.py':continue
 text=p.read_text()
 if '# FROZEN_VERSION_GUARD' not in text:
  guard="# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.\nfrom pathlib import Path as _FreezePath\nif (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():\n    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')\n"
  p.write_text(guard+text)
cfg=json.loads((r/'metadata/run_config.json').read_text());cfg.update(workflow_status='computational_freeze_ready_pending_user_review',actual_sections=dict(full_section=0,partial_ROI=1,blocked=2),actual_scvi_sections=1,core_definition_note='Core chosen before expression-marker conflict review; MAD-flagged cells and low-information excluded from structure definition, later uncertain conflicts not asserted confirmed doublets');(r/'metadata/run_config.json').write_text(json.dumps(cfg,indent=2))
files=[p for p in r.rglob('*') if p.is_file() and not p.is_symlink()];env=sum(p.stat().st_size for p in files if '.venv' in p.relative_to(r).parts);out=sum(p.stat().st_size for p in files if '.venv' not in p.relative_to(r).parts);disk=dict(environment_bytes=env,deliverable_and_audit_bytes=out,total_bytes=env+out,removed_regenerable_bytes=sum(x['bytes'] for x in cleanup),measurement='apparent file bytes before report/checksum marker; filesystem allocation may differ');(r/'metadata/output_disk_usage.json').write_text(json.dumps(disk,indent=2))
statuses=pd.DataFrame([dict(状态=k,细胞数=v,占保留细胞=f'{v/len(o):.2%}') for k,v in s['annotation_status'].items()]);support=pd.DataFrame([dict(亚群=k,支持细胞数=v,解释='少于 5 个，逐个复核，不据此宣称群体稳定存在' if v<5 else '多细胞表达支持，仍待审阅') for k,v in s['subtype_supported'].items()]);validation=pd.read_csv(r/'metadata/reference_depth_validation.tsv',sep='\t');validation=validation[(validation.label=='macro avg')&validation.model.isin(['major','subtype'])][['model','condition','f1-score']];validation['f1-score']=validation['f1-score'].round(4)
report=f'''# Visium HD — frozen_subtypes_v1 计算注释冻结报告

生成日期：2026-10-21（Asia/Shanghai）。状态：**计算注释完成，待用户审阅；未经过病理专家审定。**

## 1. 完成范围

三个候选中：**0 张 full_section；1 张 partial_ROI（Cervilla）；2 张 blocked（D1_Johnson、UC_SCC）**。本轮到标签冻结为止，未运行 MP10/胆固醇评分、物理邻域/niche、CellCharter、spottedpy、空间阳性检验或细胞通讯。

| 候选 | bin 输入 | cell 重建状态 | 重建 / 保留 / 低信息 | 完成状态与缺口 |
| --- | --- | --- | --- | --- |
| Cervilla | 2 µm，4,722,911 bins，18,085 genes | 在现有已配准图像覆盖范围内可用，但原生图像细节受限 | 106,934 / 97,847 / 12,798 | partial_ROI；原始 Bladder_HD_he.tif 未恢复，1,274 个组织内 bins 超出图像 |
| D1_Johnson | 2 µm，6,362,240 bins | 不合格：预览图原生 5.59 µm/pixel 不足以核分割 | NA / NA / NA | blocked；需原始 JohnsonA1_D1_Scan1.qptiff 及可靠坐标变换，或可验证原始 cell counts＋分割映射 |
| UC_SCC | 仅 8 µm，695,646 bins | 不合格：无 2 µm 或合格作者细胞输入 | NA / NA / NA | blocked；需 2 µm＋坐标＋合格图像，或可验证原始 cell counts＋分割映射 |

partial_ROI 指覆盖**整幅现有已配准图像画布**，不是仅对三个展示 ROI 计算，也不代表整张切片完全恢复。D1 的 2/8/16 µm 是同一张切片，未重复计数。三个候选 donor 均 unknown，不推定为三个独立患者。详细状态见 `metadata/section_status.tsv`；blocked 的标签数量是 NA，不能解释为零细胞。

## 2. 输入恢复、来源与限制

原件和既有普通 Visium 结果只读。解压数据已位于 `highres_VisiumHD/data/`，旧位置的兼容链接不算重复数据。本轮只补回 Cervilla 的 probe_set.csv、aligned_fiducials.jpg、detected_tissue_image.jpg，共 **10,464,191 bytes（10.46 MB）**；校验了 ZIP 成员 CRC32 和独立 SHA256，没有宣称验证整份远端 ZIP 的 MD5。未下载 FASTQ/BAM 或重复整包。

[Cervilla Zenodo 17999961](https://zenodo.org/records/17999961) 的当前版本列出 processed archive；[作者预处理代码](https://github.com/scervilla/SpatialBenchmarking/blob/main/preprocessing/load_visiumhd.R) 使用 8 µm，并非原始细胞结果。Dataverse 镜像 DOI 10.34810/DATA3063 返回验证 HTML，未能枚举，故不声称镜像一定没有原图。[D1 的 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9066380) 及既有 188 文件 archive 索引未提供所需原始 QPTIFF；作者 Thor 推断/平滑表达不作为 raw counts。[UC_SCC 官方 GEO 页面](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9040999) 当前只列出 8 µm processed 文件。公开检索快照、时间、URL 和哈希见 `metadata/public_recovery_manifest.json`；本地定向检索见 `targeted_local_search.json`。

Cervilla 来源为 Space Ranger 3.0.0，Human Transcriptome Probe Set v2.0，Visium HD v1 FFPE。原始 positions 的 mpp=0.8636122748，hires 比例=0.52548605；现有 6000×4711 RGB 图像原生 **1.643454 µm/pixel**。使用提供的已配准坐标与比例，不按长宽猜配准；无附加翻转、旋转或裁剪。模型图放大到约 0.5 µm/pixel 仅便于模型运行，**不会恢复原始细胞结构**。变换含像素中心偏移，记录在 `metadata/Cervilla_transform.json`。barcode 完整匹配，所有最终细胞坐标已从原始 positions 独立复算一致。

## 3. 重建与宽松 QC

固定 seed=123。StarDist `2D_versatile_he`：按 RGB 通道 1/99.8 百分位归一化，prob_thresh=0.05，nms_thresh=0.3；整幅可用画布分块推断。代表性致密、间质和边缘区域保留 native 图、核边界和扩张叠图。Bin2cell 0.3.5 最大扩张 2 个 array bins（4 µm），k=4，PCA30 用于冲突 tie-breaking；destripe 调整只用于临时扩张表达视图。**最终 X 完全由未归一化原始整数 2 µm counts 按最终映射求和**。没有把调整后表达取整为 raw，也没有把未分配 counts 插补成细胞。

172,320 个核标签中，106,934 个获得组织内 bins，65,386 个未获得 bins，仅保存可追溯核几何信息，不计为重建细胞。未使用 GEX-only 补救：当前致密组织及有限原生图像细节不足以可靠确认额外补救对象。

| 计数 / QC 项目 | 数值 |
| --- | --- |
| 最低保留条件 | 原始 n_counts ≥10 且 n_genes ≥10，无已确认严重重建错误 |
| 保留 | 97,847 / 106,934（{97847/106934:.2%}） |
| 排除 | 9,087：9,086 低于表达底线；1 个确认图像边缘截断 |
| low_information | 12,798 / 97,847（{12798/97847:.2%}）；10–30 genes，全部保留 |
| 建立表达结构核心 | 84,985；>30 genes，初始技术/形态 MAD flags 均无 |
| 分配 bins / 未分配 bins | 3,039,039 / 1,683,872 |
| 全部原始 counts | 64,431,796 |
| 已分配 / 未分配 counts | 45,682,594 / 18,749,202 |
| 已分配比例 | {45682594/64431796:.2%}；未分配不能解释为表达缺失 |
| 保留 / 排除对象 counts | 45,632,235 / 50,359 |
| 保留细胞中位数 | counts=271；genes=227；bins=27 |

n_counts、n_genes、面积、核面积在 log1p 尺度上使用逐切片 5×1.4826 MAD 标记异常，不以面积/高 counts 单独删细胞，不使用额外 200/500 genes 门槛，也未使用 Scrublet/DoubletFinder。6 个最大核对象在原生图中呈多叶/浅染等不确定形态，无法区分合并核与大型异型核：保留并标记 mixed_or_ambiguous，记录逐对象图像框和理由；这些对象原已不在结构核心内。其他 marker 冲突是疑似混合，不宣称已证实 doublet。

mt fraction 仅代表实际检测到的 11 个 MT probe genes 的计数比例，覆盖不完整，不套用统一 mt% 删除阈值。各谱系/亚群 flag 与低信息分布见 `metadata/flags_by_final_identity.tsv`；被排除对象保留紧凑 metadata 和可重建映射，没有第二套排除矩阵。

## 4. 参考版本与分类规则

采用 **cell2location_v2_frozen_20260916__scRNA32**，标签字段 `reference_celltype_v2`。实际只读对象：`{ref['object_path']}`。该文件名含 30type，但其保留的 32 类 scRNA 标签与已冻结 barcode mapping 完全一致；后来 30 类空间合并报告仍标注 PENDING_REVIEW，未据此替换已确认的 scRNA 词表。已逐谱系核对 B/TNK/Myeloid/Fibroblast/Endothelial 冻结决策及 marker 框架，来源 SHA256 见 `source_provenance_checksums.json`。

参考 300,033 cells、17,987 genes、32 labels、5 datasets、58 orig.ident specimens。**58 个 specimens 不等于 58 个已核实 donor**。保守合并疑似同源/配对 specimen 为 source groups 后留出，供识别能力检查；不是已验证患者留出。源组 provenance 表只复用既有样本映射，不使用任何旧的 MP10 评分/结果。训练与留出 source groups、cell IDs 均无交叉。

基因优先 exact Ensembl，再用唯一 exact symbol 对齐；原始特征全部留在主对象中。建模共有 14,327 个参考基因，HD 的 14,329 个匹配 features 中有两组重复身份，只在建模视图求和。未测/未匹配基因不当作有检测能力的零表达。

自定义 CellTypist：46,140 个分层训练细胞，4,745 个留出细胞；按 source group×subtype 限额抽样，保留稀有群。稀疏 counts 在共有基因上 normalize_total=10,000 后 log1p；SGD，with_mean=False，alpha=1e-4，balanced class weights，max_iter=500，seed=123。保存 8-major、32-subtype、B/Plasma 谱系 4 类模型；没有把 scVI latent/插补表达作为 CellTypist 输入。保留 individual/top3、sigmoid score 和 top1−top2 logit margin；未使用空间多数投票。

接受规则是 major → 谱系内候选 → 实际观测 marker 支持。按类在匹配 HD 深度的参考下采样留出集上选择 logit margin：至少 20 个接受预测、经验 precision≥0.8；full-depth 可区分性仅作可行性门槛。普通细胞和低信息细胞分别用共有 counts 中位深度 359、21 校准。**同一留出集用于选阈值和描述性能，估计偏乐观；不能视为 HD 跨平台准确率或病理真值。** 模型分数本身也不是经验证准确率。

{table(validation)}

冻结前修正过一次深度域不匹配：不能把 full-depth 的较高 margin 阈值直接叠加到下采样 HD 条件；只改校准逻辑，未为追求某种生物学结果调整 QC 或重训模型。理由保存在 `pre_freeze_rule_correction.json`。所有最终结果使用修正后的匹配深度规则。

GC B 必须同时有 B major（≥2 个观测 B marker 且含 pan-B anchor）和 ≥2 个 BCL6/AICDA/RGS13 的 GC 区分 marker，并通过谱系分类校准。单一 BCL6、cycling、空间邻居、TLS 预期或代谢信号不能替代此证据。B 参考没有独立 Cycling B 类；IFN-stimulated B 被既有冻结决策排除，不在本次可预测词表内。

## 5. 最终标签与未解决问题

以下比例均以保留的 97,847 个细胞为分母：

{table(statuses)}

low_information 独立于 annotation_status：12,388 Unknown、409 major_only、1 mixed、0 subtype_supported。低置信与 Unknown 不被删除。`major_label` 在疑似 mixed 对象中可保留候选谱系，因此各大类数量应结合 annotation_status 读取，不能把混合对象当作已确定单一类型。

{table(support)}

主要可解释细分证据集中于 Matrix_CAF、Pericyte、SMC、iCAF。FOLR2 TAM 为 5 个细胞，仍需审阅；1–4 个细胞的稀有类别保留但视为孤立候选，不能直接推断稳定的细胞群体。C1QC TAM、vCAF、T/NK 细分、Naive/Memory B、GC B 等仍无法可靠区分。完整 32 类候选数、低信息比例、未解决比例和支持数见 `metadata/subtype_support_summary.tsv`；“无支持标签”不等于生物学不存在。

GC B 全谱系模型有 83 个候选，经过 B/Plasma 谱系模型后有 312 个层级候选，**0 个通过联合证据规则**。B 大类最终仅 2 个（保存 unresolved_B）。在全部保留对象的观测 counts 中，CD19 仅 29 个细胞检出、MS4A1 125、CD79A 276、AICDA 12、RGS13 44；BCL6 较广泛检出，但单独不足以判定 GC。T/NK major 仅 45，CD3D/CD3E 也极稀疏。模型候选往往缺少同一对象中的联合谱系标记，不能强行细分，也不能因一个 marker 未检出就宣称排除该谱系。检测率和原始计数总和见 `observed_marker_diagnostics.tsv`。

Epithelial、Plasma、Mast 的参考身份为 broad 类别，因此保存 major_only / unresolved_相应大类；**全部 epithelial 的 malignant_status 均未推断**。没有独立肿瘤状态证据。

## 6. scVI 的真实输入与检查

**scVI 输入仅 1 个 section：Cervilla。84,985 个核心细胞进入训练数据集（内部 train_size=0.9），97,847 个保留细胞获得 latent。** D1_Johnson 与 UC_SCC 未进入；Multi2 参考没有与 HD 全量联合训练。没有跨 section batch_key，也没有可声称的跨切片 batch correction 成果。

scvi-tools 1.5.1 / torch 2.8.0+cpu，2000 HVGs、NB likelihood、gene dispersion、latent20、hidden128、2 layers、batch512、lr=0.001、max_epochs50、early-stopping patience10、seed123。只保存最终权重，不保存逐 epoch checkpoints。日志报告 best validation ELBO 118.341；logger=False，未保存完整 epoch history，不能据此画出不存在的训练曲线。实际模型在用户提出 GPU 加速偏好前已用 8 个 CPU 线程完成；已记录今后必要训练优先评估 GPU，未为改变设备重复训练。

UMAP/表达近邻（k15）/Leiden（resolution0.6，7 clusters）只由核心集定义；低信息和其他非核心对象 model_role=inference_only，保留 latent，UMAP 为 NA。模型核心在后续 marker 冲突复核前确定，部分之后被标记为疑似混合的表达对象仍在训练结构中；它们不是已确认 doublet，不能作为 supported subtype 定义依据。

与 log-normalized PCA 基线进行同一组固定查询/参考的表达近邻检查：Epithelial 同 major 邻居比例 0.891→0.847、Fibroblast 0.748→0.789、Plasma 0.727→0.730、Endothelial 0.339→0.462；Myeloid 0.322→0.272。Mast 和 T/NK 在两种嵌入中聚集均弱，B 仅 2 个未作量化。此检查仅描述分类支持子集的结构，不是独立标签准确率，不足以证明稀有群保持或未来 niche 适用性。未据此删除细胞或利用邻居提高标签置信。数值见 `expression_structure_review.tsv`。

## 7. 验收、磁盘与复现

独立验收 **32 项通过**：原始三份输入矩阵 SHA256 未变；cell_id 唯一且 h5ad/标签表所有 obs 对齐；18,085 个输入特征完整；raw X 稀疏非负整数；按每个最终 bin-to-cell 映射从原始 2 µm 矩阵独立重聚合，**所有保留细胞×所有基因完全一致**，全部重建对象 QC、坐标和每基因计数守恒一致；低信息角色、Unknown 保留和参考留出隔离均验证。详见 `metadata/delivery_validation.tsv`。

清理仅限本轮可再生 pip 缓存和小型留出临时矩阵，释放约 {disk['removed_regenerable_bytes']/1024**2:.1f} MiB。交付数据、模型、代码、审计材料约 **{out/1024**2:.1f} MiB**；独立任务环境约 **{env/1024**3:.2f} GiB**；合计约 {(env+out)/1024**3:.2f} GiB（表观大小，最终 SHA 清单/报告略增）。未复制整套 Multi2，未保存多个全量细胞矩阵。共享文件系统资源与 GPU 快照见 `resources_final.json`，环境锁定清单见 `environment_versions.txt`。

主交付路径均相对于本报告目录 `{r}`：

- `cell_annotation.h5ad`：唯一主对象，X=原始整数计数，18,085 genes，97,847 cells，obs/var、raw 空间坐标、scVI/UMAP。
- `frozen_cell_labels.parquet`：全部保留细胞标签与质量/证据字段，cell_id 与主对象完全一致。
- `excluded_cell_metadata.parquet`：9,087 个排除对象及原因；`inputs/unassigned_nuclei.parquet` 保存未获 bins 的核标签。
- `inputs/Cervilla_bin_to_cell.parquet`、`Cervilla_nuclei.npz`、`Cervilla_nucleus_geometry.npz` 和 `metadata/Cervilla_transform.json`：重建映射、mask、几何与变换。
- `models/`：三个 CellTypist 模型、scVI 最终权重、StarDist 所用权重。
- `metadata/`：输入/参考 manifests、固定词表、gene maps、label_rules、marker 表、配置、公开检索、验收、SHA256 与 FROZEN marker。
- `code/` 与 `logs/`：模块化运行代码与执行日志。失败尝试日志保留为审计，最终分类来自 `classify_cells_final.log`。冻结保护阻止在 v1 目录直接重写；需要修改应建立 v2 并列出差异。

建议按顺序审阅以下图：

1. [QC 与排除原因](figures/Cervilla_QC_review.png)
2. [全图覆盖和未分配计数](figures/Cervilla_coverage_unassigned_review.png)
3. [代表性核分割及细胞扩张](figures/Cervilla_expansion_review.png)
4. [六个不确定大核对象](figures/Cervilla_largest_nuclei_review.png)
5. [表达 UMAP：sample、major、subtype、uncertainty](figures/Cervilla_expression_UMAP_review.png)
6. [全保留细胞空间标签](figures/Cervilla_spatial_labels_review.png)
7. [观测 marker 检测率与 GC B 候选](figures/Cervilla_observed_markers_GC_review.png)

分割和表达审核是计算层面的定性/定量复核；原图分辨率限制、局部核分割误差、未分配表达、跨平台参考偏差和稀有群不确定性仍存在。冻结的是带证据和这些限制的计算标签版本。**到此停止，等待用户审阅后再决定下一阶段。**

方法核对参考：[Bin2cell 官方流程](https://bin2cell.readthedocs.io/en/latest/notebooks/demo.html)、[10x segmented outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/segmented-outputs)、[CellTypist 官方训练 API](https://celltypist.readthedocs.io/en/latest/celltypist.train.html)、[scvi-tools 官方文档](https://docs.scvi-tools.org/en/stable/)。本项目 counts/genes≥10 的宽松底线是用户指定项目规则，不称为通用金标准。
'''
(r/'ANNOTATION_FREEZE_REPORT.md').write_text(report);print('REPORT_READY',len(report),disk,flush=True)
