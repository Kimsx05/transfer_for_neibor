# FROZEN_VERSION_GUARD
from pathlib import Path
import json,hashlib,datetime,ast,pandas as pd,anndata as ad
r=Path(__file__).resolve().parents[1]
if (r/'metadata/FROZEN.json').exists():raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
checks=pd.read_csv(r/'metadata/delivery_validation.tsv',sep='\t');assert len(checks)==32 and checks.passed.all()
assert (r/'metadata/delivery_validation.tsv').stat().st_mtime>(r/'cell_annotation.h5ad').stat().st_mtime
for p in (r/'code').glob('*.py'):ast.parse(p.read_text(),filename=str(p))
o=pd.read_parquet(r/'frozen_cell_labels.parquet');a=ad.read_h5ad(r/'cell_annotation.h5ad',backed='r');assert a.obs_names.equals(o.index)
assert (o.subtype_marker_genes_detected==o.observed_subtype_markers.map(lambda x:0 if not x else len(x.split(';')))).all()
# Freeze both the audited reference vocabulary and actual output-label semantics.
d=pd.read_csv(r/'metadata/reference_label_dictionary.tsv',sep='\t');pairs=[]
for label,g in o.groupby('subtype_label',observed=True):
 majors=sorted(g.major_label.astype(str).unique());pairs.append(dict(subtype_label=label,major_mapping=';'.join(majors),annotation_status=';'.join(sorted(g.annotation_status.astype(str).unique())),n=len(g),meaning='Uncertainty status; consult per-cell major_label' if label in ['Unknown','mixed_or_ambiguous'] else 'Broad lineage unresolved' if label.startswith('unresolved_') else 'Reference subtype supported by current acceptance rules'))
pd.DataFrame(pairs).to_csv(r/'metadata/frozen_label_dictionary.tsv',sep='\t',index=False);d[['reference_celltype_v2','major_label']].to_csv(r/'metadata/reference_subtype_to_major.tsv',sep='\t',index=False)
p=r/'metadata/label_rules.json';rules=json.loads(p.read_text());rules.update(major_acceptance='matching-depth per-class calibrated margin AND >=2 observed major markers; B additionally requires CD79A/MS4A1/CD19',fine_acceptance='major accepted AND candidate major agrees AND candidate matching-depth calibrated margin AND >=2 observed subtype markers AND >=1 specific anchor when defined; broad Epithelial/Plasma/Mast remain major_only',b_plasma_override='For major-model B or Plasma predictions use B/Plasma four-class candidate and its matching-depth calibration; retain global individual prediction separately',mixed_expression_rule='Epithelial >=3 distinct observed markers; other lineage >=3 distinct markers and >=5 marker counts; min(epithelial marker counts, other marker counts)/max >=0.5. Flag unconfirmed mixture/contamination, retain. Does not assert doublet.',morphology_override='Six reviewed large nuclear objects retained mixed_or_ambiguous; see morphology_review_objects.json',field_semantics='prediction_score/margin and second/third candidate refer to hierarchical prediction; individual_prediction_score is global model; global top3/margin remain in HD_predictions_subtype.parquet; observed marker evidence uses shared reference genes, excludes unmeasured/unmatched genes',freeze_note='Rules frozen before downstream analysis; later changes require v2.');p.write_text(json.dumps(rules,indent=2))
review=dict(reviewed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reviewer='computational_agent_visual_review_not_pathologist',native_resolution_um_per_pixel=1.6434542359679203,findings=['Representative native H&E shows recognizable nuclear outlines with residual missed/merged objects; enlarged image adds no structural detail.','Final expansion is limited to 4um from nuclear labels, retains unassigned bins, and follows the supplied in-tissue mask in reviewed ROIs.','Available histology canvas includes tissue outside the supplied captured in-tissue bin footprint; no expression inference is made for image-only regions.','Whole coverage review confirms counts are restricted to the supplied capture footprint; 1274 supplied in-tissue bins lie outside image bottom.','Largest six nuclei remain morphologically ambiguous, retained with explicit flags.','Expression UMAP shows broad epithelial/fibroblast/plasma patterns, but rare lineage separation limited; no confidence upgrades from clustering.','Marker heatmap uses observed raw counts only; GC candidates lack combined B/GC evidence.'],figures=['Cervilla_QC_review.png','Cervilla_coverage_unassigned_review.png','Cervilla_expansion_review.png','Cervilla_largest_nuclei_review.png','Cervilla_expression_UMAP_review.png','Cervilla_spatial_labels_review.png','Cervilla_observed_markers_GC_review.png'])
(r/'metadata/final_visual_review.json').write_text(json.dumps(review,indent=2))
p=r/'ANNOTATION_FREEZE_REPORT.md';report=p.read_text();needle='D1 的 2/8/16 µm 是同一张切片，未重复计数。';report=report.replace(needle,'图像画布还包含不在所提供捕获 bins 范围内的组织；没有对应 counts 的图像区域不被解释为零表达，也不生成细胞表达。'+needle);report=report.replace('label_rules、marker 表、配置','label_rules、最终标签词典 frozen_label_dictionary.tsv、reference_subtype_to_major.tsv、marker 表、配置');report=report.replace('最终分类来自 `classify_cells_final.log`。','最终分类来自 `classify_cells_final.log`。冻结前将 Unknown 候选的 subtype marker 数量字段改为对全部细胞独立统计实际检出数，避免把未评估值写成零；标签、raw counts、坐标不变，详见 `marker_field_semantics.json`。');p.write_text(report)
p=r/'metadata/run_config.json';cfg=json.loads(p.read_text());cfg['workflow_status']='frozen_computational_pending_user_review';p.write_text(json.dumps(cfg,indent=2))
for name in review['figures']:assert (r/'figures'/name).stat().st_size>1000
# Hash final immutable deliverables (exclude environment, pycache, active finalization log).
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
entries=[]
for p in sorted(r.rglob('*')):
 if not p.is_file() or p.is_symlink():continue
 rel=p.relative_to(r)
 if any(x in rel.parts for x in ['.venv','.pip_cache','__pycache__']):continue
 if p.name in ['deliverable_checksums.json','FROZEN.json','finalize_freeze.log']:continue
 entries.append(dict(path=str(rel),bytes=p.stat().st_size,sha256=sha(p)))
cp=r/'metadata/deliverable_checksums.json';cp.write_text(json.dumps(entries,indent=2));hashed={x['path']:x['sha256'] for x in entries}
marker=dict(annotation_version='frozen_subtypes_v1',status='computational_annotation_frozen_pending_user_review',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),workspace=str(r.parents[1]),output_directory=str(r),section_completion=dict(full_section=0,partial_ROI=['Cervilla'],blocked=['D1_Johnson','UC_SCC']),n_cells=len(o),n_genes=a.n_vars,annotation_counts=o.annotation_status.value_counts().to_dict(),scvi_sections=1,scvi_core_cells=int(o.scvi_training_core.sum()),reference_version='cell2location_v2_frozen_20260916__scRNA32',seed=123,validation_checks_passed=32,cell_matrix_sha256=hashed['cell_annotation.h5ad'],labels_sha256=hashed['frozen_cell_labels.parquet'],report_sha256=hashed['ANNOTATION_FREEZE_REPORT.md'],checksum_manifest_sha256=sha(cp),immutable_rule='No overwrite of v1. Changes require a new v2 and differences report.',next_action='Await user label review; no downstream scoring or niche analysis authorized in this run.')
with (r/'metadata/FROZEN.json').open('x') as f:json.dump(marker,f,indent=2)
print(json.dumps(marker,indent=2),flush=True)
