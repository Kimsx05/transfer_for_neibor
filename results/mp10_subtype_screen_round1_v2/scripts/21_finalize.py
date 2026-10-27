from common import *
import subprocess,shutil,zipfile
# Join the independently rendered plot metadata, retaining actual node identities.
for kind in ['edges','nodes']:
 frames=[read(f'circle_{kind}_{mode}.tsv') for mode in ['primary','secondary']];save(pd.concat(frames,ignore_index=True),f'circle_{kind}_v2.tsv')
subprocess.run(['/usr/bin/python3',str(O/'scripts/19_report.py')],check=True)
subprocess.run(['/usr/bin/python3',str(O/'scripts/13_validate.py')],check=True)
# Small figure inventory covers the COMPLETE directory; the review copy is a curated subset.
files=sorted((O/'figures').rglob('*.pdf'));index=[]
for fp in files:
 rel=fp.relative_to(O);parts=rel.parts;ds='ALL' if parts[1]=='overall' else parts[2];index.append(dict(dataset=ds,name=fp.stem,path=str(rel),bytes=fp.stat().st_size))
save(index,'figure_inventory.tsv')
review=O/'review';review.mkdir(exist_ok=True)
for name in ['README.md','methods.md','round1_v2_report.md']:shutil.copy2(O/name,review/name)
keep_tables=['samples_v2.tsv','old_vs_new_samples.tsv','subtype_parent.tsv','frozen_MP10_Top20.tsv','score_validation.tsv','transform_parity.tsv','proportions_v2.tsv.gz','proportion_models_v2.tsv','model_sample_membership.tsv','B_marker_dataset_support.tsv','B_marker_dataset_gene_coverage.tsv','B_Plasma_target_coverage.tsv','main_candidate_evidence_v2.tsv','all_groups_evidence_summary.tsv','all_TME_state_LR_descriptive.tsv','selected_review_axes.tsv','B_Plasma_all_candidate_ranking.tsv','selected_axes_all_network_context.tsv','candidate_dataset_evidence.tsv','candidate_patient_differences.tsv','candidate_sample_expression.tsv.gz','candidate_gene_expression.tsv.gz','candidate_matched_coverage.tsv','matched_B_Plasma_LR_sensitivity.tsv','matched_B_Plasma_sensitivity_comparison.tsv','matched_B_Plasma_sensitivity_specimens.tsv','circle_edges_v2.tsv','circle_nodes_v2.tsv','network_plan.tsv','network_status.tsv','network_scaling.tsv','network_incoming_outgoing.tsv','network_cell_contributions.tsv.gz','matched_B_Plasma_network_status.tsv','stacked_bar_summary_v2.tsv','stacked_bar_boundaries_v2.tsv','existing_niche_B_Plasma_context.tsv','validation_checks.tsv','figure_inventory.tsv','package_versions.tsv','LR_definitions_v2.tsv','cell_colors.tsv']
(review/'tables').mkdir(exist_ok=True)
for name in keep_tables:shutil.copy2(T/name,review/'tables'/name)
# Each dataset receives the same figure selection, including reverse and NE panels.
core={'old_vs_new_eligibility','B_Plasma_MP10_points','B_composition_with_coverage','B_sensitivity_effects','B_marker_comparability','effects_major_TME','effects_subtype_parent','B_Plasma_LR_contrast_heatmap','B_Plasma_paired_expression','B_Plasma_patient_direction_summary','B_Plasma_common_specimen_sensitivity','major_incoming_outgoing_epi_state_within_specimen'}
for co in ['Normal_Tumor','T1_T2','sample_MP10_group_allEpi']:
 for pa in ['TME','B']:core.add('stacked_'+co+'_'+pa)
core.add('stacked_sample_MP10_group_allEpi_Atlas')
for co in ['epi_state_within_specimen','sample_MP10_group_allEpi','Normal_Tumor','T1_T2']:core.add('major_overview_'+co)
for axis in read('selected_review_axes.tsv').axis_id:core.add('circle_'+axis+'_epi_state_within_specimen')
for axis in ['AX01','AX02','AX03','AX04']:core.add('circle_'+axis+'_sample_MP10_group_allEpi')
for co in ['Normal_Tumor','T1_T2']:core.add('circle_AX01_'+co)
selected=[]
for fp in files:
 if fp.stem not in core:continue
 rel=fp.relative_to(O);dest=review/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(fp,dest);selected.append(rel)
lines=['# 逐队列核心图索引','', 'PDF均为透明背景；完整目录还保存300dpi PNG和其他条件/亚群图。所有队列使用相同选图规则，包括无支持、反向和NE结果。','', '顺序：覆盖与原始比例 → 模型与敏感性 → 主状态真实LR圈图 → 配对表达 → 次级条件背景。内部CellChat P不是患者差异P。','']
for ds in ['ALL']+sorted(read('samples_v2.tsv').dataset.unique()):
 lines+=['## '+ds,'']
 for rel in selected:
  parts=rel.parts;scope='ALL' if parts[1]=='overall' else parts[2]
  if scope==ds:lines.append(f'- [{rel.stem}]({rel.as_posix()})')
 lines.append('')
(review/'FIGURE_INDEX.md').write_text('\n'.join(lines)+'\n')
# Index links in report are portable within this review package.
report=(review/'round1_v2_report.md').read_text().replace('`review/FIGURE_INDEX.md`','`FIGURE_INDEX.md`');(review/'round1_v2_report.md').write_text(report)
(review/'README.md').write_text('# Round1 v2 轻量阅读包\n\n先读 [中文报告](round1_v2_report.md)，再看 [逐队列图索引](FIGURE_INDEX.md) 与 [方法](methods.md)。\n\n包含五个dataset及ALL的核心PDF、关键表。保留非显著、反向、单队列和NE结果；没有原始表达矩阵或Seurat/CellChat输入。完整PNG、数值网络和脚本在独立v2目录。\n')
zip_path=O/'review.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for fp in sorted(review.rglob('*')):
  if fp.is_file():z.write(fp,fp.relative_to(O))
print('FINAL',len(files),'PDF figures;',len(selected),'review PDFs;',round(zip_path.stat().st_size/1024**2,2),'MB zip')
(O/'logs/COMPLETED.txt').write_text('Modules 0–3 completed; v1 read-only; 68 main/secondary networks plus 6 matched-coverage sensitivity networks.\n')
