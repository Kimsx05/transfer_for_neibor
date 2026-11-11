from pathlib import Path
import pandas as pd,json,zipfile
from datetime import datetime
from zoneinfo import ZoneInfo
r=Path(__file__).resolve().parents[1];tabs=r/'tables'
def table(d):
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))
status=pd.read_csv(tabs/'completion_status.tsv',sep='\t');qc=pd.read_csv(tabs/'input_QC_summary.tsv',sep='\t');now=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(timespec='seconds')
lines=['# 8 μm BANKSY extension: status, evidence and limits','',f'Updated: {now} (Asia/Shanghai)','', 'All executable branches are complete and reread-validated. Original two-section results are reused; three new sections have full-bin results. Missing tissue/image inputs remain explicit limitations.','',table(status),'','## Input and QC','',table(qc),'','MPBC uses supplied-expression QC only; original tissue membership remains NA. It cannot be pooled as a tissue-verified denominator with Cervilla/D1/UC_SCC. Original matrices and old official results are preserved.','', '## D1 supplementary evidence','']
synthesis=r/'cross_section_review.md'
if synthesis.exists():lines[6:6]=['## 本轮结论','',synthesis.read_text(),'','## 完成状态','']
p=r/'D1_Johnson/supplement_marker_statistics.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t',dtype={'group':str});z=d[(d.context=='l08')&d.group.isin(['1','2'])][['group','symbol','pct_positive','mean_log','log2_mean_effect']].copy();z['pct_positive']=z.pct_positive.map(lambda x:f'{x:.2%}');z['mean_log']=z.mean_log.round(3);z['log2_mean_effect']=z.log2_mean_effect.round(3);lines +=[table(z),'','These raw-expression measurements supplement the old epithelial/collagen conflict interpretation. Keratin detection and relative enrichment are separated; collagen detection does not by itself prove fibroblast identity, cell admixture, EMT or malignancy. Old region labels remain unchanged.']
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 pre=r/s/'input_visual_review.md'
 if pre.exists():lines+=['','### '+s+' input review','',pre.read_text()]
 p=r/s/'visual_review.md';lines+=['','## '+s,'',p.read_text() if p.exists() else 'Visual evidence review pending; do not treat automated cluster summaries as a completed histological interpretation.']
p=tabs/'reference_disagreement_strata.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');z=d[(d.information_band=='All')&(d.denominator=='both_candidates')].copy();z=z[['section','level','n','disagree_fraction','disagree_UMI_fraction']];z[['disagree_fraction','disagree_UMI_fraction']]=z[['disagree_fraction','disagree_UMI_fraction']].map(lambda x:f'{x:.2%}');lines+=['','## Same-denominator reference comparison','',table(z),'','Per-information-band bin/UMI comparisons are in `tables/reference_disagreement_strata.tsv`; expanded raw-panel support is in `tables/information_strata.tsv` and `tables/reference_conflict_strata.tsv`. `information_standardized_sensitivity.tsv` gives a single bounded common-information weighted comparison; no models/labels are changed. Histological composition and assay/input differences prevent ranking sample quality from these percentages alone.']
p=r/'cross_section_review.md'
lines+=['','## Cross-section interpretation','','见报告开头的本轮结论；分母和分层原表随交付保留。','', '## Reading and rerun','', '```r',f'root <- "{r}"','source(file.path(root,"code/load_result.R"))','x <- load_extension(root,"UC_SCC") # or MPBC1, MPBC2, Cervilla, D1_Johnson','```','', 'The loader attaches the completed incremental results by unique bin ID. Final object validation is recorded separately in each `validation.json`. New sections reference native UC_SCC H5 or their sole `counts8.h5`; `normalized.rds` is the sole normalized layer. Old sections depend on sibling `hd_bin8_comparison_v1` and `hd_bin8_banksy_v1`.','', '`bash rerun.sh` reuses completed models and runs light postprocessing directly. If model stages are missing, it uses the night supervisor (23:00–07:00 CST), processes UC_SCC then MPBC1/2, resumes caches, and updates comparisons/report. Active heavy child groups pause outside that window. Inspect `logs/night_status.json`, `logs/night_supervisor.log` and `logs/run_SECTION.log`. No stage-complete claim is made merely because a job is queued.','', 'Required interpretation tables: `cluster_domain_interpretation.tsv`, `cluster_summary.tsv`, `cluster_markers.tsv.gz`, `panel_marker_statistics.tsv.gz`, `cluster_SingleR.tsv`; predictions retain both references, pruning and scoring. Supplementary 14-gene evidence is in `supplement_marker_statistics.tsv` and `supplement_gene_coverage.tsv`.','', 'All figures are white300dpi PNG plus PDF/legend/panel variants where appropriate. Review archive selects key PNGs and evidence tables; all complete figure formats remain on server. No expression/reference/large embedding data are in review.zip.','', 'This round does not perform MP10, subtyping, CNV/malignancy, communication, formal neighborhood ecology, new segmentation or RCTD.','']
report=r/'REPORT.md';report.write_text('\n'.join(lines)+'\n')
(r/'WORK_STATUS.md').write_text('# Actual stage status\n\nUpdated: '+now+'\n\n'+table(status)+'\n\nRunning state is recorded in logs/night_status.json; completion requires the corresponding results and readback validation, not a queued script.\n')
small_names={'input_summary.json','input_load_validation.json','input_check.json','annotation_summary.json','blca_prediction_summary.json','hpca_prediction_summary.json','supplement_complete.json','supplement_gene_coverage.tsv','supplement_marker_statistics.tsv','fixed_rois.tsv','neighbor_summary.json','spatial_input_summary.json','spatial_components.tsv','cluster_summary.tsv','cluster_reference_composition.tsv','cluster_domain_interpretation.tsv','cluster_SingleR.tsv','cluster_SingleR_scores.tsv','cluster_markers.tsv.gz','panel_marker_statistics.tsv.gz','model_versions.json','sessionInfo.txt','validation.json','transform.json','visual_review.md','input_visual_review.md','spatial_structure_metrics.tsv','conflict_information_summary.tsv','coarsening_comparison.tsv','coarse_composition.tsv','coarse_crosstabs.tsv','gene_mapping.tsv','ROI_hires_panel_identity.tsv','geometry_decision.md'}
core={'ROI_single_genes_conflict_focus','spatial_expanded_original_panels','input_information_all_supplied','spatial_coarse_references','spatial_interpreted_three_lambda','ROI_clusters_typical_information','ROI_single_genes_typical_information','information_composition','reference_disagreement_by_information','raw_panel_codetection_by_information','recommended_reviewed_regions','ROI_hires_fixed'}
def files():
 for p in r.rglob('*'):
  if not p.is_file():continue
  rel=p.relative_to(r)
  if any(x in rel.parts for x in ['tmp','__pycache__','logs']):continue
  if len(rel.parts)==1 and p.name!='review.zip':yield p
  elif rel.parts[0]=='code' or rel.parts[0]=='tables':yield p
  elif p.name in small_names:yield p
  elif 'figures' in rel.parts and p.suffix=='.png' and p.stem in core:yield p
def pack():
 with zipfile.ZipFile(r/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files():z.write(p,p.relative_to(r))
pack();size=sum(p.stat().st_size for p in r.rglob('*') if p.is_file());archive=(r/'review.zip').stat().st_size;allocated=sum(p.stat().st_blocks*512 for p in r.rglob('*') if p.is_file());report.write_text(report.read_text()+f'\nIncremental files: {size/1024**3:.3f} GiB logical, {allocated/1024**3:.3f} GiB allocated, including review.zip ({archive/1024**2:.1f} MiB). This value includes current temporary caches, if any; it is updated after rerun.\n');pack();print('REPORT_UPDATED',now)
