from base import *
from zipfile import ZipFile,ZIP_DEFLATED
from load_result import load_unit
from source import kinds
for s in sections:
 for k in kinds:load_unit(s,k,with_counts=True)
paths=[R/n for n in ['REPORT.md','METHODS.md','FIGURE_INDEX.md','rerun.sh','input_manifest.tsv','roi_manifest.tsv','method_config.yaml','reference_validation.tsv','confusion.tsv','reference_internal_cv.tsv','historical_unit_baselines.tsv','holdout_mixture_predictions.tsv','unit_comparison.tsv','cell_assignment_audit.tsv','FOXP3_capture_localization.tsv','candidate_marker_summary.tsv','source_overlap_and_coverage.tsv','validation.tsv','versions.tsv','figures/ROI_atlas.pdf','figures/reference_confusion.png','logs/model_status.json']]
paths+=list((R/'code').glob('*.py'))+[R/s/'count_accounting.tsv' for s in sections]
with ZipFile(R/'review.zip','w',ZIP_DEFLATED,compresslevel=6)as z:
 for p in paths:z.write(p,p.relative_to(R))
with ZipFile(R/'review.zip')as z:assert z.testzip() is None
print('archive validated',len(paths),'files', (R/'review.zip').stat().st_size,'bytes')
