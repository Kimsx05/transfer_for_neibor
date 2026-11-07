#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")" && pwd)
project=$(cd -- "$root/../.." && pwd)
py="$project/analysis/cell_annotation_freeze_v1/.venv/bin/python"
rr="$root/code/run_R.sh"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$py" "$root/code/cache_guard.py"
if [[ ! -f "$root/Cervilla/diagnosis.tsv.gz" || ! -f "$root/D1_Johnson/diagnosis.tsv.gz" ]]; then "$py" "$root/code/diagnose_rules.py" > "$root/logs/diagnose.log" 2>&1; fi
for ref in blca hpca; do
 if [[ ! -f "$root/reference/${ref}_trained.rds" ]]; then "$rr" "$root/code/train_${ref}.R" "$root" > "$root/logs/train_${ref}.log" 2>&1; fi
done
for section in Cervilla D1_Johnson; do
 if [[ ! -f "$root/$section/blca_prediction_manifest.json" || ! -f "$root/$section/hpca_prediction_manifest.json" ]]; then "$rr" "$root/code/predict_section.R" "$root" "$section" > "$root/logs/predict_${section}.log" 2>&1; fi
done
if [[ ! -f "$root/summary.tsv" ]]; then "$py" "$root/code/compare.py" > "$root/logs/compare.log" 2>&1; fi
if [[ ! -f "$root/D1_Johnson/unknown_candidate_summary.tsv" ]]; then "$py" "$root/code/extra_summary.py" > "$root/logs/extra_summary.log" 2>&1; fi
if [[ ! -f "$root/D1_Johnson/focused_cell_evidence.tsv" ]]; then "$py" "$root/code/review_evidence.py" > "$root/logs/review_evidence.log" 2>&1; fi
for section in Cervilla D1_Johnson; do
 if [[ ! -f "$root/$section/candidate_marker_statistics.tsv" ]]; then "$rr" "$root/code/marker_summary.R" "$root" "$section" > "$root/logs/markers_${section}.log" 2>&1; fi
 if [[ ! -f "$root/$section/focused_groups.tsv" ]]; then "$rr" "$root/code/focused_evidence.R" "$root" "$section" > "$root/logs/focused_${section}.log" 2>&1; fi
 if [[ ! -f "$root/$section/roi_comparison.tsv" ]]; then "$rr" "$root/code/plots.R" "$root" "$section" > "$root/logs/plots_${section}.log" 2>&1; fi
 if [[ ! -f "$root/$section/validation.json" ]]; then "$rr" "$root/code/validate.R" "$root" "$section" > "$root/logs/validate_${section}.log" 2>&1; fi
done
"$py" "$root/code/report.py"
