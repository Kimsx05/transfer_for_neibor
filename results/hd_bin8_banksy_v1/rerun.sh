#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
hd=$(cd "$root/../.." && pwd)
py="$hd/analysis/cell_annotation_freeze_v1/.venv/bin/python"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=8
mkdir -p "$root/logs"
"$py" "$root/code/check_cache.py"
for section in Cervilla D1_Johnson; do
 "$py" "$root/code/prepare_spatial.py" "$section" >> "$root/logs/spatial_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/prepare.R" "$root" "$section" >> "$root/logs/prepare_${section}.log" 2>&1
 "$py" "$root/code/types_rois.py" "$section" >> "$root/logs/types_${section}.log" 2>&1
 "$py" "$root/code/conflict_rois.py" "$section" >> "$root/logs/conflict_${section}.log" 2>&1
 if [[ ! -f "$root/$section/l08_complete.json" || ! -f "$root/$section/l02_complete.json" || ! -f "$root/$section/l00_complete.json" ]]; then
  "$root/code/run_R.sh" "$root/code/neighbors.R" "$root" "$section" >> "$root/logs/neighbors_${section}.log" 2>&1
  "$root/code/run_R.sh" "$root/code/models.R" "$root" "$section" >> "$root/logs/models_${section}.log" 2>&1
 fi
 "$root/code/run_R.sh" "$root/code/evidence.R" "$root" "$section" >> "$root/logs/evidence_${section}.log" 2>&1
 "$py" "$root/code/summarize.py" "$section" >> "$root/logs/summary_${section}.log" 2>&1
 "$py" "$root/code/interpret.py" "$section" >> "$root/logs/interpret_${section}.log" 2>&1
 "$py" "$root/code/rctd_context.py" "$section" >> "$root/logs/rctd_context_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/plots.R" "$root" "$section" >> "$root/logs/plots_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/region_plots.R" "$root" "$section" >> "$root/logs/region_plots_${section}.log" 2>&1
 if [[ ! -f "$root/$section/validation.json" ]]; then
  "$root/code/run_R.sh" "$root/code/validate.R" "$root" "$section" >> "$root/logs/validation_${section}.log" 2>&1
 fi
done
"$py" "$root/code/report.py"
