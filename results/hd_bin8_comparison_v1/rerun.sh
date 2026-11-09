#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
hd=$(cd "$root/../.." && pwd)
py="$hd/analysis/cell_annotation_freeze_v1/.venv/bin/python"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=8
cd "$hd"
# Sequential sections; fixed project parameters. Retains original cell results.
if [[ ! -f "$root/Cervilla/native8/spatial/tissue_positions.parquet" ]]; then
 "$py" "$root/code/recover_native8.py" >> "$root/logs/native8_recovery.log" 2>&1
fi
for section in Cervilla D1_Johnson; do
 "$py" "$root/code/prepare.py" "$section" >> "$root/logs/prepare_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/structure.R" "$root" "$section" >> "$root/logs/structure_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/annotate.R" "$root" "$section" >> "$root/logs/annotate_${section}.log" 2>&1
 "$py" "$root/code/roi_expression.py" "$section" >> "$root/logs/roi_expression_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/rctd.R" "$root" "$section" >> "$root/logs/rctd_${section}.log" 2>&1
 if [[ "$section" == Cervilla ]]; then
  "$root/code/run_R.sh" "$root/code/author_sensitivity.R" "$root" "$section" >> "$root/logs/author_sensitivity.log" 2>&1
 fi
 if [[ ! -f "$root/$section/validation.json" ]]; then
  "$root/code/run_R.sh" "$root/code/finalize.R" "$root" "$section" >> "$root/logs/finalize_${section}.log" 2>&1
 fi
 "$py" "$root/code/compare.py" "$section" >> "$root/logs/compare_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/plots.R" "$root" "$section" >> "$root/logs/plots_${section}.log" 2>&1
 "$root/code/run_R.sh" "$root/code/coverage_plots.R" "$root" "$section" >> "$root/logs/coverage_${section}.log" 2>&1
done
"$py" "$root/code/report.py"
