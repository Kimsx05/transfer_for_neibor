#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
hd=$(cd "$root/../.." && pwd)
py="$hd/analysis/cell_annotation_freeze_v1/.venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in UC_SCC MPBC1 MPBC2; do "$py" "$root/code/check_cache.py" "$s"; done
"$py" "$root/code/attach_review.py"
"$py" "$root/code/compare.py"
"$py" "$root/code/final_tables.py"
Rscript "$root/code/reference_panel_check.R" "$root" reference_check
Rscript "$root/code/comparison_plots.R" "$root" comparison
Rscript "$root/code/uc_roi_histology.R" "$root" UC_SCC
for s in Cervilla D1_Johnson UC_SCC MPBC1 MPBC2; do
 if [[ -e "$root/$s/supplement_complete.json" ]]; then
  Rscript "$root/code/recommended_plot.R" "$root" "$s"
  Rscript "$root/code/validate.R" "$root" "$s"
 fi
done
"$py" "$root/code/report.py"
