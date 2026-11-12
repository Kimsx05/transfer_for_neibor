#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
py="$root/.venv/bin/python"
"$py" "$root/code/check_cache.py"
for stage in zero_distance_review projection_diagnostics finalize final_summary semantic_comparison quality_plots comparison_figures review_rois composition_figures ordinary_correspondence; do "$py" "$root/code/$stage.py"; done
Rscript "$root/code/validate_loader.R" "$root"
"$py" "$root/code/report.py"
