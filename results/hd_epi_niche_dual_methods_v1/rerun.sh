#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=4 NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py="$root/.venv/bin/python"
"$py" "$root/code/check_cache.py"
# Active model jobs are not duplicated.
if pgrep -f "python.*$root/code/model.py" >/dev/null; then echo 'Existing model process present; inspect logs/runtime.json'; exit 1; fi
for stage in prepare support coverage_geometry compare_composition; do "$py" "$root/code/$stage.py"; done
need_models=0
for ref in blca hpca; do for r in 50 200; do [[ -f "$root/${ref}_r${r}/complete.json" ]] || need_models=1; done; done
if (( need_models )); then
hour=$(TZ=Asia/Shanghai date +%H)
if (( 10#$hour >= 7 && 10#$hour < 23 )); then echo 'Heavy fitting waits until 23:00 Asia/Shanghai'; while :; do hour=$(TZ=Asia/Shanghai date +%H); if (( 10#$hour >= 23 || 10#$hour < 7 )); then break; fi; sleep 30; done; fi
"$py" "$root/code/model.py" &
model_pid=$!
"$py" "$root/code/night_watch.py" "$model_pid" &
wait "$model_pid"
fi
for stage in zero_distance_review projection_diagnostics finalize final_summary semantic_comparison quality_plots comparison_figures review_rois composition_figures ordinary_correspondence; do "$py" "$root/code/$stage.py"; done
if [[ ! -s "$root/tables/loader_validation.tsv" ]]; then Rscript "$root/code/validate_loader.R" "$root"; fi
"$py" "$root/code/report.py"
