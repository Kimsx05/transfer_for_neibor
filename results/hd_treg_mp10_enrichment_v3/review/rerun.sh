#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT=$(cd -- "$ROOT/../../.." && pwd)
PY="$ROOT/../hd_epi_niche_multimethod_v2/.venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1
cd -- "$PROJECT"
"$PY" "$ROOT/code/11_delivery.py" --check
hour=$(TZ=Asia/Shanghai date +%H)
if [[ ${HD_ALLOW_DAYTIME:-0} != 1 && $((10#$hour)) -ge 7 && $((10#$hour)) -lt 23 ]]; then
  echo 'Default compute window is23:00–07:00 Asia/Shanghai. Use HD_ALLOW_DAYTIME=1 only with existing direct-compute authorization.' >&2
  exit 2
fi
run_py(){ /usr/bin/time -v "$PY" "$ROOT/code/$1.py" > "$ROOT/logs/$1.log" 2>&1; }
run_py 01_audit
[[ -f "$ROOT/calibration_fix.yaml" && ${FORCE:-0} != 1 ]] || run_py 02_calibration_fix
[[ -f "$ROOT/MPBC2/neighborhoods.parquet" && ${FORCE:-0} != 1 ]] || run_py 03_geometry
run_py 04_enrichment
[[ -f "$ROOT/tables/mixing_scores.tsv" && ${FORCE:-0} != 1 ]] || run_py 05_mixing
[[ -f "$ROOT/tables/continuous_associations.tsv" && ${FORCE:-0} != 1 ]] || run_py 06_continuous
[[ -f "$ROOT/tables/standard_cross_KL.tsv" && ${FORCE:-0} != 1 ]] || run_py 07_standard_crossK
run_py 08_collect_validate
if [[ ! -f "$ROOT/tables/figure_index.tsv" || ${FORCE:-0} == 1 ]]; then run_py 09_plot; fi
if [[ ! -f "$ROOT/tables/label_measurement_diagnostics.tsv" || ${FORCE:-0} == 1 ]]; then run_py 09b_diagnostics; fi
if [[ ! -f "$ROOT/tables/figure_index_R.tsv" || ${FORCE:-0} == 1 ]]; then Rscript --vanilla "$ROOT/code/09c_distribution.R" > "$ROOT/logs/09c_distribution.log" 2>&1; fi
run_py 10_report
run_py 11_delivery
