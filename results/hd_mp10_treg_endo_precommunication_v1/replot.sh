#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
O="$PWD/highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1"
PY="$PWD/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
for step in 13_replot 17_supplement_plots 20_final_figures; do "$PY" "$O/code/$step.py" > "$O/logs/$step.log" 2>&1; done
