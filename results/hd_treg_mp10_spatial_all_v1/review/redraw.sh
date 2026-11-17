#!/usr/bin/env bash
set -euo pipefail
RESULT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PY="$RESULT_ROOT/../hd_epi_niche_multimethod_v2/.venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
"$PY" "$RESULT_ROOT/code/10_plot.py"
"$PY" "$RESULT_ROOT/code/13_aux_plots.py"
"$PY" "$RESULT_ROOT/code/14_report.py"
"$PY" "$RESULT_ROOT/code/16_package.py"
