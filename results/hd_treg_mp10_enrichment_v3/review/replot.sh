#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT=$(cd -- "$ROOT/../../.." && pwd)
PY="$ROOT/../hd_epi_niche_multimethod_v2/.venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1
cd -- "$PROJECT"
"$PY" "$ROOT/code/11_delivery.py" --check
unset PLOT_UPDATE_NAMES
"$PY" "$ROOT/code/09_plot.py"
"$PY" "$ROOT/code/09b_diagnostics.py"
Rscript --vanilla "$ROOT/code/09c_distribution.R"
"$PY" "$ROOT/code/10_report.py"
"$PY" "$ROOT/code/11_delivery.py"
