#!/usr/bin/env bash
set -euo pipefail
OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$OUT/../.."
export OUT PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1 RCPP_PARALLEL_NUM_THREADS=1 ARROW_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PY="$PWD/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
export R_LIBS="$PWD/highres_VisiumHD/analysis/hd_bin8_comparison_v1/Rlib:/refdir/Rlib"
export NUMBA_CACHE_DIR="$OUT/logs/numba" MPLCONFIGDIR="$OUT/logs/matplotlib"
cpus=$(python3 -c 'import os; print(",".join(map(str,sorted(os.sched_getaffinity(0))[:16])))')
taskset -pc "$cpus" $$ >/dev/null

