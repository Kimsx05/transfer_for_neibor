#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=8
ulimit -v 188743680
exec nice -n 10 Rscript "$@"
