#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=4
ulimit -v 167772160
exec nice -n 10 Rscript "$@"
