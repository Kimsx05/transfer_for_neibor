#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 RCPP_PARALLEL_NUM_THREADS=1
ulimit -v 185000000
cpu=$(awk '/Cpus_allowed_list/ {split($2,x,/[-,]/);print x[1]}' /proc/self/status)
exec /usr/bin/time -v -o "${V3_MODEL_OUT:-$PWD}/logs/time_composition.log" taskset -c "$cpu" Rscript --vanilla scripts/composition.R

