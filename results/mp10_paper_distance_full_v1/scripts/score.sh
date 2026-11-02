#!/usr/bin/env bash
set -eu
out=${MP10_OUTPUT:-$(dirname "$(dirname "$(realpath "$0")")")}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
workers=${SCORE_WORKERS:-6}
pids=()
for ((i=1;i<=workers;i++)); do
 Rscript "$(dirname "$0")/score.R" "$out" "$i" "$workers" > "$out/logs/score_$i.log" 2>&1 &
 pids+=($!)
done
status=0
for child in "${pids[@]}"; do wait "$child" || status=1; done
exit "$status"
