#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
out="results/mp10_p95_front_programs_v1"
for i in 1 2 3 4 5 6; do Rscript "$out/scripts/score.R" "$out" "$i" 6 > "$out/logs/score_$i.log" 2>&1 & done
wait
