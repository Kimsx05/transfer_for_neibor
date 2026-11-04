#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=results/mp10_section_specificity_pearson_v1
for i in 1 2 3 4 5 6; do Rscript "$out/scripts/score.R" "$i" 6 > "$out/logs/score_$i.log" 2>&1 & done
wait
