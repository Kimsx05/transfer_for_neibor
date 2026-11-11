#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for s in Cervilla D1_Johnson UC_SCC MPBC1 MPBC2; do
 Rscript "$root/code/supplement_plots.R" "$root" "$s"
 if [[ "$s" != Cervilla && "$s" != D1_Johnson ]]; then Rscript "$root/code/plots.R" "$root" "$s"; fi
 Rscript "$root/code/recommended_plot.R" "$root" "$s"
done
