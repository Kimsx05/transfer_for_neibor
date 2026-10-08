#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ -n "${1:-}" ]]; then
 while kill -0 "$1" 2>/dev/null; do sleep 5; done
fi
test -f mp_ucell_state_v1/02_ucell_scores/COMPLETE
Rscript mp_ucell_state_v1/code/02_graphs.R > mp_ucell_state_v1/logs/02_graphs.log 2>&1
Rscript mp_ucell_state_v1/code/03_summaries_compare.R > mp_ucell_state_v1/logs/03_summaries_compare.log 2>&1
Rscript mp_ucell_state_v1/code/04_figures.R > mp_ucell_state_v1/logs/04_figures.log 2>&1
