#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
plot_pid=$(cat mp_ucell_state_v1/logs/paused_plot_pid.txt)
trap 'kill -CONT "$plot_pid" 2>/dev/null || true' EXIT
Rscript mp_ucell_state_v1/code/04_figures_branch.R FULL_RAW > mp_ucell_state_v1/logs/04_figures_FULL_RAW.log 2>&1 &
pid_raw=$!
Rscript mp_ucell_state_v1/code/04_figures_branch.R BAL_Z_seed42 > mp_ucell_state_v1/logs/04_figures_BAL_Z_seed42.log 2>&1 &
pid_bz=$!
Rscript mp_ucell_state_v1/code/04_figures_branch.R BAL_RAW_seed42 > mp_ucell_state_v1/logs/04_figures_BAL_RAW_seed42.log 2>&1 &
pid_br=$!
status=0
wait "$pid_raw" || status=1
wait "$pid_bz" || status=1
wait "$pid_br" || status=1
exit "$status"
