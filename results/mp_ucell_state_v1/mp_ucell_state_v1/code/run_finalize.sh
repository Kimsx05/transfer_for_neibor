#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
while kill -0 642481 2>/dev/null; do sleep 5; done
test -f mp_ucell_state_v1/logs/FIGURES_COMPLETE
Rscript mp_ucell_state_v1/code/11_finalize_cluster_colors.R > mp_ucell_state_v1/logs/11_finalize_cluster_colors.log 2>&1
Rscript mp_ucell_state_v1/code/05_validate.R > mp_ucell_state_v1/logs/05_validate.log 2>&1
python3 mp_ucell_state_v1/code/10_validate_artifacts.py > mp_ucell_state_v1/logs/10_validate_artifacts.log 2>&1
python3 mp_ucell_state_v1/code/09_interpretation.py > mp_ucell_state_v1/logs/09_interpretation.log 2>&1
python3 mp_ucell_state_v1/code/06_report.py > mp_ucell_state_v1/logs/06_report.log 2>&1
