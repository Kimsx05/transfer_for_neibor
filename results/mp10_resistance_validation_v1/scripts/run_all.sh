#!/usr/bin/env bash
set -euo pipefail
B='/home/data/t070721/codex_workspace/Bladder Metabolism/results/mp10_resistance_validation_v1'
cd "$B"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -v 12582912
stamp=$(date +%Y%m%dT%H%M%S)
logdir="$B/logs/replay_$stamp"
mkdir -p "$logdir"
/usr/bin/time -v Rscript scripts/00_extract_clinical.R > "$logdir/00_extract_clinical.log" 2>&1
/usr/bin/time -v /usr/bin/python3 scripts/01_prepare.py > "$logdir/01_prepare.log" 2>&1
/usr/bin/time -v Rscript scripts/02_bulk_scores.R > "$logdir/02_bulk_scores.log" 2>&1
/usr/bin/time -v Rscript scripts/03_epi_mapping.R > "$logdir/03_epi_mapping.log" 2>&1
/usr/bin/time -v /usr/bin/python3 scripts/04_epi.py > "$logdir/04_epi.log" 2>&1
/usr/bin/time -v /usr/bin/python3 scripts/04b_patient_summary.py > "$logdir/04b_patient_summary.log" 2>&1
/usr/bin/time -v Rscript scripts/05_models.R > "$logdir/05_models.log" 2>&1
/usr/bin/time -v Rscript scripts/06_figures.R > "$logdir/06_figures.log" 2>&1
/usr/bin/python3 scripts/08_sources.py > "$logdir/08_sources.log" 2>&1
/usr/bin/python3 scripts/07_validate.py > "$logdir/07_validate.log" 2>&1
/usr/bin/python3 scripts/07_finalize.py > "$logdir/07_finalize.log" 2>&1
printf '%s\n' "$logdir" > "$B/logs/latest_replay.txt"
printf 'Replay completed: %s\n' "$logdir"
