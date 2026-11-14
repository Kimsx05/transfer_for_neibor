#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../../.."
O="$PWD/highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1 RCPP_PARALLEL_NUM_THREADS=1
export R_PROFILE_USER="$O/code/worker_profile.R"
PY="$PWD/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
CPUS=$(python3 -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
ulimit -v 188743680
for step in 02_reference 04_cells; do
  if [ "$step" = "02_reference" ] && [ -f "$O/tables/reference_holdout_predictions.tsv.gz" ]; then continue; fi
  if [ "$step" = "04_cells" ] && [ -f "$O/tables/Cervilla_cell_predictions.tsv.gz" ]; then continue; fi
  /usr/bin/time -v taskset -c "$CPUS" nice -n 10 Rscript --vanilla "$O/code/$step.R" > "$O/logs/$step.log" 2>&1; done
if [ ! -f "$O/tables/measurement_summary.tsv" ]; then /usr/bin/time -v nice -n 10 "$PY" "$O/code/05_measure.py" > "$O/logs/05_measure.log" 2>&1; fi
/usr/bin/time -v taskset -c "$CPUS" nice -n 10 Rscript --vanilla "$O/code/06_aucell_check.R" > "$O/logs/06_aucell_check.log" 2>&1
for unit in bin8 bin16; do
  if [ -f "$O/tables/${unit}_rctd_summary.json" ]; then continue; fi
  /usr/bin/time -v taskset -c "$CPUS" nice -n 10 Rscript --vanilla "$O/code/03_rctd.R" "$unit" > "$O/logs/03_rctd_$unit.log" 2>&1; done
