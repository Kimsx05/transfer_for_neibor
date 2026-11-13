#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
lin="$1"
scenic_wait_pid="${2:-}"
if [ -n "$scenic_wait_pid" ]; then
  while kill -0 "$scenic_wait_pid" 2>/dev/null; do sleep 30; done
fi
scenic_script_root='mp10_immune_endothelial_mechanisms_v1/scripts'
scenic_log_root='results/mp10_immune_endothelial_mechanisms_v1/logs'
scenic_python='mp10_immune_endothelial_mechanisms_v1/venv/bin/python'
"$scenic_python" "$scenic_script_root/24_grn_resume.py" "$lin" --adopt-initial-completed > "$scenic_log_root/scenic_${lin}_grn_provenance.log" 2>&1
"$scenic_python" "$scenic_script_root/13_scenic_ctx.py" "$lin" > "$scenic_log_root/scenic_${lin}_ctx.log" 2>&1
if ! Rscript --vanilla "$scenic_script_root/14_aucell_full.R" "$lin" > "$scenic_log_root/scenic_${lin}_aucell.log" 2>&1; then
  echo "AUCell R process exited nonzero; require full output verification before continuing" >> "$scenic_log_root/scenic_${lin}_aucell.log"
fi
"$scenic_python" "$scenic_script_root/29_verify_aucell.py" "$lin" > "$scenic_log_root/scenic_${lin}_aucell_verification.log" 2>&1
Rscript --vanilla "$scenic_script_root/15_scenic_summarize.R" "$lin" > "$scenic_log_root/scenic_${lin}_summary.log" 2>&1
if [ "$lin" = Epi ]; then Rscript --vanilla "$scenic_script_root/15_scenic_summarize.R" Epi M20-masked > "$scenic_log_root/scenic_Epi_masked_summary.log" 2>&1; fi
"$scenic_python" "$scenic_script_root/15b_scenic_rss.py" "$lin" > "$scenic_log_root/scenic_${lin}_rss.log" 2>&1
"$scenic_python" "$scenic_script_root/16_scenic_second_seed.py" "$lin" > "$scenic_log_root/scenic_${lin}_seed456.log" 2>&1
