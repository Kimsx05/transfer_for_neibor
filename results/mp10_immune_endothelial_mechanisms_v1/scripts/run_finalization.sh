#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_THREAD_LIMIT=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
for mp10_wait_pid in "$@"; do while kill -0 "$mp10_wait_pid" 2>/dev/null; do sleep 30; done; done
s='mp10_immune_endothelial_mechanisms_v1/scripts'
b='results/mp10_immune_endothelial_mechanisms_v1'
py='mp10_immune_endothelial_mechanisms_v1/venv/bin/python'
for lin in Epi T_lineage Endothelial; do
 "$py" "$s/29_verify_aucell.py" "$lin" > "$b/logs/scenic_${lin}_aucell_verification.log" 2>&1
 test -f "$b/tables/scenic_${lin}_connection_seed_stability.tsv.gz"
done
Rscript --vanilla "$s/23_regulatory_integration.R" > "$b/logs/23_regulatory_integration.log" 2>&1
Rscript --vanilla "$s/30_evidence_details.R" > "$b/logs/30_evidence_details.log" 2>&1
Rscript --vanilla "$s/32_axis_comparison.R" > "$b/logs/32_axis_comparison.log" 2>&1
"$py" "$s/35_reference_coverage.py" > "$b/logs/35_reference_coverage.log" 2>&1
Rscript --vanilla "$s/36_scenic_eligibility.R" > "$b/logs/36_scenic_eligibility.log" 2>&1
"$py" "$s/34_validate_scenic.py" > "$b/logs/34_validate_scenic.log" 2>&1
"$py" "$s/37_patient_relationships.py" > "$b/logs/37_patient_relationships.log" 2>&1
bash "$s/replot_all.sh" > "$b/logs/final_replot.log" 2>&1
Rscript --vanilla "$s/21_validate.R" > "$b/logs/21_validate.log" 2>&1
"$py" "$s/25_verify_inputs.py" > "$b/logs/25_verify_inputs.log" 2>&1
"$py" "$s/27_delivery_manifest.py" > "$b/logs/27_delivery_manifest.log" 2>&1
