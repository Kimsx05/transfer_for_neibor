#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/run_env.sh"
bash "$OUT/code/fetch_sources.sh"
py(){ "$PY" "$OUT/code/$1.py" "${@:2}"; }
r(){ Rscript --vanilla "$OUT/code/$1.R" "${@:2}"; }
test -s "$OUT/shared/taxonomy.tsv" || py xe_prepare
py inputs
py frozen_truth
test -s "$OUT/hd/UC_SCC/conservation.json" || py hd_prepare
for s in Cervilla D1_Johnson UC_SCC A1R2 A1R4 A1R5 A1R6 A2R2; do
 if [[ -d "$OUT/hd/$s" ]]; then p="$OUT/hd/$s"; target=epi64_centers.parquet; else p="$OUT/xenium/$s"; target=graph_nodes.parquet; fi
 test -s "$p/$target" || py geometry "$s"
done
py summarize_units frozen
for s in A2R2 A1R2 A1R4 A1R5 A1R6; do
 test -s "$OUT/logs/${s}_R_session.txt" || r xenium_fit "$s"
 py summarize_units "$s"; py review "$s"; py spatial_analysis "$s"
done
py boundary_audit; py boundary_sensitivity; py unknown_audit; r source_diagnostic; py ucell_check; py validation_prepare; r validation; py validation_metrics; r violins
for s in Cervilla D1_Johnson UC_SCC; do
 test -s "$OUT/hd/$s/fit_complete.txt" || r hd_fit "$s"
 py summarize_units "$s"; py spatial_analysis "$s"
done
test -s "$OUT/hd/Cervilla/scale_window.tsv" || py scale_check
py scale_summary; py report; py final_checks

