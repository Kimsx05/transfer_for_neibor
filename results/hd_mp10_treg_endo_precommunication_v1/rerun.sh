#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
O="$PWD/highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1"
PY="$PWD/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1 RCPP_PARALLEL_NUM_THREADS=1
export HD_RECOMPUTE=0
if [[ ${1:-} == --from-inputs ]]; then export HD_RECOMPUTE=1; fi
CPUS=$(python3 -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
ulimit -v 188743680
# Existing immutable source sizes/mtimes prevent accidental reuse against changed upstream inputs.
"$PY" - "$O" <<'PY'
from pathlib import Path
import pandas as pd,sys
p=Path(sys.argv[1])/'input_manifest.tsv'
if p.exists():
 for r in pd.read_csv(p,sep='\t').itertuples():
  f=Path(r.path);assert f.exists() and f.stat().st_size==r.size_bytes and f.stat().st_mtime_ns==r.mtime_ns,('Upstream input identity changed: use a documented new result version',r.role)
PY
run_py() { local target=$1 step=$2; if [[ $HD_RECOMPUTE == 1 || ! -f "$O/$target" ]]; then /usr/bin/time -v nice -n 10 "$PY" "$O/code/$step.py" > "$O/logs/$step.log" 2>&1; fi; }
run_r() { local target=$1 step=$2; shift 2; if [[ $HD_RECOMPUTE == 1 || ! -f "$O/$target" ]]; then /usr/bin/time -v taskset -c "$CPUS" nice -n 10 Rscript --vanilla "$O/code/$step.R" "$@" > "$O/logs/${step}_rerun.log" 2>&1; fi; }
# Recreate temporary sparse exports only if raw-stage outputs or required caches are missing.
if [[ $HD_RECOMPUTE == 1 || ! -f "$O/models/Cervilla_bin16_counts.npz" || ! -f "$O/models/fine_trained.rds" || ! -f "$O/models/fine_reference.rds" || ! -f "$O/tables/reference_holdout_predictions.tsv.gz" || ! -f "$O/tables/bin8_rctd_summary.json" || ! -f "$O/tables/bin16_rctd_summary.json" ]]; then
 "$PY" "$O/code/01_prepare.py" > "$O/logs/01_prepare.log" 2>&1
fi
run_r models/fine_reference.rds 02_reference
run_r models/fine_trained.rds 02_reference
run_r tables/reference_holdout_predictions.tsv.gz 02_reference
run_r tables/reference_target_marker_coverage.tsv 02c_reference_marker_audit
run_py tables/annotation_validation.tsv 02b_reference_validation
run_r tables/Cervilla_cell_predictions.tsv.gz 04_cells
# The native AUCell check needs bounded raw exports, intentionally removed after delivery.
if [[ ! -f "$O/tables/official_AUCell_validation.tsv" && ! -f "$O/temporary/cell_AUCell.mtx" && $HD_RECOMPUTE != 1 ]]; then
 /usr/bin/time -v nice -n 10 "$PY" "$O/code/05_measure.py" > "$O/logs/05_measure.log" 2>&1
fi
run_py tables/measurement_summary.tsv 05_measure
run_py tables/historical_score_validation.tsv 10b_historical_check
run_r tables/official_AUCell_validation.tsv 06_aucell_check
run_r tables/bin8_rctd_summary.json 03_rctd bin8
run_r tables/bin16_rctd_summary.json 03_rctd bin16
run_r tables/bin8_native_selected_weights.tsv.gz 03b_selected_weights
run_r tables/bin16_native_selected_weights.tsv.gz 03b_selected_weights
run_r tables/conservative4um_ROI_cells.tsv.gz 07b_conservative
run_py tables/D1_light_MP10_evaluability.tsv 07c_D1_light
run_py tables/Cervilla_cell_annotation.parquet 08_evaluate
# Unit choice is a documented frozen scientific decision; reruns reuse it, never choose by association P.
test -f "$O/tables/unit_selection.tsv"
run_py niche_frozen.json 09_niche
run_py tables/Cervilla_epithelial_regions.parquet 10_regional_scores
run_py tables/association_results.tsv 11_association
run_py tables/association_tissue_context.tsv 18_context_sensitivity
run_py tables/association_Endo_mixture_sensitivity.tsv 19_endothelial_mixture_sensitivity
"$PY" "$O/code/16_delivery_tables.py" > "$O/logs/16_delivery_tables.log" 2>&1
run_r tables/R_loader_validation.tsv 14_validate
run_py tables/validation_checks.tsv 15_validate
if [[ $HD_RECOMPUTE == 1 || ! -f "$O/figures/Cervilla_core_registered_overview.png" ]]; then bash "$O/replot.sh"; fi
"$PY" "$O/code/21_finalize.py" > "$O/logs/21_finalize.log" 2>&1
