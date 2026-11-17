#!/usr/bin/env bash
set -euo pipefail
RESULT_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd -- "$RESULT_ROOT/../../.." && pwd)
PY="$RESULT_ROOT/../hd_epi_niche_multimethod_v2/.venv/bin/python"
cd -- "$PROJECT_ROOT"
export HD_ALLOW_DAYTIME=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1 FORCE="${FORCE:-0}"
# CPU affinity applies to all descendant R/Python workers as a total budget.
if [[ ${HD_AFFINITY_SET:-0} != 1 ]]; then
  cpus=$("$PY" -c 'import os; print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
  export HD_AFFINITY_SET=1
  exec taskset -c "$cpus" bash "$RESULT_ROOT/rerun.sh"
fi
run(){ /usr/bin/time -v "$PY" "$RESULT_ROOT/code/$1.py" "${@:2}" > "$RESULT_ROOT/logs/rerun_$1.log" 2>&1; }
need(){ [[ "$FORCE" == 1 || ! -f "$RESULT_ROOT/$1" ]]; }
# Only identity metadata, not full checksums or full-matrix verification on every resume.
"$PY" - "$RESULT_ROOT" <<'PY'
import json,sys
from pathlib import Path
r=Path(sys.argv[1])
for p in r.glob('*/prepared.json'):
 d=json.loads(p.read_text());f=Path(d['counts_path']);assert f.exists() and f.stat().st_size==d['counts_bytes'] and f.stat().st_mtime_ns==d['counts_mtime_ns'],f
PY
for s in Cervilla D1_Johnson UC_SCC MPBC1 MPBC2; do
  if need "$s/prepared.json"; then run 01_prepare "$s"; fi
done
if [[ ! -f "$RESULT_ROOT/models/reference_prepared.json" ]]; then run 02_reference; fi
if [[ ! -f "$RESULT_ROOT/frozen_config.yaml" ]]; then run 03_calibrate_freeze; fi
if need 'Cervilla/T32_regions.parquet' || [[ ! -f "$RESULT_ROOT/MPBC2/Epi64_regions.parquet" ]]; then run 15_restore_frozen_scores; fi
if need 'tables/targeted_T_depth_validation.tsv'; then run 04b_targeted_depth; fi
if need 'tables/native_score_validation.tsv'; then
 run 04_export_ranks
 /usr/bin/time -v Rscript --vanilla "$RESULT_ROOT/code/04_native_UCell.R" > "$RESULT_ROOT/logs/rerun_native_UCell.log" 2>&1
 run 09_native_validation
fi
if [[ ! -f "$RESULT_ROOT/tables/roi_manifest.tsv" ]]; then run 04c_freeze_roi; fi
for s in Cervilla D1_Johnson UC_SCC MPBC1 MPBC2; do
 if need "$s/T32_neighborhoods.parquet"; then run 05_neighbors "$s"; fi
done
if need 'tables/association_all_results.tsv'; then run 06_associations; fi
if need 'tables/spottedpy_hotspot_distances.tsv.gz'; then run 07_candidates_spotted; fi
if need 'tables/cell_region_association_comparison.tsv'; then run 08_cell_comparison; fi
if need 'tables/cell_bin_common_effects.tsv'; then run 08b_common_cell_comparison; fi
if need 'tables/validation_checks.tsv'; then run 11_validate; fi
run 12_delivery_tables
if need 'figures/MPBC2_control_2_atlas.png'; then run 10_plot; fi
if need 'figures/ROI_raw_marker_evidence.png'; then run 13_aux_plots; fi
run 14_report
run 16_package
printf '%s\n' "Completed/resumed: $RESULT_ROOT"
