#!/usr/bin/env bash
set -uo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
hd=$(cd "$root/../.." && pwd)
py="$hd/analysis/cell_annotation_freeze_v1/.venv/bin/python"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=8
run_section() {
 local section="$1"
 "$py" "$root/code/inputs.py" "$section" &&
 "$root/code/run_R.sh" "$root/code/prepare.R" "$root" "$section" &&
 "$py" "$root/code/check_cache.py" "$section" &&
 "$py" "$root/code/prepare_spatial.py" "$section" &&
 "$py" "$root/code/select_rois.py" "$section" &&
 "$root/code/run_R.sh" "$root/code/annotate.R" "$root" "$section" &&
 "$py" "$root/code/merge_annotations.py" "$section" &&
 "$py" "$root/code/coarse_types.py" "$section" &&
 "$root/code/run_R.sh" "$root/code/neighbors.R" "$root" "$section" &&
 "$root/code/run_R.sh" "$root/code/models.R" "$root" "$section" &&
 "$root/code/run_R.sh" "$root/code/evidence.R" "$root" "$section" &&
 "$py" "$root/code/summarize.py" "$section" &&
 "$py" "$root/code/interpret.py" "$section" &&
 "$root/code/run_R.sh" "$root/code/supplement.R" "$root" "$section" &&
 "$root/code/run_R.sh" "$root/code/plots.R" "$root" "$section" &&
 "$root/code/run_R.sh" "$root/code/supplement_plots.R" "$root" "$section" &&
 "$root/code/run_R.sh" "$root/code/validate.R" "$root" "$section"
}
status=0
"$root/code/run_R.sh" "$root/code/reference_panel_check.R" "$root" reference_check >> "$root/logs/reference_panel_check.log" 2>&1 || status=1
for section in UC_SCC MPBC1 MPBC2; do
 echo "SECTION_START $section $(date -Is)"
 if run_section "$section" >> "$root/logs/run_${section}.log" 2>&1; then
  echo "SECTION_COMPLETE $section $(date -Is)"
 else
  echo "SECTION_FAILED $section $(date -Is); continuing other sections";status=1
 fi
 "$py" "$root/code/compare.py" >> "$root/logs/comparison.log" 2>&1 || true
 "$py" "$root/code/report.py" >> "$root/logs/report.log" 2>&1 || true
done
bash "$root/code/postprocess.sh" || status=1
exit "$status"
