#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 RAYON_NUM_THREADS=1
py=../hd_epi_niche_multimethod_v2/.venv/bin/python
phase="${1:-verify}"
case "$phase" in
  inputs) "$py" code/audit.py ;;
  nuclei) bash code/gpu.sh cal; bash code/gpu.sh full; "$py" code/geometry.py; "$py" code/panel.py ;;
  model) "$py" code/classifier.py ;;
  native)
   bash code/install_binary.sh
   for s in Cervilla D1_Johnson; do
    if ! grep -q 'Exit status: 0' "$s/native_interface/native.log" 2>/dev/null; then "$py" code/proseg_native.py "$s"; fi
    if ! grep -q 'Exit status: 0' "$s/native_full/native.log" 2>/dev/null || { [[ ! -s "$s/native_full/registered/registration_audit.parquet" ]] && [[ ! -f "$s/native_full/mask.npy" ]]; }; then "$py" code/proseg_native.py full "$s"; fi
    if [[ ! -s "$s/native_full/registered/registration_audit.parquet" ]]; then "$py" code/repair_edge.py "$s"; fi
    bash code/native_downstream.sh "$s"
   done ;;
  conservative) "$py" code/identity.py; "$py" code/validate.py; "$py" code/cards.py both index; "$py" code/freeze.py; "$py" code/spatial.py; "$py" code/figures.py ;;
  finish|verify)
   export HD6_BRANCH=native
   if [[ "$phase" == finish ]]; then "$py" code/validate.py; "$py" code/freeze.py; "$py" code/spatial.py; "$py" code/lineage_report.py; "$py" code/figures.py; fi
   "$py" code/verify.py; "$py" code/summary.py; if [[ -d proseg_observed/D1_Johnson/card_archive ]]; then "$py" code/archive_check.py; "$py" code/final_disposition.py; fi ;;
  cards)
   export HD6_BRANCH=native
   for s in Cervilla D1_Johnson; do "$py" code/archive_cards.py "$s" 12; done ;;
  *) echo 'Phases: inputs, model, nuclei, native, conservative, finish, cards, verify. Visual notes must correspond to the current masks and model output.' >&2; exit 2 ;;
esac
