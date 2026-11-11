#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
hd=$(cd "$root/../.." && pwd)
complete=1
for s in UC_SCC MPBC1 MPBC2; do
 for f in l00_complete.json l02_complete.json l08_complete.json evidence_complete.json supplement_complete.json validation.json; do
  [[ -e "$root/$s/$f" ]] || complete=0
 done
done
if [[ $complete == 1 ]]; then exec bash "$root/code/postprocess.sh"; fi
exec "$hd/analysis/cell_annotation_freeze_v1/.venv/bin/python" "$root/code/night_gate.py" bash "$root/code/night_worker.sh"
