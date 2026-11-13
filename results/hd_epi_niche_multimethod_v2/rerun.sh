#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
exec 9>"$root/logs/run.lock"
flock -n 9 || { echo 'An execution already holds the task lock'; exit 1; }
if [[ "${1:-}" == "--report-only" ]]; then
  "$root/.venv/bin/python" "$root/code/report.py"
  exit 0
fi
"$root/.venv/bin/python" "$root/code/freeze.py"
"$root/.venv/bin/python" "$root/code/check_inputs.py"
"$root/.venv/bin/python" "$root/code/prepare_roi_links.py"
"$root/.venv/bin/python" "$root/code/night_supervise.py"
