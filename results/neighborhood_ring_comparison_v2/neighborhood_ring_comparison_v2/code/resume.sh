#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/neighborhood_ring_comparison_v2'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/usr/bin/python3.12 "$ROOT/code/01_prepare.py"
/usr/bin/python3.12 "$ROOT/code/02_cluster.py"
/usr/bin/python3.12 "$ROOT/code/03_source_audit.py"
/usr/bin/python3.12 "$ROOT/code/03b_pseudospot.py"
if [[ "${1:-}" == '--review' ]]; then
 /usr/bin/python3.12 "$ROOT/code/04_review.py" --review "$2"
fi
/usr/bin/python3.12 "$ROOT/code/07_diagnostics.py"
/usr/bin/python3.12 "$ROOT/code/05_mp10.py"
/usr/bin/python3.12 "$ROOT/code/09_validate.py"
/usr/bin/python3.12 "$ROOT/code/08_figures.py"
/usr/bin/python3.12 "$ROOT/code/06_report.py"
