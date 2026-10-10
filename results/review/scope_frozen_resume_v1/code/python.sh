#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism'
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$ROOT/cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
exec /usr/bin/python3.12 -S "$@"
