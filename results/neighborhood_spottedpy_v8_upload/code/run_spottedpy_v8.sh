#!/usr/bin/env bash
set -euo pipefail
v8_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
exec "$v8_dir/env/bin/python" "$v8_dir/code/run.py" "$@"
