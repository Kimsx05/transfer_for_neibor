#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OMP_THREAD_LIMIT=1
exec "$ROOT/../hd_epi_niche_multimethod_v2/.venv/bin/python" "$ROOT/code/run_pipeline.py" "$@"
