#!/usr/bin/env bash
set -euo pipefail
v2_dir="$(cd "$(dirname "$0")/.." && pwd)"
base_py="$v2_dir/../cell_annotation_freeze_v1/.venv/bin/python"
export PYTHONPATH="$v2_dir/.gpu_packages${PYTHONPATH:+:$PYTHONPATH}"
extra_libs="$("$base_py" -c 'from pathlib import Path; import sys; print(":".join(str(p) for p in Path(sys.argv[1]).glob("nvidia/*/lib")))' "$v2_dir/.gpu_packages")"
export LD_LIBRARY_PATH="$extra_libs${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$v2_dir/.gpu_packages/nvidia/cuda_nvcc/bin:$PATH"
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8
exec "$base_py" "$@"
