#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 ARROW_NUM_THREADS=1 BLIS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
py="${VISIUM_PY:-$HOME/codex_singularity/python-envs/bladder-spatial-py312/bin/python}"
ulimit -v 185000000
logroot="${V3_MODEL_OUT:-$PWD}/logs"
mkdir -p "$logroot"
tag="$1"
if [[ "$1" == curves || "$1" == geometry_sensitivity ]]; then tag="$1_${2:-analyze}"; fi
exec /usr/bin/time -v -o "$logroot/time_$tag.log" taskset -c "$(awk '/Cpus_allowed_list/ {split($2,x,/[-,]/);print x[1]}' /proc/self/status)" "$py" "scripts/$1.py" "${@:2}"
