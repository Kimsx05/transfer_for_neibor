#!/bin/bash
set -eu
cd "$(dirname "$0")/.."
py="../../highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 POLARS_MAX_THREADS=1 RAYON_NUM_THREADS=1 TF_NUM_INTEROP_THREADS=1 TF_NUM_INTRAOP_THREADS=1 VECLIB_MAXIMUM_THREADS=1 CUDA_VISIBLE_DEVICES="" PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg NUMBA_CACHE_DIR="$PWD/logs/numba" MPLCONFIGDIR="$PWD/logs/mpl" CELLTYPIST_FOLDER="$PWD/shared/ct_cache" TF_ENABLE_ONEDNN_OPTS=0
exec systemd-run --user --slice=blca-retry27.slice --scope --quiet -p TasksMax=17 -p MemoryMax=275000000000 -p MemorySwapMax=0 taskset -c 2-5 /usr/bin/time -a -o logs/resources.tsv -f '%C\t%e\t%M\t%x' "$py" -u "$@"
