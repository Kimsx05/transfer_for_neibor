#!/usr/bin/env bash
export OPENBLAS_NUM_THREADS=${VP_BLAS_THREADS:-3} OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
export R_LIBS="$PWD/cache/Rlib:/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4:/refdir/Rlib"
# R rejects whitespace in its temporary directory. Temporary data stay inside the result directory.
link="/tmp/mp10_vp_t070721_rtmp"
if [[ ! -e "$link" ]]; then ln -s "$PWD/tmp" "$link"; fi
[[ "$(readlink -f "$link")" == "$PWD/tmp" ]] || exit 3
export TMPDIR="$link" PYTHONDONTWRITEBYTECODE=1
export VP_PY="../../mp10_immune_endothelial_mechanisms_v1/venv/bin/python"

