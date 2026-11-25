export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export R_LIBS="$PWD/../mp10_srebf_nectin_virtual_perturbation_v1/cache/Rlib:/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4:/refdir/Rlib"
export VP_PY="../../mp10_immune_endothelial_mechanisms_v1/venv/bin/python"
export TMPDIR="/tmp/mp10_v2_r_t070721"
if [[ ! -e "$TMPDIR" ]]; then ln -s "$PWD/tmp" "$TMPDIR"; fi
[[ "$(readlink -f "$TMPDIR")" == "$PWD/tmp" ]] || exit 3

