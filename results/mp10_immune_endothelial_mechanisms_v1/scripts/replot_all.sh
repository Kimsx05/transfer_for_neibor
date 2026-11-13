#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_THREAD_LIMIT=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
Rscript --vanilla mp10_immune_endothelial_mechanisms_v1/scripts/11_replot.R
Rscript --vanilla mp10_immune_endothelial_mechanisms_v1/scripts/22_replot_extended.R
Rscript --vanilla mp10_immune_endothelial_mechanisms_v1/scripts/33_replot_integration.R
mp10_immune_endothelial_mechanisms_v1/venv/bin/python mp10_immune_endothelial_mechanisms_v1/scripts/26_render_audit.py
