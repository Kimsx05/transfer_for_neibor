#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round3'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$ROOT"
Rscript scripts/01_extract.R
/usr/bin/python3 scripts/02_summarize.py
Rscript scripts/03_circles.R
/usr/bin/python3 scripts/04_figures.py
/usr/bin/python3 scripts/05_treg_figures.py
/usr/bin/python3 scripts/06_supplement.py
/usr/bin/python3 scripts/08_finalize.py
