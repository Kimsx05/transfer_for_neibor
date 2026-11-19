#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
py=../hd_epi_niche_multimethod_v2/.venv/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [[ ! -f D1_Johnson/historical_candidates.parquet || ! -f Cervilla/historical_candidates.parquet ]]; then "$py" code/screen.py; fi
"$py" code/reconstruct.py
"$py" code/domain_run.py
if [[ ! -f D1_Johnson/zoom_index.parquet || ! -f Cervilla/zoom_index.parquet ]]; then "$py" code/pictures.py details; fi
"$py" code/score_stats.py
"$py" code/pictures.py second
"$py" code/focus.py
"$py" code/deliver.py
"$py" code/review_sheets.py
"$py" code/candidate_neighbors.py
"$py" code/finish.py
