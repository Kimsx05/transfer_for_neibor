#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/usr/bin/python3.12 neighborhood_ring_harmony_v3/code/run_niche_only.py
if [[ -f neighborhood_ring_harmony_v3/code/finish.py ]]; then /usr/bin/python3.12 neighborhood_ring_harmony_v3/code/finish.py; fi

if [[ -f neighborhood_ring_harmony_v3/reports/review.zip ]]; then /usr/bin/python3.12 neighborhood_ring_harmony_v3/code/final_integrity.py; fi
