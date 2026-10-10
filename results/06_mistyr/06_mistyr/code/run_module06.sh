#!/usr/bin/env bash
set -euo pipefail
module_root='/home/data/t070721/codex_workspace/Bladder Metabolism/spatial_ecology_top20_v2/06_mistyr'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=''
cd "$module_root"
/usr/bin/python3.12 -S code/prepare.py > logs/prepare.log 2>&1
/usr/bin/python3.12 -S code/run_all.py > logs/run_all.log 2>&1
Rscript code/correct_importance.R "$module_root" > logs/correct_importance.log 2>&1
/usr/bin/python3.12 -S code/summarize.py > logs/summarize.log 2>&1
Rscript code/plots.R "$module_root" > logs/plots.log 2>&1
/usr/bin/python3.12 -S code/finalize.py > logs/finalize.log 2>&1
/usr/bin/python3.12 -S code/seal_manifest.py
