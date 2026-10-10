#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd '/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3'
[[ -f 04_pseudobulk_markers/COMPLETE ]]
if [[ ! -f 08_review/TABLES_COMPLETE ]]; then /usr/local/bin/Rscript code/05_review.R > logs/05_review.log 2>&1; fi
if [[ ! -f 06_stability/candidate_Jaccard_and_profile_concordance.tsv ]]; then /usr/local/bin/Rscript code/05c_expression_stability.R > logs/05c_expression_stability.log 2>&1; fi
if [[ ! -f 07_figures/COMPLETE ]]; then /usr/local/bin/Rscript code/06_figures.R > logs/06_figures.log 2>&1; fi
python3 code/07_report.py > logs/07_report.log 2>&1
/usr/local/bin/Rscript code/08_validate.R > logs/08_validate.log 2>&1
python3 code/09_package.py > logs/09_package.log 2>&1
