#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd '/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_followup_v2'
if [[ ! -f 00_audit/COMPLETE ]]; then /usr/local/bin/Rscript code/01_audit_usage.R > logs/01_audit_usage.log 2>&1; fi
/usr/local/bin/Rscript code/02_aggregate.R > logs/02_aggregate.log 2>&1
/usr/local/bin/Rscript code/03_pseudobulk.R > logs/03_pseudobulk.log 2>&1
if [[ ! -f 03_marker_review/review_tables_COMPLETE ]]; then /usr/local/bin/Rscript code/04_review.R > logs/04_review.log 2>&1; fi
/usr/local/bin/Rscript code/06_decisions.R > logs/06_decisions.log 2>&1
if [[ ! -f 04_figures_by_dataset/COMPLETE ]]; then /usr/local/bin/Rscript code/05_figures.R > logs/05_figures.log 2>&1; fi
python3 code/07_report.py > logs/07_report.log 2>&1
/usr/local/bin/Rscript code/08_validate.R > logs/08_validate.log 2>&1
python3 code/09_package.py > logs/09_package.log 2>&1
