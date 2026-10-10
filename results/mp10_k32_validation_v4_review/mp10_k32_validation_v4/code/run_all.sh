#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_validation_v4'
cd "$ROOT"
python3 code/00_check_inputs.py
/usr/local/bin/Rscript code/00_common.R
if [[ ! -f 00_audit/COMPLETE ]]; then /usr/local/bin/Rscript code/01_audit_scores.R > logs/01_audit_scores.log 2>&1; fi
if [[ ! -f 02_rna_DEG/C14_AGGREGATE_COMPLETE ]]; then /usr/local/bin/Rscript code/02_C14_aggregate.R > logs/02_C14_aggregate.log 2>&1; fi
if [[ ! -f 02_rna_DEG/COMPLETE ]]; then /usr/local/bin/Rscript code/03_models.R > logs/03_models.log 2>&1; fi
if [[ ! -f 06_supplement/SUMMARIES_COMPLETE ]]; then /usr/local/bin/Rscript code/04_summaries.R > logs/04_summaries.log 2>&1; fi
if [[ ! -f 06_supplement/FIGURES_COMPLETE ]]; then /usr/local/bin/Rscript code/05_figures.R > logs/05_figures.log 2>&1; fi
python3 code/06_audit_manifest.py > logs/06_audit_manifest.log 2>&1
python3 code/07_report.py > logs/07_report.log 2>&1
/usr/local/bin/Rscript code/08_validate.R > logs/08_validate.log 2>&1
python3 code/09_package.py > logs/09_package.log 2>&1
