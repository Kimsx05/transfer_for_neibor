#!/usr/bin/env bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd '/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3'
if [[ ! -f 01_partitions/COMPLETE ]]; then echo 'Partition scan not complete' >&2; exit 1; fi
if [[ ! -f 02_coverage/COMPLETE ]]; then /usr/local/bin/Rscript code/02_coverage_candidates.R > logs/02_coverage_candidates.log 2>&1; fi
/usr/local/bin/Rscript code/03_aggregate.R > logs/03_aggregate.log 2>&1
/usr/local/bin/Rscript code/04_models.R > logs/04_models.log 2>&1
if [[ ! -f 08_review/TABLES_COMPLETE ]]; then /usr/local/bin/Rscript code/05_review.R > logs/05_review.log 2>&1; fi
