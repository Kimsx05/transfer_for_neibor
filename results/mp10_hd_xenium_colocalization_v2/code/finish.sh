#!/bin/bash
set -u
cd "$(dirname "$0")/.."
samples=(Cervilla A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2)
job(){ bash code/run.sh "$@" || echo "failed $*"; }
job code/check.py
job code/evidence.py
job code/niche.py cellcharter A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2
job code/reverse.py
job bash code/actual_R.sh code/crawdad.R "${samples[@]}"
job bash code/actual_R.sh code/formal_crawdad.R "${samples[@]}"
job code/target_tools.py
job code/niche_results.py
job code/figures.py
job bash code/actual_R.sh code/violin.R "${samples[@]}"
job code/kprep.py
job bash code/actual_R.sh code/kcross.R "${samples[@]}"
job bash code/actual_R.sh code/spatstat.R "${samples[@]}"
job code/report.py
job code/compare.py
job code/report.py
job code/import_review.py review/Cervilla/human_review.tsv
