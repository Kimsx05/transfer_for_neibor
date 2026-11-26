#!/bin/bash
set -u
cd "$(dirname "$0")/.."
samples=(Cervilla A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2)
job(){ bash code/run.sh "$@" || echo "failed $*"; }
job bash code/actual_R.sh code/crawdad_debug.R
job bash code/actual_R.sh code/crawdad.R A1R5 A1R6 A2R2
job bash code/actual_R.sh code/formal_crawdad.R "${samples[@]}"
job bash code/actual_R.sh code/kcross.R "${samples[@]}"
job bash code/actual_R.sh code/spatstat.R "${samples[@]}"
job code/extras.py
job code/pack.py
job code/more_figures.py
job bash code/actual_R.sh code/versions.R
job bash code/actual_R.sh code/native_primary.R "${samples[@]}"
job code/final_audit.py
job code/report.py
job code/compare.py
job code/report.py
