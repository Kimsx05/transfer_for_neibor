#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
samples=(Cervilla A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2)
stage(){ bash code/run.sh "$@"; }
stage code/audit.py
stage code/review.py
stage code/score.py
stage code/graph.py
stage code/windows.py
stage code/evidence.py
stage code/niche.py mender "${samples[@]}"
stage bash code/actual_R.sh code/build.R "${samples[@]}"
stage code/niche.py cellcharter "${samples[@]}"
stage code/spatial.py
stage code/reverse.py
stage code/squid.py "${samples[@]}"
stage code/export_tools.py
stage bash code/actual_R.sh code/spatstat.R "${samples[@]}"
stage bash code/actual_R.sh code/crawdad.R "${samples[@]}"
stage g++ -shared -fPIC -O2 -I/usr/local/lib/R/include -o code/crawdad_counts.so code/crawdad_counts.cpp
stage bash code/actual_R.sh code/formal_crawdad.R "${samples[@]}"
stage code/target_tools.py
stage code/kprep.py
stage bash code/actual_R.sh code/kcross.R "${samples[@]}"
stage code/niche_results.py
stage code/check.py
stage code/pack.py
stage code/more_figures.py
stage bash code/actual_R.sh code/native_primary.R "${samples[@]}"
stage code/figures.py
stage bash code/actual_R.sh code/violin.R "${samples[@]}"
stage code/extras.py
stage bash code/actual_R.sh code/versions.R
stage code/final_audit.py
stage code/captions.py
stage code/report.py
stage code/compare.py
stage code/report.py
stage code/interpretation.py
