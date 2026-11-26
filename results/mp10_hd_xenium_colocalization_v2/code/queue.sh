#!/bin/bash
set -u
cd "$(dirname "$0")/.."
while ! grep -q 'squid.py_.*[[:space:]]0[[:space:]]' audit/resources.tsv && ! grep -q 'squid.py_.*[[:space:]]1[[:space:]]' audit/resources.tsv; do sleep 3; done
bash code/run.sh bash code/actual_R.sh code/install.R
bash code/run.sh code/export_tools.py
bash code/run.sh bash code/actual_R.sh code/build.R Cervilla A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2
bash code/run.sh code/windows.py
bash code/run.sh code/niche.py cellcharter A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2
bash code/run.sh code/reverse.py
bash code/run.sh bash code/actual_R.sh code/spatstat.R Cervilla A1R2 D1_Johnson UC_SCC A1R4 A1R5 A1R6 A2R2
