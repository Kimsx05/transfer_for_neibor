#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ -z "${DUAL_AFFINITY:-}" ]]; then
 cpus=$(/usr/bin/python3.12 -c 'import os; print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
 export DUAL_AFFINITY=1
 exec taskset -c "$cpus" bash ./run.sh
fi
source code/env.sh
timed() { local stage=$1; shift; /usr/bin/time -v -o "logs/$stage.resources" "$@" >"logs/$stage.log" 2>&1; }
for s in prepare_B supplement functional detection; do timed "$s" /usr/bin/python3.12 "code/$s.py"; done
timed map Rscript --vanilla -e 'source("code/map.R")'
for s in analyze compare figures; do timed "$s" /usr/bin/python3.12 "code/$s.py"; done
timed violins Rscript --vanilla code/violins.R
timed report /usr/bin/python3.12 code/report.py
timed verify /usr/bin/python3.12 code/verify.py
/usr/bin/python3.12 code/cleanup.py

