#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=results/mp10_p95_front_programs_v1
while ! rg -q 'summary complete' "$out/logs/summarize.log"; do sleep 5; done
/usr/bin/python3.12 "$out/scripts/diagnostics.py" > "$out/logs/diagnostics.log" 2>&1
/usr/bin/python3.12 "$out/scripts/plots.py" > "$out/logs/plots.log" 2>&1
