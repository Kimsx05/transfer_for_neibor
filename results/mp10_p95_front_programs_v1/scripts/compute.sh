#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=results/mp10_p95_front_programs_v1
while [ "$(find "$out/scores" -name '*_programs.tsv.gz' | wc -l)" -lt 63 ]; do sleep 15; done
/usr/bin/python3.12 "$out/scripts/analyze.py" > "$out/logs/analyze.log" 2>&1
/usr/bin/python3.12 "$out/scripts/extra.py" > "$out/logs/extra.log" 2>&1
/usr/bin/python3.12 "$out/scripts/summarize.py" > "$out/logs/summarize.log" 2>&1
