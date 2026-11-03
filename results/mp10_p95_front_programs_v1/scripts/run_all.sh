#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")/../../.."
out=results/mp10_p95_front_programs_v1
/usr/bin/python3.12 "$out/scripts/prepare.py"
bash "$out/scripts/score.sh"
/usr/bin/python3.12 "$out/scripts/analyze.py" > "$out/logs/analyze.log" 2>&1
/usr/bin/python3.12 "$out/scripts/extra.py" > "$out/logs/extra.log" 2>&1
/usr/bin/python3.12 "$out/scripts/bands.py" > "$out/logs/bands.log" 2>&1
/usr/bin/python3.12 "$out/scripts/summarize.py" > "$out/logs/summarize.log" 2>&1
for i in 1 2 3 4 5 6; do Rscript "$out/scripts/detection.R" "$i" > "$out/logs/detect_$i.log" 2>&1 & done
wait
if [ ! -f "$out/tables/sc_sample_scores.tsv" ]; then
 /usr/bin/python3.12 "$out/scripts/sc_prepare.py" > "$out/logs/sc_prepare.log" 2>&1
 Rscript "$out/scripts/sc_score.R" > "$out/logs/sc_score.log" 2>&1
fi
/usr/bin/python3.12 "$out/scripts/diagnostics.py"
/usr/bin/python3.12 "$out/scripts/plots.py" > "$out/logs/plots.log" 2>&1
/usr/bin/python3.12 "$out/scripts/supplement_plots.py" > "$out/logs/supplement_plots.log" 2>&1
/usr/bin/python3.12 "$out/scripts/finalize.py"
