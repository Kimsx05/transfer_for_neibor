#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=results/mp10_section_specificity_pearson_v1
/usr/bin/python3.12 "$out/scripts/prepare.py"
bash "$out/scripts/score.sh" > "$out/logs/scoring_driver.log" 2>&1 &
score_pid=$!
/usr/bin/python3.12 "$out/scripts/dispatch.py" > "$out/logs/dispatch.log" 2>&1 &
stat_pid=$!
/usr/bin/python3.12 "$out/scripts/context_dispatch.py" > "$out/logs/context.log" 2>&1 &
context_pid=$!
/usr/bin/python3.12 "$out/scripts/report_dispatch.py" > "$out/logs/reports.log" 2>&1 &
report_pid=$!
if [ ! -f "$out/tables/sc_sample_program_source.tsv.gz" ]; then
 /usr/bin/python3.12 "$out/scripts/sc_prepare.py" > "$out/logs/sc.log" 2>&1
 Rscript "$out/scripts/sc_score.R" >> "$out/logs/sc.log" 2>&1
fi
wait "$score_pid" "$stat_pid" "$context_pid" "$report_pid"
/usr/bin/python3.12 "$out/scripts/aggregate.py" > "$out/logs/aggregate.log" 2>&1
/usr/bin/python3.12 "$out/scripts/context_summary.py" >> "$out/logs/aggregate.log" 2>&1
/usr/bin/python3.12 "$out/scripts/overview.py" > "$out/logs/overview.log" 2>&1
/usr/bin/python3.12 "$out/scripts/discussion.py" > "$out/logs/discussion.log" 2>&1
/usr/bin/python3.12 "$out/scripts/finalize.py" > "$out/logs/finalize.log" 2>&1
