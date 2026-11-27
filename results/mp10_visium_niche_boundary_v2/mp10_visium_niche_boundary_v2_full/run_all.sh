#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
bash run.sh preflight
for stage in base metrics blind_maps discover classify measure audit mask_audit analyze associate supplement influence major_lineages report_tables review; do bash run.sh "$stage" >"logs/$stage.log" 2>&1; done
bash replot.sh

bash run.sh validate >logs/validate.log 2>&1
python3 scripts/package.py
