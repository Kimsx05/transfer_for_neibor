#!/usr/bin/env bash
set -euo pipefail
source code/env.sh
export MP10_MONITOR_STOP="$PWD/tmp/stop.monitor.$$"
"$VP_PY" code/monitor.py &
mon=$!
trap 'touch "$MP10_MONITOR_STOP"; wait "$mon" || true' EXIT
bash code/run_pipeline.sh

