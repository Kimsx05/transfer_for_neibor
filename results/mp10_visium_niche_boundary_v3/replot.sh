#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ $# -eq 1 ]]; then export V3_MODEL_OUT="$PWD/versions/$1"; fi
bash run.sh curves plot
bash run.sh review
bash run.sh group_plots
bash composition.sh
bash run.sh index3

