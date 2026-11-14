#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../../.."
O="$PWD/highres_VisiumHD/analysis/hd_mp10_treg_endo_precommunication_v1"
PY="$PWD/highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python"
for step in 09_niche 10_regional_scores; do /usr/bin/time -v nice -n 10 "$PY" "$O/code/$step.py" > "$O/logs/$step.log" 2>&1; done
