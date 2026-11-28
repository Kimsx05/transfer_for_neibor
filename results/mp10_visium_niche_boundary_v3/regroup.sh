#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ $# -lt 1 || ! "$1" =~ ^[a-zA-Z0-9_-]+$ ]]; then echo 'Usage: bash regroup.sh NEW_NAME [grid] [T|I] [seed]' >&2; exit 2; fi
export V3_MODEL_OUT="$PWD/versions/$1" V3_GRID="${2:-500}" V3_AXIS="${3:-T}" V3_SEED="${4:-20261127}"
if [[ -e "$V3_MODEL_OUT/classification_lock.json" ]]; then echo 'Output version already locked. Use a new version name.' >&2; exit 2; fi
mkdir -p "$V3_MODEL_OUT/logs"
for stage in cluster review lock down auxiliary geometry_sensitivity paired_curves extra_tests group_plots; do bash run.sh "$stage" >"$V3_MODEL_OUT/logs/$stage.log" 2>&1; done
bash composition.sh >"$V3_MODEL_OUT/logs/composition.log" 2>&1
bash run.sh index3
