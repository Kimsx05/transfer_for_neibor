#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
if [[ $# != 2 || "$1" == primary || ! "$1" =~ ^[a-zA-Z0-9_-]+$ ]]; then echo 'Usage: bash regroup.sh NEW_VERSION path/to/rules.json' >&2; exit 2; fi
bash run.sh classify --version "$1" --config "$2"
for stage in analyze associate supplement influence major_lineages report_tables review; do bash run.sh "$stage" --version "$1"; done
bash replot.sh "$1"

