#!/bin/bash
set -euo pipefail
base="$(cd "$(dirname "$0")/.." && pwd)"
plan="$(mktemp "$base/logs/review_plan_XXXXXX")"
export MP_V2_PLAN_FILE="$plan"
bash "$base/code/run.sh" code/import_review.py "$@"
if [[ -s "$plan" ]]; then
 IFS= read -r entry < "$plan"
 exec /bin/bash "$entry"
fi
