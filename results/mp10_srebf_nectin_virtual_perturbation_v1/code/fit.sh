#!/usr/bin/env bash
set -euo pipefail
source code/env.sh
ulimit -v 37000000
if [[ -f "logs/time_$1.tsv" ]]; then mv "logs/time_$1.tsv" "logs/time_$1.previous.tsv"; fi
exec /usr/bin/time -v -o "logs/usage_$1.txt" Rscript --vanilla code/fit.R "$1"

