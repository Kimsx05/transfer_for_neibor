#!/usr/bin/env bash
set -euo pipefail
scripts=$(dirname "$(realpath "$0")")
export MP10_OUTPUT=${MP10_OUTPUT:-$(dirname "$scripts")}
export MP10_PROJECT=${MP10_PROJECT:-$(dirname "$MP10_OUTPUT")}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
cd "$MP10_PROJECT"
python3.12 "$scripts/register.py"
bash "$scripts/score.sh"
python3.12 "$scripts/distance.py"
python3.12 "$scripts/context.py"
python3.12 "$scripts/extended.py"
python3.12 "$scripts/summarize.py"
python3.12 "$scripts/render.py" core --workers 8
python3.12 "$scripts/render.py" program --workers 8
python3.12 "$scripts/overview.py"
python3.12 "$scripts/package.py"
