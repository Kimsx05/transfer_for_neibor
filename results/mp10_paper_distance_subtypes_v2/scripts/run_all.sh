#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$(dirname "$(dirname "$(realpath "$0")")")")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
p=mp10_paper_distance_subtypes_v2
python3.12 -B "$p/scripts/run.py"
python3.12 -B "$p/scripts/niche_extra.py"
python3.12 -B "$p/scripts/ranks.py"
Rscript "$p/scripts/enrich.R"
python3.12 -B "$p/scripts/sc_prepare.py"
Rscript "$p/scripts/sc_score.R"
python3.12 -B "$p/scripts/summarize.py"
python3.12 -B "$p/scripts/plots.py"
python3.12 -B "$p/scripts/support.py"
python3.12 -B "$p/scripts/comparisons.py"
python3.12 -B "$p/scripts/axis_composition.py"
python3.12 -B "$p/scripts/finalize.py"
