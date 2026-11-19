#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
py=../hd_epi_niche_multimethod_v2/.venv/bin/python
"$py" code/audit.py
[[ -s Cervilla/count_accounting.tsv && -s D1_Johnson/count_accounting.tsv ]] || "$py" code/units.py
"$py" code/source.py
"$py" code/localization.py
"$py" code/historical.py
[[ -s D1_Johnson/existing_cell8_evidence.parquet && -s reference_validation.tsv ]] || "$py" code/worker.py
"$py" code/summary.py
"$py" code/atlas.py
"$py" code/finish.py
"$py" code/package_result.py
