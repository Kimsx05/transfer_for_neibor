#!/bin/bash
set -eu
export HD6_BRANCH=native OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
py=../hd_epi_niche_multimethod_v2/.venv/bin/python
s="$1"
test -s "$s/native_full/out/cells.geojson.gz"
"$py" code/proseg_observed.py "$s"
"$py" code/identity.py "$s"
"$py" code/cards.py "$s" index
