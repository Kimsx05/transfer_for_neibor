#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
out=results/mp10_section_specificity_pearson_v1
while [ "$(find "$out/tables" -maxdepth 1 -name '*_status.tsv' | wc -l)" -lt 63 ] || [ "$(find "$out/tables" -maxdepth 1 -name '*_context_done.tsv' | wc -l)" -lt 63 ]; do sleep 10; done
/usr/bin/python3.12 "$out/scripts/aggregate.py"
/usr/bin/python3.12 "$out/scripts/context_summary.py"
while [ "$(find "$out/reports" -maxdepth 1 -name '*_figure_index.tsv' | wc -l)" -lt 63 ]; do sleep 10; done
/usr/bin/python3.12 "$out/scripts/overview.py"
