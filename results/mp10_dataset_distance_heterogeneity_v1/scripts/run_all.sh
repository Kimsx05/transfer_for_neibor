#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
out=results/mp10_dataset_distance_heterogeneity_v1
force=${FORCE:-0}
if [ "$force" = 1 ] || [ ! -s "$out/tables/unit_background.tsv.gz" ]; then
 /usr/bin/python3.12 "$out/scripts/prepare.py" > "$out/logs/prepare.log" 2>&1
 /usr/bin/python3.12 "$out/scripts/state_supplement.py" > "$out/logs/state.log" 2>&1
fi
if [ "$force" = 1 ] || [ ! -s "$out/tables/search_scope.tsv" ]; then
 /usr/bin/python3.12 "$out/scripts/search.py" > "$out/logs/search.log" 2>&1
fi
if [ "$force" = 1 ] || [ ! -s "$out/tables/stability_summary.tsv" ]; then
 /usr/bin/python3.12 "$out/scripts/heterogeneity.py" > "$out/logs/heterogeneity.log" 2>&1
 /usr/bin/python3.12 "$out/scripts/representative_sensitivity.py" > "$out/logs/representative_sensitivity.log" 2>&1
fi
if [ "$force" = 1 ] || [ ! -s "$out/tables/ecology_sensitivity.tsv.gz" ]; then
 /usr/bin/python3.12 "$out/scripts/ecology.py" > "$out/logs/ecology.log" 2>&1
fi
/usr/bin/python3.12 "$out/scripts/select_results.py" > "$out/logs/selection.log" 2>&1
/usr/bin/python3.12 "$out/scripts/plots.py" > "$out/logs/plots.log" 2>&1
/usr/bin/python3.12 "$out/scripts/finalize.py" > "$out/logs/finalize.log" 2>&1
