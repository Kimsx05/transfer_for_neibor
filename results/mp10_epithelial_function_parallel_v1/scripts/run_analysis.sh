#!/usr/bin/env bash
set -euo pipefail
src=$(cd -- "$(dirname -- "$0")/.." && pwd)
root=$(cd "$src" && pwd | sed 's@/results/mp10_epithelial_function_parallel_v1.*@@')
cd "$root"
base="$root/results/mp10_epithelial_function_parallel_v1"
dest=$(mktemp -d "$base/replay_$(date +%Y%m%d_%H%M%S)_XXXX")
cp -a "$src/scripts" "$dest/"
cp "$src/PREREGISTRATION.md" "$src/METHODS.md" "$src/FIGURE_INDEX.md" "$dest/"
mkdir -p "$dest"/{logs,figures,plot_data,tables}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
cpu=$(/usr/bin/python3 -c 'import os; print(min(os.sched_getaffinity(0)))')
stage() {
 local name=$1
 shift
 flock results/mp10_parallel_resources/heavy_object.lock taskset -c "$cpu" "$@" >"$dest/logs/$name.log" 2>&1
}
stage cache /usr/bin/python3 "$dest/scripts/cache_audit.py"
stage ctrp /usr/bin/python3 "$dest/scripts/ctrp.py"
stage depmap /usr/bin/python3 "$dest/scripts/depmap.py"
stage project Rscript "$dest/scripts/project.R" "$dest"
stage footprint_members Rscript "$dest/scripts/footprint_members.R" "$dest"
stage ferro /usr/bin/python3 "$dest/scripts/ferro_reference.py"
stage acute Rscript "$dest/scripts/acute.R" "$dest"
stage orthology Rscript "$dest/scripts/orthology.R" "$dest"
stage ctl /usr/bin/python3 "$dest/scripts/ctl.py"
stage membrane /usr/bin/python3 "$dest/scripts/membrane.py"
stage figures Rscript "$dest/scripts/figures.R" "$dest"
stage summary /usr/bin/python3 "$dest/scripts/summarize.py"
stage verification /usr/bin/python3 "$dest/scripts/verify.py"
echo "$dest"

