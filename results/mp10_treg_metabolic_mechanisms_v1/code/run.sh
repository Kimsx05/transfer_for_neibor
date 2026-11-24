#!/usr/bin/env bash
# Run from the original project root. Each replay gets its own dated directory.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
export TZ=Asia/Shanghai
base=results/mp10_treg_metabolic_mechanisms_v1
test -r "$base/code/project.R"
export B_OUT
B_OUT=$(mktemp -d "$base/replay_$(date +%Y%m%d)_XXXXXX")
mkdir -p "$B_OUT"/{tables,figures,logs} results/mp10_parallel_resources
cp -R "$base/code" "$B_OUT/code"
cp "$base/registration.md" "$B_OUT/"
cp "$base/tables/"{supply_evidence,branch_status}.tsv "$B_OUT/tables/"
python3 "$B_OUT/code/fetch.py" > "$B_OUT/logs/download.log" 2>&1
python3 "$B_OUT/code/prepare.py" > "$B_OUT/logs/prepare.log" 2>&1
for stage in reference local project disjoint audit figures; do
  cmd=(Rscript --vanilla "$B_OUT/code/$stage.R")
  if [[ $stage == local || $stage == audit ]]; then
    cmd=(flock results/mp10_parallel_resources/heavy_object.lock "${cmd[@]}")
  fi
  "${cmd[@]}" > "$B_OUT/logs/$stage.log" 2>&1
done
# Release the frozen interpretation only after checking that its numerical tables reproduce.
for t in model_status reference_summary recovered_sender_receiver uptake_coverage; do
  cmp "$base/tables/$t.tsv" "$B_OUT/tables/$t.tsv"
done
cp "$base/"{METHODS,RESULTS}.md "$B_OUT/"
mp10_immune_endothelial_mechanisms_v1/venv/bin/python "$B_OUT/code/validate.py" > "$B_OUT/logs/validation.log"
echo "Completed $B_OUT"
