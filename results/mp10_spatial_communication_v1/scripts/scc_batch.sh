#!/usr/bin/env bash
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -v 50331648
for section in GSE171351_A1 Inhousedata_T4469479; do
  if [[ -f "mp10_spatial_communication_v1/spatialcellchat/${section}/main/COMPLETED.txt" ]]; then
    printf 'Reuse complete main: %s\n' "$section" >> mp10_spatial_communication_v1/logs/scc_cache_reuse.log
    continue
  fi
  /usr/bin/time -v /usr/local/bin/Rscript mp10_spatial_communication_v1/scripts/scc_run.R "$section" 2 > "mp10_spatial_communication_v1/logs/scc_${section}_main.log" 2>&1
  status=$?
  printf '%s\t%s\n' "$section" "$status" >> mp10_spatial_communication_v1/logs/scc_batch_exit.tsv
done
if [[ -f mp10_spatial_communication_v1/spatialcellchat/GSE171351_A1/range3dNN/COMPLETED.txt ]]; then
  printf 'Reuse complete first-slice 3dNN sensitivity\n' >> mp10_spatial_communication_v1/logs/scc_cache_reuse.log
else
  /usr/bin/time -v /usr/local/bin/Rscript mp10_spatial_communication_v1/scripts/scc_run.R GSE171351_A1 3 > mp10_spatial_communication_v1/logs/scc_GSE171351_A1_range3dNN.log 2>&1
fi
for destination in GSE171351_A1/main Inhousedata_T4469479/main GSE171351_A1/range3dNN; do
  if [[ ! -f "mp10_spatial_communication_v1/spatialcellchat/${destination}/COMPLETED.txt" ]]; then
    printf 'Incomplete SpatialCellChat run: %s\n' "$destination" >&2
    exit 1
  fi
done
/usr/local/bin/Rscript mp10_spatial_communication_v1/scripts/scc_distance_QA.R > mp10_spatial_communication_v1/logs/scc_distance_QA.log 2>&1
/usr/local/bin/Rscript mp10_spatial_communication_v1/scripts/scc_sensitivity_QA.R > mp10_spatial_communication_v1/logs/scc_sensitivity_QA.log 2>&1
/usr/bin/python3 mp10_spatial_communication_v1/scripts/scc_finalize.py
/usr/bin/python3 mp10_spatial_communication_v1/scripts/scc_summary_notes.py
