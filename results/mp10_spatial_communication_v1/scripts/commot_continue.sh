#!/bin/bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
ROOT='mp10_spatial_communication_v1'
while [ ! -f "$ROOT/commot/GSE171351_A1/primary2dNN/COMPLETED.json" ]; do
 if ! kill -0 1958738 2>/dev/null; then echo 'First section stopped without completion; inspect log'; exit 1; fi
 sleep 10
done
# First measured peak ~1 GiB; both continuations together capped at 24 GiB each.
(
 ulimit -v 25165824
 /usr/bin/time -v /usr/bin/python3 "$ROOT/scripts/commot_run.py" Inhousedata_T4469479 > "$ROOT/logs/commot_Inhousedata_T4469479_primary.log" 2>&1
 /usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT Inhousedata_T4469479 > "$ROOT/logs/plot_COMMOT_Inhousedata_T4469479.log" 2>&1
) & task_inhouse=$!
(
 ulimit -v 25165824
 /usr/bin/time -v /usr/bin/python3 "$ROOT/scripts/commot_run.py" GSE171351_A1 --wide > "$ROOT/logs/commot_GSE171351_A1_wide.log" 2>&1
 /usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT GSE171351_A1 --variant wide3dNN > "$ROOT/logs/plot_COMMOT_GSE171351_A1_wide.log" 2>&1
) & task_wide=$!
/usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT GSE171351_A1 > "$ROOT/logs/plot_COMMOT_GSE171351_A1.log" 2>&1
wait "$task_inhouse"
wait "$task_wide"
