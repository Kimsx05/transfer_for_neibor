#!/bin/bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
ROOT='mp10_spatial_communication_v1'
/usr/bin/python3 "$ROOT/scripts/prepare_inputs.py"
for section in GSE171351_A1 Inhousedata_T4469479; do
 if [ "$section" = Inhousedata_T4469479 ]; then
  /usr/bin/python3 "$ROOT/scripts/commot_run.py" "$section" --parallel-components
 else
  /usr/bin/python3 "$ROOT/scripts/commot_run.py" "$section"
 fi
 /usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT "$section"
done
/usr/bin/python3 "$ROOT/scripts/commot_run.py" GSE171351_A1 --wide
/usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT GSE171351_A1 --variant wide3dNN
bash "$ROOT/scripts/scc_batch.sh"
for section in GSE171351_A1 Inhousedata_T4469479; do
 /usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" SpatialCellChat "$section"
done
/usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" SpatialCellChat GSE171351_A1 --variant range3dNN
/usr/bin/python3 "$ROOT/scripts/spatial_finalize.py"
