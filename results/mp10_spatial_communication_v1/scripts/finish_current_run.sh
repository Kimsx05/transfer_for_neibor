#!/bin/bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
ROOT='mp10_spatial_communication_v1'
while kill -0 1993775 2>/dev/null; do sleep 10; done
if [ ! -f "$ROOT/commot/Inhousedata_T4469479/primary2dNN/COMPLETED.json" ]; then echo 'Inhouse not complete; inspect inference log'; exit 1; fi
/usr/bin/python3 "$ROOT/scripts/spatial_postprocess.py" COMMOT Inhousedata_T4469479 > "$ROOT/logs/plot_COMMOT_Inhousedata_T4469479.log" 2>&1
/usr/bin/python3 "$ROOT/scripts/commot_distance_audit.py" > "$ROOT/logs/commot_distance_audit.log" 2>&1
/usr/bin/python3 "$ROOT/scripts/spatial_finalize.py" > "$ROOT/logs/spatial_finalize.log" 2>&1
