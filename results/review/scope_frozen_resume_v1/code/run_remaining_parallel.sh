#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
PY='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/code/python.sh'
while [[ ! -f "$ROOT/03_reference/NO_P07T/REFERENCE_DONE.json" ]]; do sleep 10; done
bash "$PY" "$ROOT/code/04_mapping.py" NO_P07T DEV > "$ROOT/logs/04_mapping_NO_P07T_DEV.log" 2>&1
bash "$PY" "$ROOT/code/05_freeze_DEV.py" > "$ROOT/logs/05_freeze_DEV.log" 2>&1
# Both DEV fits and threshold freeze precede both TEST fits. Each full mapping
# was measured at about7.2GiB; the two paired fits fit the23GiB GPU.
free_mib=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits)
if (( free_mib >= 17000 )); then
 bash "$PY" "$ROOT/code/04_mapping.py" MAIN TEST > "$ROOT/logs/04_mapping_MAIN_TEST.log" 2>&1 &
 main_pid=$!
 bash "$PY" "$ROOT/code/04_mapping.py" NO_P07T TEST > "$ROOT/logs/04_mapping_NO_P07T_TEST.log" 2>&1 &
 no_pid=$!
 wait "$main_pid"
 wait "$no_pid"
else
 for ref in MAIN NO_P07T; do bash "$PY" "$ROOT/code/04_mapping.py" "$ref" TEST > "$ROOT/logs/04_mapping_${ref}_TEST.log" 2>&1; done
fi
