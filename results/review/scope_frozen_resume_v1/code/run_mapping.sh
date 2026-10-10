#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
PY='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/code/python.sh'
while [[ ! -f "$ROOT/03_reference/MAIN/REFERENCE_DONE.json" || ! -f "$ROOT/04_pseudospots/TEST/READY.json" ]]; do sleep 10; done
bash "$PY" "$ROOT/code/04_mapping.py" MAIN DEV --smoke > "$ROOT/logs/04_mapping_SMOKE.log" 2>&1
bash "$PY" "$ROOT/code/04_mapping.py" MAIN DEV > "$ROOT/logs/04_mapping_MAIN_DEV.log" 2>&1
while [[ ! -f "$ROOT/03_reference/NO_P07T/REFERENCE_DONE.json" ]]; do sleep 10; done
bash "$PY" "$ROOT/code/04_mapping.py" NO_P07T DEV > "$ROOT/logs/04_mapping_NO_P07T_DEV.log" 2>&1
bash "$PY" "$ROOT/code/05_freeze_DEV.py" > "$ROOT/logs/05_freeze_DEV.log" 2>&1
for ref in MAIN NO_P07T; do bash "$PY" "$ROOT/code/04_mapping.py" "$ref" TEST > "$ROOT/logs/04_mapping_${ref}_TEST.log" 2>&1; done
