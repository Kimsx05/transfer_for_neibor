#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
PY='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/code/python.sh'
while [[ ! -f "$ROOT/02_inputs/INPUTS_READY.json" ]]; do sleep 10; done
bash "$PY" "$ROOT/code/02_reference.py" MAIN > "$ROOT/logs/02_reference_MAIN.log" 2>&1
bash "$PY" "$ROOT/code/02_reference.py" NO_P07T > "$ROOT/logs/02_reference_NO_P07T.log" 2>&1
