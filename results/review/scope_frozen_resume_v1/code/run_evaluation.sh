#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
PY='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/code/python.sh'
while [[ ! -f "$ROOT/05_deconvolution/MAIN_TEST/MAPPING_DONE.json" || ! -f "$ROOT/05_deconvolution/NO_P07T_TEST/MAPPING_DONE.json" ]]; do sleep 10; done
bash "$PY" "$ROOT/code/06_evaluate.py" > "$ROOT/logs/06_evaluate.log" 2>&1
bash "$PY" "$ROOT/code/07_figures.py" > "$ROOT/logs/07_figures.log" 2>&1
