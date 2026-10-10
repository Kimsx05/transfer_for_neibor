#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1'
PY='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_k32_reference_identifiability_v5/code/python.sh'
run(){ bash "$PY" "$ROOT/code/$1" "${@:2}"; }
run 01_prepare.py
run 03_pseudospots.py
run 00_validate_inputs.py
run 09_evidence_tables.py
run 02_reference.py MAIN
run 02_reference.py NO_P07T
run 04_mapping.py MAIN DEV --smoke
run 04_mapping.py MAIN DEV
run 04_mapping.py NO_P07T DEV
run 05_freeze_DEV.py
run 04_mapping.py MAIN TEST
run 04_mapping.py NO_P07T TEST
run 06_evaluate.py
run 07_figures.py
if [[ -f "$ROOT/08_review/ASSESSMENT.json" ]]; then run 08_report_package.py; fi
