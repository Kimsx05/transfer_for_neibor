#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
AUDIT_PYTHON=/home/data/t070721/codex_singularity/python-envs/bladder-spatial-py312/bin/python
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2
"$AUDIT_PYTHON" -u code/download_inputs.py >> logs/download_inputs.log 2>&1
"$AUDIT_PYTHON" -u code/audit_inputs.py >> logs/audit_inputs.log 2>&1
"$AUDIT_PYTHON" -u code/finalize_audit.py >> logs/finalize.log 2>&1
