#!/bin/bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONPATH="$root/.packages${PYTHONPATH:+:$PYTHONPATH}"
exec "$root/../cell_annotation_v2/code/run_gpu.sh" "$@"
