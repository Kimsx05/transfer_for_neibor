#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
unset V3_MODEL_OUT V3_SEED V3_GRID V3_AXIS
bash run.sh metadata >logs/metadata.log 2>&1
bash run.sh curves prepare >logs/bins.log 2>&1
bash run.sh curves plot >logs/sections.log 2>&1
for stage in side_support features input_audit coverage_extra cluster review lock morphology down auxiliary; do bash run.sh "$stage" >"logs/$stage.log" 2>&1; done
bash run.sh geometry_sensitivity prepare >logs/geometry_bins.log 2>&1
for stage in geometry_sensitivity paired_curves extra_tests group_plots; do bash run.sh "$stage" >"logs/$stage.log" 2>&1; done
bash composition.sh >logs/composition.log 2>&1
bash run.sh index3 >logs/index.log 2>&1
bash run.sh validate3 >logs/validate.log 2>&1
