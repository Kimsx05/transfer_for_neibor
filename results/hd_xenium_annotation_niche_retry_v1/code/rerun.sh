#!/bin/bash
set -eu
cd "$(dirname "$0")/.."
systemctl --user start blca-retry27.slice
systemctl --user set-property --runtime blca-retry27.slice MemoryMax=285000000000 TasksMax=20 CPUQuota=1800%
stage(){ bash code/run.sh "$@"; }
stage code/audit.py
stage code/reconstruct.py pilot
test -f shared/reconstruction_freeze.md
stage code/reconstruct.py full
stage code/image_support.py
stage code/geometry_qc.py
stage code/primary_check.py
stage code/components.py
for platform in xenium hd; do
 stage code/embed.py "$platform"
 stage code/classify.py "$platform"
 stage code/annotate.py "$platform"
 stage code/niche.py "$platform"
 stage code/deliver.py "$platform"
 stage code/qc_summary.py "$platform"
done
stage code/roi_evidence.py
stage code/old_gates.py
stage code/validate.py
stage code/lineage_markers.py
stage code/index_figures.py
stage code/status.py
stage code/report.py

