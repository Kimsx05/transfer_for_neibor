#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
# flock in scheduler.py prevents two schedulers from owning the same queue.
# If a stopped scheduler left live inference children, inspect those PIDs first;
# the scheduler intentionally refuses to duplicate live workers.
exec /usr/bin/python3 mp10_spatial_communication_v2_all_sections/scripts/scheduler.py "$@"
