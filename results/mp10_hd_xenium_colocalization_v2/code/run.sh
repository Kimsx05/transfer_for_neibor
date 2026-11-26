#!/bin/bash
set -eu
cd "$(dirname "$0")/.."
exec systemd-run --user --scope --quiet -p TasksMax=7 -p MemoryHigh=155000000000 -p MemoryMax=165000000000 -p MemorySwapMax=0 python3 code/monitor.py "$@"
