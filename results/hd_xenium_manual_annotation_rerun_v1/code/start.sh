#!/bin/bash
set -eu
cd "$(dirname "$0")/.."
exec systemd-run --user --unit=blca-manual-expression --collect --property=TasksMax=10 --property=MemoryMax=195000000000 --property=MemorySwapMax=0 --working-directory="$PWD" /usr/bin/python3 "$PWD/code/monitor.py" "$@"

