#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Respect the existing allocation and keep native libraries within the four-CPU budget.
if [[ -z "${PILOT_AFFINITY_SET:-}" ]]; then
  cpus=$(/usr/bin/python3.12 -c 'import os;print(",".join(map(str,sorted(os.sched_getaffinity(0))[:4])))')
  export PILOT_AFFINITY_SET=1
  exec taskset -c "$cpus" bash ./run.sh
fi
source code/env.sh
for stage in prepare audit; do
 /usr/bin/time -v -o "logs/$stage.resources" /usr/bin/python3.12 "code/$stage.py" >"logs/$stage.log" 2>&1
done
bash code/map.sh
for stage in analyze figures report verify; do
 /usr/bin/time -v -o "logs/$stage.resources" /usr/bin/python3.12 "code/$stage.py" >"logs/$stage.log" 2>&1
done
/usr/bin/python3.12 code/clean_scratch.py

