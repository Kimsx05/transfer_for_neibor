set -eu
source code/env.sh
for stage in analyze compare; do /usr/bin/time -v -o "logs/$stage.resources" taskset -c 0-1 /usr/bin/python3.12 "code/$stage.py" >"logs/$stage.log" 2>&1; done
/usr/bin/time -v -o logs/violins.resources taskset -c 2-3 Rscript --vanilla code/violins.R >logs/violins.log 2>&1 &
vpid=$!
/usr/bin/time -v -o logs/figures.resources taskset -c 0-1 /usr/bin/python3.12 code/figures.py >logs/figures.log 2>&1
wait "$vpid"

