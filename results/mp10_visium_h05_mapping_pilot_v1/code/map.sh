set -eu
source code/env.sh
/usr/bin/time -v -o logs/map.resources Rscript --vanilla code/map.R >logs/map.log 2>&1

