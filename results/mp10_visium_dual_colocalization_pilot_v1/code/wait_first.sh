set -eu
source code/env.sh
# First run may overlap CPU-isolated preparation; normal run.sh is sequential.
while [[ ! -f tables/gene_coverage_base.tsv ]]; do sleep 3; done
/usr/bin/time -v -o logs/map.resources taskset -c 2-3 Rscript --vanilla code/map.R >logs/map.log 2>&1

