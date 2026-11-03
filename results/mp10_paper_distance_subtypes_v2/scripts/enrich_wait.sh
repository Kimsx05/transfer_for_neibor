#!/usr/bin/env bash
set -e
while [[ $(find mp10_paper_distance_subtypes_v2/tables -name '*_gene_ranks.tsv.gz' | wc -l) -lt 61 ]]; do sleep 10; done
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript mp10_paper_distance_subtypes_v2/scripts/enrich.R > mp10_paper_distance_subtypes_v2/logs/enrich.log 2>&1
