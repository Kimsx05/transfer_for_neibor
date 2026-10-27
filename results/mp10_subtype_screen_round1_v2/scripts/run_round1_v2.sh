#!/usr/bin/env bash
set -euo pipefail
ROOT='/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_subtype_screen_round1_v2'
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MP10_WORKERS=8
cd "$ROOT"
/usr/bin/python3 scripts/01_prepare.py
Rscript scripts/02_fill_score.R
/usr/bin/python3 scripts/03_finalize_samples.py
Rscript scripts/04_models.R
Rscript scripts/22_transform_check.R
Rscript scripts/05_cellchat.R
Rscript scripts/14_coverage.R
Rscript scripts/07_network_tables.R
Rscript scripts/08_expression_support.R
/usr/bin/python3 scripts/17_measurement_annotation.py
/usr/bin/python3 scripts/10_candidate_evidence.py
Rscript scripts/16_matched_networks.R
Rscript scripts/18_matched_results.R
Rscript scripts/15_other_groups.R
/usr/bin/python3 scripts/06_abundance_figures.py
/usr/bin/python3 scripts/09_focus_figures.py
/usr/bin/python3 scripts/11_network_figures.py
/usr/bin/python3 scripts/20_supp_figures.py
MP10_CIRCLE_MODE=primary Rscript scripts/12_circles.R
MP10_CIRCLE_MODE=secondary Rscript scripts/12_circles.R
/usr/bin/python3 scripts/21_finalize.py
