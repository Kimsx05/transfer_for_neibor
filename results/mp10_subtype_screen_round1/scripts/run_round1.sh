#!/usr/bin/env bash
set -euo pipefail
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset MP10_DATASETS_ONLY MP10_ALL_ONLY MP10_SKIP_RUNNING_NORMAL
TASK_RUN='mp10_subtype_screen_round1'
/usr/bin/python3 "$TASK_RUN/scripts/01_samples.py" > "$TASK_RUN/logs/01_samples.log" 2>&1
if [[ ! -f "$TASK_RUN/objects/signaling_RNA_data.qs" ]]; then
 Rscript "$TASK_RUN/scripts/02_extract_expression.R" > "$TASK_RUN/logs/02_extract_expression.log" 2>&1
fi
Rscript "$TASK_RUN/scripts/validate_fast_trimean.R" > "$TASK_RUN/logs/validate_fast_trimean.log" 2>&1
Rscript "$TASK_RUN/scripts/03_proportion_models.R" > "$TASK_RUN/logs/03_proportion_models.log" 2>&1
Rscript "$TASK_RUN/scripts/04_cellchat.R" > "$TASK_RUN/logs/04_cellchat.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/05_abundance_figures.py" > "$TASK_RUN/logs/05_abundance_figures.log" 2>&1
Rscript "$TASK_RUN/scripts/06_model_sensitivity.R" > "$TASK_RUN/logs/06_model_sensitivity.log" 2>&1
Rscript "$TASK_RUN/scripts/07_network_tables.R" > "$TASK_RUN/logs/07_network_tables.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/08_network_overviews.py" > "$TASK_RUN/logs/08_network_overviews.log" 2>&1
Rscript "$TASK_RUN/scripts/09_candidate_support.R" > "$TASK_RUN/logs/09_candidate_support.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/10_candidate_evidence.py" > "$TASK_RUN/logs/10_candidate_evidence.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/select_review_axes.py" > "$TASK_RUN/logs/selected_axes_details.log" 2>&1
Rscript "$TASK_RUN/scripts/11_circle_figures.R" > "$TASK_RUN/logs/11_circle_figures.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/12_candidate_figures.py" > "$TASK_RUN/logs/12_candidate_figures.log" 2>&1
Rscript "$TASK_RUN/scripts/13_selected_gene_support.R" > "$TASK_RUN/logs/13_selected_gene_support.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/13b_candidate_network_context.py" > "$TASK_RUN/logs/13b_candidate_network_context.log" 2>&1
/usr/bin/python3 "$TASK_RUN/scripts/14_finalize.py" > "$TASK_RUN/logs/14_finalize.log" 2>&1
