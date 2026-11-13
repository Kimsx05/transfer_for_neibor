#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_THREAD_LIMIT=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
mp10_mode="${1:-resume}"
mp10_scripts='mp10_immune_endothelial_mechanisms_v1/scripts'
mp10_results='results/mp10_immune_endothelial_mechanisms_v1'
mp10_python='mp10_immune_endothelial_mechanisms_v1/venv/bin/python'
# Serial modules; AUCell/key-seed use at most two single-thread workers. Keep total task workers <=4.
if [ "$mp10_mode" = rebuild_pseudobulk ]; then
  Rscript --vanilla "$mp10_scripts/00_inspect_object.R"
  Rscript --vanilla "$mp10_scripts/01_pseudobulk.R"
  Rscript --vanilla "$mp10_scripts/01b_finalize_audit.R"
elif [ "$mp10_mode" = resume ]; then
  "$mp10_python" "$mp10_scripts/25_verify_inputs.py"
else
  echo 'Usage: run_all.sh resume | rebuild_pseudobulk' >&2; exit 2
fi
Rscript --vanilla "$mp10_scripts/02_models_worker.R" 0 2
Rscript --vanilla "$mp10_scripts/02_models_worker.R" 1 2
Rscript --vanilla "$mp10_scripts/12_finalize_models.R"
Rscript --vanilla "$mp10_scripts/05_prepare_networks.R"
Rscript --vanilla "$mp10_scripts/06_abundance_confounding.R"
if [ ! -f "$mp10_results/tables/epi_within_cluster_candidate_changes.tsv" ]; then Rscript --vanilla "$mp10_scripts/04b_T_training_cluster_candidates.R"; fi
Rscript --vanilla "$mp10_scripts/07_activities.R"
Rscript --vanilla "$mp10_scripts/08_gsea.R"
Rscript --vanilla "$mp10_scripts/09_nichenet.R"
Rscript --vanilla "$mp10_scripts/10_candidate_tables.R"
Rscript --vanilla "$mp10_scripts/19_prior_expression.R"
Rscript --vanilla "$mp10_scripts/17_integrate_evidence.R"
Rscript --vanilla "$mp10_scripts/20_additional_audits.R"
Rscript --vanilla "$mp10_scripts/28_nichenet_source_annotations.R"
# Training tables are frozen, source-hash-verified inputs; never overwrite a changed network silently.
if [ ! -f "$mp10_results/objects/scenic/Epi_training.tsv.gz" ]; then Rscript --vanilla "$mp10_scripts/04_scenic_training.R"; fi
for mp10_lineage in Endothelial T_lineage Epi; do
  "$mp10_python" "$mp10_scripts/24_grn_resume.py" "$mp10_lineage"
  "$mp10_python" "$mp10_scripts/13_scenic_ctx.py" "$mp10_lineage"
  if ! Rscript --vanilla "$mp10_scripts/14_aucell_full.R" "$mp10_lineage"; then echo "AUCell process exited nonzero; verifying completed outputs" >&2; fi
  "$mp10_python" "$mp10_scripts/29_verify_aucell.py" "$mp10_lineage"
  Rscript --vanilla "$mp10_scripts/15_scenic_summarize.R" "$mp10_lineage"
  if [ "$mp10_lineage" = Epi ]; then Rscript --vanilla "$mp10_scripts/15_scenic_summarize.R" Epi M20-masked; fi
  "$mp10_python" "$mp10_scripts/15b_scenic_rss.py" "$mp10_lineage"
  "$mp10_python" "$mp10_scripts/16_scenic_second_seed.py" "$mp10_lineage"
done
Rscript --vanilla "$mp10_scripts/23_regulatory_integration.R"
Rscript --vanilla "$mp10_scripts/30_evidence_details.R"
Rscript --vanilla "$mp10_scripts/31_gene_domain.R"
Rscript --vanilla "$mp10_scripts/32_axis_comparison.R"
"$mp10_python" "$mp10_scripts/35_reference_coverage.py"
Rscript --vanilla "$mp10_scripts/36_scenic_eligibility.R"
"$mp10_python" "$mp10_scripts/34_validate_scenic.py"
"$mp10_python" "$mp10_scripts/37_patient_relationships.py"
bash "$mp10_scripts/replot_all.sh"
Rscript --vanilla "$mp10_scripts/21_validate.R"
"$mp10_python" "$mp10_scripts/27_delivery_manifest.py"
