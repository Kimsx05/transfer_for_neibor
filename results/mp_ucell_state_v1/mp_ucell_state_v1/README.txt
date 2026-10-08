14 frozen MPs / UCell state clustering — mp_ucell_state_v1

New outputs only. Original NMF, annotations, embeddings, prognosis and spatial analyses are read-only.
Configuration: config.R and config.sha256. Frozen input: output_full_cc_seed5_v2.
Canonical basis: 07_EpiMP_gene_basis.qs$sample_balanced_mean, verified against CSV.
FULL: exact original input cells and 01_input_state.qs$cells; no malignant-only or atlas expansion.

Run from /home/data/t070721/codex_workspace/Bladder Metabolism:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript mp_ucell_state_v1/code/01_audit_score.R
  bash mp_ucell_state_v1/code/run_remaining.sh
  Rscript mp_ucell_state_v1/code/07_input_qc_figures.R
  Rscript mp_ucell_state_v1/code/08_covariate_figures.R
  Rscript mp_ucell_state_v1/code/11_finalize_cluster_colors.R
  Rscript mp_ucell_state_v1/code/12_finish_render_validation.R
  Rscript mp_ucell_state_v1/code/05_validate.R
  python3 mp_ucell_state_v1/code/10_validate_artifacts.py
  python3 mp_ucell_state_v1/code/09_interpretation.py
  python3 mp_ucell_state_v1/code/06_report.py

Stage scripts bootstrap the existing R 4.5 library, not the shell's R 4.4 user library.
Score once on sparse complete RNA counts, 14 plain positive top50 gene sets together.
UCell 2.14.0, maxRank1500, average ties, missing_genes=impute, chunk500, workers2.
All graph instances consume numeric cells x 14 score matrices directly; no PCA/Harmony/smoothing.
Leiden backend: leidenbase / RBConfigurationVertexPartition (configuration-model modularity).
FULL_Z is primary; balanced seed42 is primary display. All fixed seeds and resolutions retained.
BAL is cap sampling by biological specimen, not equal sample or patient weighting.
Patient mappings recorded as derived/provisional are audit-only, not verified patient identities.
FULL proportions describe FULL; balanced proportions describe only the selected cells.
MP10 high clusters are descriptive because MP10 participates in clustering.
Agreement with original usage is method linkage, not external validation.

Folders:
00_input_audit: inputs, cell identity, mappings, RNA layers, installed API, validation.
01_signatures: top50 with loading/rank/coverage/hash; frozen comparison; overlap/Jaccard.
02_ucell_scores: unique M_RAW.qs, M_Z.qs, FULL scaling, depth/detection and vector QC.
03_full_z / 04_full_raw: independent sparse graphs, labels, embeddings, summaries and figures.
05_balanced: exact manifests for seeds42/43/44 and independent Z/RAW instances.
06_comparisons: common-ID ARI/NMI, overlap, normalized overlap, best Jaccard, 14MP profiles.
07_mp10_review: all-cluster candidates, usage and depth correlations, interpretation limits.
08_figures_by_dataset: subsets of each global embedding; no dataset-specific UMAP fit.
code / logs: scripts, execution logs, sessionInfo, stage completion sentinels.

Methods documentation:
https://www.bioconductor.org/packages/release/bioc/vignettes/UCell/inst/doc/UCell_parameters.html
https://satijalab.org/seurat/reference/findneighbors
https://satijalab.org/seurat/reference/findclusters
https://satijalab.org/seurat/reference/runumap
Framework inspired by the user-specified Hu et al. 2026 description; parameters are project-prespecified.
This analysis does not claim exact reproduction of unpublished author code.
