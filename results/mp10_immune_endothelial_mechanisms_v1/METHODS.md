# MP10 immune/endothelial mechanisms: frozen methods

## Scope and input

This is an MP10 analysis, not an LLPS/C5-CAF analysis. The source is `subtype/v2/final_annotation_pre_reference/final_annotation_v2_pre_reference.qs`, with RNA Assay5 `counts` and original normalized `data` layers. The input audit confirms 90,011 features and 384,612 source cells. Only exact members of the final c2l 30-class reference and the 40 originally eligible Tumor specimens from baseline `c2lfinal_20261108T123140Z_0b5fd64f` are used. The frozen MP10 scores and Top20 are reused without rescoring. CellChat's signaling cache and capped membership are not used for pseudobulk.

Full specimen Epi score distributions define L50 <= Q50 and H50/H10/H05 > Q50/Q90/Q95, using R type-7 quantiles without splitting ties. There are 116,448 original eligible Epi cells. High pools overlap; the three contrasts and sensitivities are not independent replication.

Input file paths, sizes, timestamps and SHA256 values are in `tables/input_manifest.tsv`. Barcode, feature, score and membership checks are saved separately. Raw counts are verified nonnegative and integer-valued; no counts are reconstructed from normalized data. The original Seurat object is read in place and is not copied.

## Specimens, recorded patient units and eligibility

The confirmed T4463209a/b pair is handled by retaining a in every primary inference and replacing a with b in the prespecified replacement analysis. Failure of a to meet a receiver's cell threshold does not trigger substitution of b. Both specimens can enter descriptive coverage and the shared regulatory-network training set. Other patient identifiers remain recorded, derived/provisional identifiers; their biological independence is not assumed verified. Dataset-prefixed identifiers are never merged merely because their short names match.

Pseudobulk eligibility is **n_cells >= 50**, with **>= 30** sensitivity. Each receiver has its own eligible specimens. Each Epi contrast has its own eligible pairs; an additional common H50/H10/H05/L50 set is used for comparable trends. These thresholds differ from CellChat's node-filter boundaries. Models require at least four recorded patient units, a full-rank design and at least three residual degrees of freedom. Unavailable models are NE, not zero effects.

## Pseudobulk and DE

Raw UMI counts are summed over all eligible cells in each specimen × Epi state and specimen × final TME type. Epi_all profiles support exposure-independent sender expression checks. The reusable compact object retains counts, per-gene positive-cell counts, sums of original normalized expression and profile metadata. Overlapping Epi profiles retain their actual members; there is no equal-cell sampling or cap of 100.

edgeR 4.8.2 is used with filterByExpr on the actual design, TMM normalization, robust dispersion estimation and robust quasi-likelihood inference. Epi uses `~ specimen + state`; adding the nested dataset term would be redundant. Receiver models use `~ dataset + exposure` (or `~ exposure` when only one dataset is represented). The exposure is the Epi median; Q95 is a separate sensitivity. Both are standardized by their respective SD across the original 39 main-a specimens, frozen across every subsequent model. Coefficients therefore refer to one fixed original-specimen SD, not a dataset-specific or subset-specific SD.

Per-comparison BH FDR is retained. Joint primary BH q values are calculated separately across the three Epi contrasts and across the four receiver model families, using all testable genes. Sensitivity analyses are separate families; no minimum P over settings is selected. Numerically tiny negative edgeR F values with P=1 are retained in the original DE output; only the downstream signed ranking is defined as `sign(logFC)*sqrt(max(F,0))` to avoid NaN rankings.

The model plan includes n>=30, a/b replacement, Epi common-set contrasts, receiver Q95, leave-one-dataset-out, and descriptive Epi per-dataset models. Model-specific memberships and design rows are delivered. Cached inference is reused only after checking frozen inputs, design, parameters and gene universe. A documented migration retains completed serial inference and reconstructs its normalization matrices without repeating the statistical test.

## Expression, abundance and confounding

Candidate expression uses actual UMI, original normalized mean expression and the fraction of cells with RNA counts >0. A gene absent from the source feature universe differs from an observed zero. Receiver ligand/receptor RNA ratios are not receptor occupancy or evidence of competitive dominance.

Abundance uses uncapped cell counts. Treg and Tfh use total T/NK as the primary denominator; Arterial and Capillary EC use total Endothelial. Total TME is secondary. Sample-level quasibinomial models account for overdispersion. Absence of captured subtype cells is an actual numerator of zero when the denominator is observed; no denominator cells is NE. The abundance model is distinct from RNA and TF activity models.

Existing `epi_allgene_cluster` labels are used as expression clusters, not newly inferred malignant states or biological subtype names. Major clusters are frozen as >=5% of original Epi cells. Candidate expression is examined within these clusters with paired specimen summaries at >=50/30 cells per state and paired Wilcoxon tests when at least four pairs are available. Quality differences and within-specimen MP10–depth/features correlations are reported. No new CNV calls are made.

## Functional and TF analyses

MSigDB human 2026.1 Hallmark is the primary GSEA family; Reactome is secondary. GSEA uses all finite signed testable-gene statistics, fgseaMultilevel, minSize=15, maxSize=500, nPermSimple=1000, eps=1e-10, seed=123 and one worker. Pathway mapping counts, NES, P, BH q, leading edges and numerical warnings are retained. The same deterministic gene ordering resolves tied ranking statistics; no target genes are added to favor a hypothesis.

CollecTRI is the frozen human OmniPath protein TF–target network, retaining AP1/NFKB complexes using the decoupleR 2.9.7 convention. The downloaded network, provenance and hashes are retained; it is not silently refreshed. ULM is run through the installed official decoupleR package with at least five mapped targets. PROGENy uses the frozen human full model's 500 lowest-P targets per pathway and ULM with at least ten mapped targets. DoRothEA A/B provides a separate high-confidence sensitivity, also requiring five targets.

Activity scoring uses a fixed gene universe per cell-family, filtered from profiles with at least 30 cells, TMM logCPM (prior.count=2), and gene-wise centering across that family's profiles. This is a downstream activity input; raw counts, libraries and MP10 definitions remain unchanged. Sample-level activity models use the same specimen-paired or dataset-adjusted designs as DE, implemented with limma and robust empirical Bayes. `sample_model_p/q` are sample-level association results. `ulm_enrichment_p/q` are internal per-profile target-enrichment results and are never presented as between-patient tests. TF RNA and regulon activity are reported separately.

For Epi, Full and M20-masked pathway/network analyses are both retained. M20 masking removes the frozen Top20 from the downstream ranking/network and scoring universe only. Receiver primary results retain all genes. Hallmark, Reactome, CollecTRI, PROGENy and DoRothEA have separate correction families; primary pooled corrections and sensitivity results are explicitly labeled.

## NicheNet

NicheNet 2.2.1.1 uses the human v2 prior from Zenodo 7074291: `lr_network_human_21122021.rds` and `ligand_target_matrix_nsga2r_final.rds`. The latter has 33,354 targets × 1,226 ligands. The audit checks literal names and the NECTIN2/PVRL2 nomenclature relationship, LR entries, target-matrix columns and observed expression. A missing LR, missing ligand-target column, undetected receptor, insufficient receiver signal and insufficient sample coverage are distinct reasons.

Epi→TME targets come only from that receiver's own dataset-adjusted continuous MP10 association, within-model BH FDR<0.05, with positive, negative and union sets analyzed separately. Background is the receiver's complete testable gene universe intersected with the ligand-target matrix. Targets must be subsets of background. Empty targets are NE; 1–9 targets are explicitly low-information exploratory results. Thresholds are not relaxed to produce a favorite ligand.

Receptor-qualified sender-agnostic and sender-focused candidates are both retained. The expression rule is median per-eligible-specimen positive-cell fraction >=10%, with 5% sensitivity. Per-specimen and cell-source coverage accompany the rankings. This does not impose a post-hoc requirement for a fixed number of datasets. All suitable ligands compete; the pool is not restricted to NECTIN2 or VEGFA.

Native NicheNet ligand activities are sorted by aupr_corrected, with AUROC, AUPR and Pearson retained. Focused ranks are derived from the same activity estimates. Ligand activity is not an LR patient-level P/FDR. Prior weights and actual target logFC/FDR remain separate. Neither a high NECTIN2 rank nor high TIGIT RNA identifies the dominant receptor branch or proves immune suppression/exhaustion.

The supplementary TME→Epi analysis uses the actual TME as sender and Epi paired High–Low targets. Full and M20-masked target AND background sets are retained; M20-masked is the primary reverse-direction ranking. The same TME cells are never duplicated into invented High/Low environments. Shared ligand/TF targets are reported as overlapping evidence, not a proven ligand→TF causal path.

## pySCENIC

An isolated Python 3.10 environment contains pySCENIC 0.12.1 and pinned compatible dependencies. Three shared networks are learned: Epi, T lineage (T/NK major class with NK excluded), and Endothelial. Within each lineage, seed=123 random priorities select at most 100 cells per specimen, independent of MP10 and candidate expression. Genes must be detected in >=1% and >=3 training cells; selection is not restricted to Seurat HVGs. Original RNA data, not integrated values or Harmony coordinates, is used.

GRNBoost2 is followed by cisTarget motif pruning and AUCell on the full eligible lineage cell set. Human hg38/refseq-r80/mc9nr ranking and matching v9 motif annotations are frozen. The official resource host encountered TLS failures; same-version archived files from Zenodo 7562143 are accepted only after matching their published MD5. Official Feather format conversion does not change rankings. Training memberships, features, seeds, reference hashes, regulon targets/coverage and exact command parameters are retained.

A prespecified seed=456 checks key regulatory connections; unverified connections are labeled accordingly. RSS is subtype specificity, not an MP10 association P. Motif support is not direct binding, AUCell is not TF RNA, and network plots are not ChIP validation.

## Resources and interpretation

The requested limit is one BLAS/OMP/MKL thread per worker and four total computational workers. A runtime audit discovered that this host’s global Rprofile reset early R jobs to 24 threads despite shell settings. The deviation and affected stage are recorded in `logs/resource_constraint_incident.json`; subsequent R jobs bypass the site profile and explicitly reset thread libraries, with CPU affinity as an additional bound. The early deviation is not represented as compliance. Existing HD jobs and prior CellChat results are not changed. Public references are shared, original large data are not copied, and no archive contains original Seurat/counts objects.

Spatial HD support is pending. Repeated settings, overlapping specimens and multiple algorithms applied to the same expression data are not independent biological validation. Negative and not-evaluable modules remain in the report.

Primary method sources: [NicheNet](https://github.com/saeyslab/nichenetr), [NicheNet protocol](https://www.nature.com/articles/s41596-024-01121-9), [decoupleR TF workflow](https://saezlab.github.io/decoupleR/articles/tf_sc.html), [pseudobulk workflow](https://decoupler.scverse.org/en/stable/notebooks/scell/rna_psbk.html), [SCENIC protocol](https://github.com/aertslab/SCENICprotocol), [cisTarget format tools](https://github.com/aertslab/create_cisTarget_databases).

## Implementation validation and numerical limits

The NicheNet public ligand-activity wrapper also computes several metrics that it does not return. The delivered implementation computes exactly its four returned metrics using the same ROCR precision/recall, caTools trapezoidal AUPR, native split-AUROC and Pearson definitions, with the same background and target ordering. It was checked against the installed native function on four cached Full/masked positive/negative target sets (maximum absolute difference <=5.6e-16) and an actual two-target receiver set (exact equality). Validation tables and the native function source snapshots are delivered. This optimization does not change the prior, candidate domain, targets or ranking definition.

With very sparse targets, constant prior predictions can receive misleadingly high native trapezoidal AUPR because of tie handling. Native scores remain unchanged, but missing Pearson correlations/constant predictions and low-information target sets are flagged. Such ranks are not interpreted as mechanistic winners.

The pySCENIC second seed is a **key-connection** check: top five CollecTRI TFs per primary model by q rank, regardless of significance, define the union of their motif-pruned targets; named hypothesis molecules and MP10 Top20 are also included. GRNBoost2 seed456 uses these targets with all original TF predictors and the unchanged training matrix. It is not a second full motif-pruned network, a second AUCell experiment, or independent biological validation. Coverage and recovered connections are retained. Epi AUCell additionally has an M20-masked sensitivity with Top20 removed from regulon targets and downstream ranking background; original expression, raw counts, MP10 scores and groups remain intact.

The methodological example [Xu et al., 2025](https://doi.org/10.1002/advs.202502774) is used only for the GRNBoost2 → motif filtering → AUCell analysis hierarchy. Its CAF groups, patient counts and biological conclusions are not transferred to this MP10 project.

The complete native NicheNet ranks are preserved even when a ligand has constant predictions across the tested background. `rank_variable_prior_*` is an explicitly labeled diagnostic rank among nonconstant prior predictions, never a replacement for the native table and never a statistical significance claim. Low-information receiver target sets remain low information under either ranking.

Abundance BH families are the four receivers within each fixed denominator type, minimum denominator count, patient selection, exposure and omitted-dataset setting. A denominator threshold does not require the numerator subtype to have 50 cells; zero observed subtype counts remain valid abundance observations when the denominator is observed. This is why abundance coverage differs from receiver RNA coverage.

All testable genes appear in the DE statistics; `de_gene_filter_exclusions.tsv.gz` explicitly lists every other feature in each completed model and distinguishes all-zero model counts from nonzero features failing filterByExpr. Full per-profile positive-cell fractions can be reconstructed exactly as `pseudobulk$detected_cells / profiles$n_cells`; these are not inferred from normalized RNA.

NicheNet ligand activities are attached to candidate LR rows for navigation only. The ligand–target matrix does not identify a receptor-specific branch; `NicheNet_activity_receptor_branch_specific=FALSE` and the receptor's own detection flag are retained. A ligand qualified by TIGIT or CD226 expression cannot be treated as evidence for a low-detected PVRIG branch.

Motif pruning initially used one worker per lineage. The two incomplete Epi/T motif jobs were restarted with two single-thread workers each to use the available four-worker budget; completed GRNBoost2 outputs were unchanged. The actual execution worker counts are in each ctx metadata file and the scheduling-change log. Worker count is treated as an execution parameter when validating scientifically identical ctx inputs for resumption.

The completed Endothelial AUCell/seed456 jobs used one worker. Epi/T AUCell and seed456 use at most two single-thread workers per lineage within their assigned two-core affinities. At most two lineage jobs run together. Official AUCell/GRNBoost2 interfaces and the frozen seeds, genes, cells and regulons are unchanged; worker counts are execution metadata, not extra biological replicates.

Second-seed summaries report edge presence/recovery and the descriptive Spearman correlation of positive GRN importance values available in both seeds. No post-hoc recovery cutoff defines a validated connection; this technical consistency measure is separate from patient-level associations and motif support.

A regulon may retain its TF itself as an AUCell signature member. Such a self-target is not a GRNBoost2-inferred self-regulatory edge. Raw signatures and AUCell remain unchanged; self-target members and any edge lacking an original seed123 GRN weight are excluded from the seed-recovery denominator and from displayed regulatory arrows. Their reasons are retained separately.
