# Actual methods — K32 diagnostic v6

This run adds one S30/G0 NB reference, one fullDEV S30 map and four small paired DEV maps. Old S31 MAIN/NO_P07T DEV/TEST are read-only. No new TEST, G1 or real spatial fitting is performed.

## Scope and configuration

The unchanged scope is281,851 cells:123,699 frozen epithelial cells (3,282 K32;120,417 other) plus all old non-epithelial members. The18,182 excluded cells remain outside_frozen_discovery_scope, not MP10-low/QC-fail/K32-negative. Legitimate Cycling-derived members remain. The prior BLOCKED record reflects a scope-rule conflict, not biological failure. MAIN TRAIN contains135,860 cells including2,395 K32. Original conservative source groups25/8/8 are unchanged; patients are not independently verified. P07T/P07N stay TRAIN. No balancing or reference downsampling.

S30 merges only Epi_K32/Epi_other into Epithelial in a new label column; other29 labels and cells are unchanged. K32 identity stays in truth metadata, never a prediction covariate. G0 is the exact old11,230-gene order. Sparse raw integer RNA caches are reused by path/SHA256; no normalization is substituted. Inhouse MT remains unreliable.

CONFIG.yaml inherits the effective v5 DEV-frozen settings. Its simulation block describes old fullDEV; new controls use CONTROL_FROZEN.yaml (old repeats0/1, not10 new repeats). Actual label orders are labels_S30/S31.json. The later user-authorized RUNTIME_CONCURRENCY_OVERRIDE replaces only the sequential-GPU sentence. NEXT_G1_CONFIG is future-only.

Master seed123; derived seed=int.from_bytes(SHA256('123|purpose')[:4],'little') modulo(2^31-1). Purposes/numerical values are saved. S30 NB uses old reference_fit seed; fullDEV uses old mapping_DEV; all four controls share v6_control_mapping and v6_control_posterior seeds. Bag seeds are frozen before fitting. No seed search. Existing cell2location0.1.5/scvi1.3.3/torch2.6 environment is reused; environment.json is authoritative. Up to3 GPU jobs after resource checks and explicit user authorization; other tasks untouched. Each process uses1 CPU thread.

## NB and signature scale

RegressionModel: every MAIN TRAIN cell, label_S30, sample batch, no extra covariates;250 initial epochs,125 extensions, cap1000; batch2500, lr0.002, train_size1. Adjacent50-epoch loss means must differ by<1%. S30 actually stopped at250. Reconstruction uses a deterministic1000-cell TRAIN subset for QC only. No earlier signature checkpoint exists, so interstage signature stability is not claimed.

Signatures=posterior mean per_cluster_mu_fg.T times TRAIN-cell-weighted mean batch detection_mean_y_e, matching the old project. This module computes detection_y_c deterministically as sample-design times detection_mean_y_e. Installed export scales only if detection_y_c appears among exported sites; the project explicitly reconstructs that mean. Saved signatures are checked against this formula, with actual cell versus equal-batch weights. No truth-based rescaling. Remaining29 signature correlations/sums are descriptive.

## Exact paired controls

216 unique bags from frozen DEV, old repeats0/1: S1 epithelial10/K32=0,1,2,5,10; S3 epithelial20/K32=0,10,20; S2 epithelial10/K32=0 with C14,C7,C17,UCell-high non-K32,NMF-high non-K32. All feasible sources retained; unsupported combinations listed. Every bag has20 distinct cells from one sample. No cross-sample or TRAIN/TEST borrowing. Both rules/depths share exact cell IDs; across-bag reuse recorded. Frozen scores, clusters/high flags are unchanged.

Full RNA matches old per-cell raw_UMI, not just G0. Unique-bag median198,515.5 UMI. Anchors3146/7005/14967 (low/medium/high). TARGET_UMI p=min(1,T/S_b). FIXED_CAPTURE p=0.015847629026448817/0.03528691714249013/0.07539461654127763, constant across sources/bags within depth. All feasible; no clipping/upscaling/post-thinning renormalization. Each rule has648 spots.

Within a bag, probabilities are ordered descending; first draw at largestp, then nested binomial thinning by adjacent probability ratios. This preserves Binomial(n,p) marginals and couples comparisons. Sparse G0 counts and the non-G0 per-cell full-RNA remainder are sampled separately, giving an equivalent remainder total without inventing gene columns. Cell/subtype retained RNA is saved. Truth counts stay unchanged even if a cell retains zero RNA. Newly generated counts/smaller mapping collection do not numerically reproduce old predictions. S31/S30 use byte-identical input per rule. Input obs is sample only: no truth/scenario/depth/p/MP score covariates.

## Mapping and posterior

N_cells_per_location20 is a prior, not a hard total. detection_alpha20, sample batch, lr0.002, full batch, train_size1. Initial3000, extensions2000, minimum5000, cap30000. Loss adjacent500 means relative change<0.005 AND abundance stability:256-draw checkpoint means, abs(new-old)/(old+0.1), median<0.05 andq95<0.20. S31 checks K32/other/sum; S30 checks Epithelial with identical numerical cutoffs. Actual epochs/statuses in MODEL_INDEX.tsv. Nonconvergence is not biological failure.

Final1000 w_sf draws plus mean/median/q05/q95. Sums and ratios computed within draw before quantiles. For K32/epithelial, exclude draws with denominator<0.1; if>5%draws fail, ratio for the spot is NA; quantity metrics retain the spot. S30 K32 predictions are inapplicable, never zero or truth-allocated. No normalization to20 or truth-based calibration.

Per-location detection_y_s has a batch hypermean; it is not one fixed sample coefficient. Technical summaries use100 posterior means. Expected RNA=((w_sf @ signature.T)*m_g + batch s_g_gene_add)*detection_y_s, checked against installed reconstruction on256 locations. Plug-in means are not joint posterior expectations/intervals. No spot×gene×draw array is saved; old checkpoints are not retrained.

## Evaluation

44 old quantity/ratio checks reproduce v5 sample-equal metrics. New primary weighting: equal source groups, equal samples within group, equal spots within sample. Weighted RMSE=sqrt(weighted squared error); slope from weighted covariance/variance, not averaged sample slopes. Sample-equal alternative retained. Ratio weights are renormalized over evaluable spots. Conditional FPR/recall/precision average eligible samples/groups; undefined denominators explicitly counted. Simulated prevalence accompanies PR/AUPRC/precision and micro/defined-source macro results.

MAE/RMSE/bias, calibration,90%coverage/width, counts versus ratios, all29 non-epithelial labels/major mappings, and quantity/depth/scenario/sample/group/dataset strata are exported. Constant truth, including total20, has no slope/correlation. Differences paired per bag/spot before source aggregation. No per-spot significance tests or false patient-level confidence intervals. Historical threshold is exactly8.367393, strict posterior-mean>; never refit. Zero q05 is not a probability criterion. K32 recovery differs from MP10-high coverage. Existing S4 uses old fullDEV/TEST only; RNA-depth confounding retained, no new core definition.

Old TEST was already viewed: historical descriptive evidence, not fresh independent validation. No TEST-driven genes/bags/probabilities/priors/stopping selection.

## Features and G1

Actual preflight inventories, local processing scripts and official GEO metadata are separate evidence layers. Six small GSE285715 deposited feature lists were fetched; no spatial expression matrices. HSD17B7/KRT18 already absent there, but deposited feature absence does not prove assay impossibility. Public fresh-frozen/mm10 metadata conflicts with human feature IDs. GSE238145 B5/B9 document distinct FFPE probe versions. Official probe CSV downloads failed; raw uncertainty is not called an ID bug. URLs/access dates and per-section limitations retained.

G1 candidates use exact target observed features grouped by documented chemistry where available, MAIN TRAIN-only permissive filter, canonical IDs and MT exclusion. No absent gene expression is zero-padded. No substitute MP genes, outcome-based panel choice or G1 fitting.

## Limits and package

K32 came from the existing full discovery set: this is frozen-candidate transfer, not independent rediscovery. Source groups are not verified patients. Pseudospots cannot fully model scRNA-to-Visium effects. No real slices fitted. Technical completion, convergence, recovery and spatial validity are separate conclusions. Review.zip contains results/figures/plot_data/configs/code/provenance/index; counts/checkpoints/large w_sf draws remain on server.
