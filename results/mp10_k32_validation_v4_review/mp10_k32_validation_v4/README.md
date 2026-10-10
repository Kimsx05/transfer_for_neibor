# K32 validation v4

Open REPORT.html first. K32=C18 of frozen FULL_Z_r1.00 (3282/123699 cells).

- 01_score_comparison: all27 cluster summaries/violins, raw-scale sample-paired tests.
- 02_rna_DEG: full rest/C17/C7 reused; full C14 and exclude_P07T rest/C14 newly fitted. Complete DEG, ordered top20/top50 and full GSEA in every model directory.
- 03_sample_robustness: all specimens/source coverage, eligibility, QC and fixed-panel sample expression.
- 04_14MP_profiles: raw and frozen-global-Z means, ranks and sizes.
- 05_umap: frozen FULL_Z/RNA, ALL and every dataset; per-cell plot data.
- 06_supplement: figure index, English captions, fixed panel and functional summaries.
- 07_review: concise full/exclusion and model-source tables.

Execute code/run_all.sh to resume with matching inputs/configuration. Do not source v3 scripts. Historical source objects and large aggregate caches are not part of the review ZIP; paths/checksums are recorded in 00_audit/INPUT_MANIFEST.tsv. All per-cell UMAP data and complete DEG/enrichment tables needed for review are included.
