# Reproduction and resume

Run from the Bladder Metabolism project root. Use `bash mp10_k32_reference_diagnostic_v6/code/python.sh <script> [arguments]`; the wrapper uses the existing project environment, not system packages. All large source paths are in00_manifest/INPUT_MANIFEST.tsv. Source v5 files remain read-only.

1. 01_prepare.py validates inputs and labels. 13_index_reused_audits.py records reusable scope audit copies and seeds.
2. 02_merged_reference.py fits/resumes the one S30 NB. A valid REFERENCE_DONE and matching fingerprint prevents duplicate work.
3. 03_paired_control.py freezes/generates the216 exact DEV bags and both648-spot counts; CONTROL_FROZEN is authoritative. Do not overwrite frozen generated counts unless deliberately reconstructing and checking every hash.
4. 04_mapping.py takes `S30 FULL_DEV`, `S31 TARGET_UMI`, `S31 FIXED_CAPTURE`, `S30 TARGET_UMI`, `S30 FIXED_CAPTURE`. S31 fullDEV reuses old final output. Model/input/signature/seed fingerprints are checked. `run_remaining.py` is the resource-aware resume scheduler; it detects existing processes and completed tags. Do not launch duplicate tags manually.
5. 05_feature_audit.py audits local inventories. 05b_fetch_annotations.py optionally obtains small public annotation resources;05c_close_provenance.py incorporates their actual availability and official metadata. Download failures remain explicit; no G1 model is fitted. NEXT_G1_CONFIG is only a candidate plan.
6. Once all five maps have final DONE records:06_evaluate.py,07_reference_and_history_audit.py,09_additional_evidence.py,08_figures.py,10_extra_figures.py,11_validate_delivery.py. `finish_analysis.py` waits for existing jobs and runs that sequence. It does not create extra models.
7. 08_review/INTERPRETATION.json records the result-based scientific review. 14_slim_plot_data.py retains only plotted variables and exact group-point values; full core per_spot_all remains. After delivery checks pass,12_report_package.py generates REPORT.html and review.zip. This step never changes thresholds or fits models.

For direct checkpoint reconstruction, read the matching `run_fingerprint.json`; reload its counts path and exact gene/signature order. S30 reference: recreate label_S30 by replacing only Epi_K32/Epi_other with Epithelial; setup_anndata(labels_key='label_S30',batch_key='sample'). Spatial models: setup_anndata(batch_key='sample'); load cell2location checkpoint with this same registered data. No scenario/truth covariates. w_sf draws remain on the server with hashes; no other large latent draws are required.

The initial inherited CONFIG simulation block describes historical fullDEV. New control repeats/scenarios/probabilities are in CONTROL_FROZEN. GPU parallel authorization is in RUNTIME_CONCURRENCY_OVERRIDE; it overrides only resource scheduling. Existing checkpoint segments retain the old training behavior. No real slice or TEST fitting is authorized by these scripts.
