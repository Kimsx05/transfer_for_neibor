# Reproduction order

From this derived run directory, use `bash code/python.sh code/SCRIPT.py ...`; the included wrapper is a byte-identical copy of the historical v5 launcher. Existing run scripts retain the absolute historical launcher path used during execution. This reuses the actual installed cell2location environment without modifying historical outputs.

1. `01_prepare.py`: verify raw input SHA256, original split and revised scope; recompute scoped TRAIN filter and cache sparse inputs. Reuse only if fingerprints agree.
2. `02_reference.py MAIN`, then `02_reference.py NO_P07T`: Reference models can run sequentially; mapping of an available reference may overlap another reference fit under the later explicit GPU-concurrency authorization. Checkpoints are saved at training-segment boundaries; validate matching fingerprints before loading.
3. `03_pseudospots.py`, then `00_validate_inputs.py`: original scene/repeat scheme; preserves truth, source isolation and RNA bookkeeping.
4. `04_mapping.py MAIN DEV --smoke`; then MAIN DEV and NO_P07T DEV. Full scenarios follow the smoke check.
5. `05_freeze_DEV.py`: wait for both DEV fits, then freeze the source-group-balanced difficult-negative threshold and effective mapping settings before TEST.
6. `04_mapping.py MAIN TEST`, and NO_P07T TEST, concurrently only if measured free GPU memory supports both.
7. `06_evaluate.py`, `07_figures.py`: evaluate every generated scene, keep source-level results and plot data.
8. Review actual TEST outputs and record `08_review/ASSESSMENT.json` with conclusion, summary, A–G answers and `reviewed_actual_TEST_results=true`; then `08_report_package.py` validates outputs and builds report/ZIP. Do not fabricate this assessment before results exist.

GPU concurrency is explicitly authorized in RUNTIME_CONCURRENCY_OVERRIDE.md; measured memory thresholds and the actual remaining-run scheduler are in run_remaining_parallel.sh. Do not relabel excluded epithelial cells, reselect core K32, adjust threshold with TEST truth, or launch real spatial fitting. The seed manifest uses SHA256 purpose-derived values from master123. Resuming a partially completed training segment restarts that segment from its last saved model checkpoint with the recorded seed; completed fits are reused only when fingerprints agree. Model optimizer state is reset between separately invoked training segments by the underlying training API; this is recorded rather than claimed identical to uninterrupted Adam optimization.

DEV convergence amendment: the original3000epoch MAIN DEV stage is retained as ARCHIVE_MAIN_DEV_3000. The5000epoch probe showed that loss plateau alone was insufficient. The effective rule requires loss plateau plus mean-abundance stability for K32, Epi_other and their sum across segments, with256 diagnostic draws, median relative change<5%,q95<20% and denominator(previous mean+0.1), minimum5000epochs and cap30000. See06_evaluation/DEV_CONVERGENCE_REVIEW.json and DEV_FROZEN_CONFIG.yaml. No TEST truth selects the stopping point. The5000epoch probe script validates its exact archived3000epoch baseline and is not a second parameter scan.

`09_evidence_tables.py` rebuilds current support and MP input retention tables. `07_figures.py --reference-only` may run while GPU mapping is ongoing. `10_resource_monitor.py` records GPU usage until both TEST maps finish. Formal outputs always come from active MAIN_DEV/NO_P07T_DEV/MAIN_TEST/NO_P07T_TEST directories; archive/probe outputs never determine final metrics directly.
