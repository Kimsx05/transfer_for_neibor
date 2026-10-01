# Module 04 execution and review

This module depends on Module 00 only. The current explicit user request authorizes Module 04; the historical Module 00 contract and its 00/01 authorization metadata remain unchanged.

Code lives in `04_wang_niche/code/`. Run from the workspace with the existing `/usr/bin/python3.12` interpreter; `prepare.py` imports the accepted project Python libraries and isolated runtime without altering the base environment. R uses installed scran/igraph. All dataset jobs run serially; four candidate workers run within that one dataset job. BLAS/OpenMP and graph threads are one, with no nested parallelism; GPU is disabled. The total CPU budget is eight. Identical serial-versus-parallel memberships and silhouettes are verified on the two completed baseline datasets.

Execution order:

1. `python3.12 code/run_module04.py --phase prepare`: validate input SHA, preserve all frozen spot keys, normalize all 30 accepted posterior means, save the fixed evaluation sets and parameters.
2. `python3.12 code/run_module04.py --phase smoke`: run the complete GSE171351 dataset first, including all 45 candidates, annotations and figures.
3. Inspect the smoke thumbnails. Then `python3.12 code/run_module04.py --phase all`: reuse only signature- and checksum-matching completed outputs; run the remaining datasets serially.
4. Record visual checks for each dataset in `visual_review.tsv`.
5. `python3.12 code/finalize.py`: independently verify selection rules, complete spot keys, unchanged inputs/upstream contract, numerical tests and PDF/600 dpi PNG files; create the tables, review bundle, manifest and handoff.

The command examples above are relative to this module directory. From the workspace, prepend `spatial_ecology_top20_v2/04_wang_niche/` to script paths.

The review archive contains summary tables, code, logs, thumbnails and representative PDFs. Full-resolution figures and full spot matrices remain outside the compact archive. `review.html` provides a gallery. `figure_index.tsv` points to each actual PDF/PNG/thumbnail; `output_manifest.tsv` provides file SHA-256 values.

Statistical boundaries and deviations from the publication are in `methods_parity.tsv`. Spot-level Wilcoxon tests define the clusters' composition; they are not independent biological validation or patient-level inference. No MP–niche tests or later modules are run.
