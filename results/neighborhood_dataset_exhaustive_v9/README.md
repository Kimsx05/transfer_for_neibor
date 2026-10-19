# MP10 neighborhood exhaustive dataset analysis v9

This run enumerates all 99 subsets of the frozen seven dataset identifiers with sizes 3–7 (57 with inhouse, 42 without). The all-seven pilot is combination `1111111` and is counted once. Bit positions follow `config.json`; subsequent order is dataset count then bitmask, independent of P/q. The immutable registration remains REGISTERED; live statuses are in `combination_summary.tsv` and each combination's stage manifests.

All prior projects are read-only. Shared inputs are hashed references to v2 original float32 q1/q2 Hellinger features, v7 float64 ring arrays and scores, v8 full-tissue geometry and cell-target hotspots. No raw expression or cell2location objects are copied. v3 representation code runs with its existing numerical runtime; SpottedPy uses the existing v8 environment. The locked source commit is 9dc364f203b7617b6701aa8fb579cbbb35357b26; esda 2.9.0 and libpysal 4.13.0. Pure audited functions are loaded via AST without prior projects' import-time output paths.

## Execution commands

Run/resume the complete serial workflow (default one combination and one section at a time):

```bash
cd '/home/data/t070721/codex_workspace/Bladder Metabolism'
PYTHONDONTWRITEBYTECODE=1 neighborhood_spottedpy_v8/env/bin/python neighborhood_dataset_exhaustive_v9/code/run.py --resume
```

Resume one registered combination, retaining the same scientific settings:

```bash
PYTHONDONTWRITEBYTECODE=1 neighborhood_spottedpy_v8/env/bin/python neighborhood_dataset_exhaustive_v9/code/run.py --combo 1111111
```

Refresh summaries without fitting:

```bash
PYTHONDONTWRITEBYTECODE=1 neighborhood_spottedpy_v8/env/bin/python neighborhood_dataset_exhaustive_v9/code/run.py --report-only
```

After all 99 main analyses complete, generate/continue the selected review (also runs automatically):

```bash
PYTHONDONTWRITEBYTECODE=1 neighborhood_spottedpy_v8/env/bin/python neighborhood_dataset_exhaustive_v9/code/run.py --review-only
```

Regenerate a specified combination's reference representation for display; frozen labels and main results are untouched:

```bash
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 neighborhood_dataset_exhaustive_v9/code/refit.py 1111111 --rebuild
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 neighborhood_dataset_exhaustive_v9/code/umap_display.py 1111111
```

Only when the controller is not running, clear registered regenerable temporary exchange files:

```bash
PYTHONDONTWRITEBYTECODE=1 neighborhood_spottedpy_v8/env/bin/python neighborhood_dataset_exhaustive_v9/code/run.py --clean-temp
```

The cleanup checks real paths, ownership, exact registered file paths and a controller lock. It never recursively removes old projects or follows file symlinks. R requires a path without spaces; `/tmp/neighborhood_v9_rtemp_t070721_20261018` is only a symlink alias into this run's `tmp` directory, where all actual R temporary data remain. A newly started server session can restore this alias using the same exact target after verifying that an existing alias is correct.

## Scientific definitions

REFIT independently resamples at most 20,000 references by the original dataset/patient-stratum/section hierarchy, fits minimum-95% PCA without whitening, then joint SAMPLE/WEAK Harmony on all member anchors (theta 1, lambda 5, sigma .1, internal K100, seed 123, 20 iterations with the audited finite 50-iteration continuation). Euclidean self-including 30NN SNN uses Jaccard >1/15; weighted Leiden RB, 10 iterations. The finite resolution grid and audited nondegeneracy/positive-separation/moderate-K priorities select niches without MP10. Stability seeds 124–128 use balanced 80% references; a fixed spatial-block perturbation is additional. Transfer uses 15NN, reference native override, support ≥.6 and reference leave-self-out distance Q99. All combination labels and numeric results are retained.

The original manual v7 139 niche/type list is not transferable to redefined niches. Every combination therefore freezes all K×29 non-Epi pairs before accessing scores, plus 29 all-eligible-Epi association queries. R1/R2/R3 use common 4/8/12 supports; raw/excess/paired-ring definitions and rank-residualized correlations follow v7. V7 numerical zero tolerance is exactly zero, while v8 signed-rank uses 1e-12. Fisher-z section correlation averaging precedes unit-level tests.

Source E0 is reconstructed from geometric/score eligibility and the *current* accepted mask. Source cache keys include ordered UIDs, actual E0, scores, inputs and parameters. Prior source hotspots are reusable only after equality of E0, supported rows, scores and sparse W; otherwise source G is recomputed. Cell targets retain the complete-tissue v8 definitions. SOURCE_P05 and LOCAL_BH05 are separate tracks, with local BH over doubled native tails across all calculable positions. Formal hotspots require R1 components ≥5. O/D/D_out support thresholds are unchanged; unreachable target distances remain NA. No N×N tissue distances or retained permutation simulations are created.

Section effects are averaged equally within `analysis_unit_uid` across datasets before signed-rank inference. Inhouse identities remain the documented independent-sample policy, not seven verified patients. GSE238145 B17/B18 retain the traceable conflict policy. Full 63-section identity table and source limitations are in `audit`.

FROZEN_LABEL always plans nine old niches and uses original v8 E0 section effects. Its 99 primary-style control families contain exactly 3,564 relations. It is not independent validation and its source coverage can differ from REFIT.

Primary family: REFIT × SOURCE_P05 × k10 × ALL_POLICY_UNITS × both layers × all frozen niches × O/D, totaling sum(4 K_c). Within-combination BH and global BH/BY use raw P, with successfully evaluated but scientifically untestable planned relations occupying P=1 slots internally; their reported P/q remain NA. Technical failures are unresolved, never ordinary NA. Global values are PROVISIONAL until all required main stages complete. LOCAL_BH05, D_out, cell targets, verified-only, frozen labels and ring endpoints have separate named families.

## Files and disk budget

`primary_relations.parquet` and auxiliary Parquet tables contain complete tested and untestable planned results, not just discoveries. Each combination retains compact indices/labels, profiles, cross-label correspondence, coverage, section and unit effects, raw P and within-family q. Only final aggregate tables contain cross-combination q. No redundant large CSV/Parquet copies are created. Hit TSVs are small derived views.

Budget: 3 GiB maximum total new data; cache ≤25%, temporary exchange ≤15%, review ≤20%, with 2 GiB external free-space reserve. Quota is checked separately from filesystem capacity. Guards execute before heavy stages and source computations. `resource_usage.tsv`, `cleanup_manifest.tsv`, `RESUME.json` and stage hashes describe actual usage and restart state. Stage outputs are written atomically before completion markers. Input/config/required-artifact hashes are checked on resume; code changes after launch require an explicit technical repair audit.

Reviews are generated only after the entire primary scan. At most three combinations and three relations per combination are selected by the prespecified global/within/non-significant fallback, with actual inhouse support, unit coverage, structure and composition diversity. Leave-one-dataset-out filters section effects then reaggregates shared units. k6/k18 and technical sensitivities are descriptive and do not replace primary P. Ordinary selected-effect displays have no winner-adjusted uncertainty claim.

The live serial controller was launched with bounded recovery supervision. PIDs are recorded in `audit/supervisor.json` and the live `RESUME.json` snapshot. The supervisor waits for the original controller, makes at most two same-configuration recovery attempts, never bypasses RESOURCE_PAUSED, and runs `final_verify.py` only after main analysis and review completion. `audit/final_validation.json` then independently verifies all required artifact hashes, sum(4K), the 3,564 frozen control relations, SciPy BH/BY equality, member/reference axes and actual allocated disk use. A final validation failure changes the saved delivery state to SEARCH_INCOMPLETE. No chat notification or separate scheduled task was created.
