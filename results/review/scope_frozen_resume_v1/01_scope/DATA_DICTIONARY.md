# Scope and label fields

`reference_scope_manifest.tsv.gz` retains300033 old-reference cells for audit;281851 are included in this run and18182 are explicitly excluded. Inclusion is the union of all old non-Epithelial members and123699 frozen epithelial-analysis members.

- `old_reference_label`: final label of the historical30-class reference.
- `in_frozen_scope`: exact membership in the frozen epithelial discovery set. Old non-epithelial cells need not belong to that set to remain eligible.
- `frozen_cluster`: the retained v4 cluster annotation, with missing values for outside-frozen epithelial members. Frozen K32 is exactly the3282 members of FULL_Z_r1.00 C18.
- `new_reference_label`: the ONLY model label column. It contains Epi_K32,Epi_other and the unchanged other29 old subtype names. Excluded epithelial cells have no new reference label.
- `included_in_this_run`: eligibility in the derived candidate reference universe, before TRAIN/DEV/TEST separation. It does not mean a DEV/TEST cell trained a reference.
- `source_group` and `split`: unchanged conservative source bundle and TRAIN/DEV/TEST assignment. Groups are not independently verified patients.
- `exclusion_reason`: outside_frozen_discovery_scope for the18182 excluded epithelial cells. The18120 earlier-cleaning reasons remain unverified;62 smaller-cluster removals retain their documented stage. Neither category is called QC-fail or K32-negative.
- `new_subtype`, `new_major` and `label_status`: retained legacy v5 audit provenance, not model inputs. The legacy audit name Epi_rest is replaced by Epi_other only in the authoritative new_reference_label column. Legacy new_major preserves original source major, including Cycling.

`new_subtype_major_mapping.tsv` is the label-level reference major map; both new epithelial labels inherit Epithelial. `new_subtype_original_major_provenance.tsv` separately preserves original source annotations. All4344 frozen Cycling-origin epithelial-reference members remain included (74 Epi_K32 and4270 Epi_other); no major==Epithelial filter was used.

`heldout_cell_simulation_use_manifest.tsv.gz` records whether each eligible DEV/TEST cell was actually selected by the frozen simulation draws. DEV uses593/634 eligible K32 and TEST253/253. The41 unsampled DEV K32 remain eligible members, all from HRA003620_BC-04; they are not excluded or relabeled. Counts per original cell, source and subtype are in actual_heldout_cell_simulation_coverage.tsv. Pseudospot replicates are technical simulations, not independent patients.
