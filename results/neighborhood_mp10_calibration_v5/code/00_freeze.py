from common import *
import shutil,subprocess

def main(out):
 p=out/'00_audit_and_spec';attach=Path('/home/data/t070721/.codex/attachments/01e40887-8d63-4a35-83e7-3d62ae944484')
 for q in attach.glob('*.md'):shutil.copy2(q,safe(out,p/q.name))
 paths=[V4/'04_niche_comparison/anchor_MP10_all_layers.parquet',V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet',OLD/'00_preflight/MP10_top20.tsv',ROOT/'ST/mapping/07_EpiMP_basis_sample_balanced.csv',ROOT/'ST/mapping/08_cell_EpiMP_scores.qs']
 paths+=list((V4/'05_source_validation/counts_cache').glob('*/axis_manifest.json'));paths+=list((V4/'05_source_validation/counts_cache').glob('*/*.npz'))
 frozen=p/'input_checksums.tsv'
 if frozen.exists():
  for r in pd.read_csv(frozen,sep='\t').itertuples():assert sha(r.path)==r.sha256
 else:write_tsv(out,frozen,[dict(path=str(q),sha256=sha(q),bytes=q.stat().st_size) for q in paths])
 spec=dict(seed=123,status='RESULT_INFORMED_EXPLORATORY_CALIBRATION',fixed_pairs=PAIRS,min_n_each=20,min_unit_cells_fidelity=30,min_CI_units=5,bootstrap=2000,matching={'epi_caliper':.05,'logUMI_caliper_sd':.25,'replacement':False},near_zero='1/original_AUC_denominator; additionally paired Monte Carlo noise from fixed-state repeats',reference_background_selection='first lexicographic full20 section CanDis_2026_Sample1; same30 mapping',new_synthetic_models_max=2,new_synthetic_models_planned=1,synthetic_bags_max=2000,extra_panel_max=3000,optimizer_checkpoints=[5000,10000,20000],convergence='finite loss; abs difference last500 vs previous500 median / abs previous median <0.01; no truth score selection',threads=4,synthetic_fit_device='CPU while existing GPU jobs active',official_AUCell_seeds=[123,124,125],primary_score='unchanged v4 average-ties; full model genes; Top20',source_interpretation='replay distinct from known-mixture distinct from new synthetic fit')
 if not (p/'frozen_analysis_spec.yaml').exists():write_json(out,p/'frozen_analysis_spec.yaml',spec)
 d=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet');write_tsv(out,p/'identity_and_independence.tsv',d[['section_id','dataset_id','patient_uid','analysis_unit_uid','analysis_unit_type','source_patient_identity_known']].drop_duplicates());write_tsv(out,p/'resolved_inputs.tsv',[dict(role='frozen_asset',path=str(q),exists=q.exists()) for q in paths[:5]]);safe(out,p/'execution_scope.md').write_text('v4/v3/v2 and real models read-only. All work is outcome-informed exploration. One CPU synthetic refit planned, maximum two. No pathology restart or new niche. All output writes guarded beneath v5.\n');log('specification frozen')
if __name__=='__main__':main(OUT)
