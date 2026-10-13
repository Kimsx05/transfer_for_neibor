from common import *
import argparse,subprocess
p=argparse.ArgumentParser();p.add_argument('--verify-only',action='store_true');p.add_argument('--refresh-reports',action='store_true');a=p.parse_args()
for r in pd.read_csv(OUT/'00_audit_and_spec/input_checksums.tsv',sep='\t').itertuples():assert sha(r.path)==r.sha256,r.path
assert sha(OUT/'00_audit_and_spec/frozen_analysis_spec.yaml')==json.loads((OUT/'00_audit_and_spec/specification_hash.json').read_text())['sha256']
def run(n):subprocess.run(['/usr/bin/python3.12','-u',str(OUT/'code'/n)],cwd=ROOT,check=True)
if a.verify_only:
 for r in pd.read_csv(OUT/'00_audit_and_spec/final_checksums.tsv',sep='\t').itertuples():assert sha(OUT/r.relative_path)==r.sha256,r.relative_path
 log('Verified specification, immutable inputs and final artifacts. No analysis/fitting launched.');exit()
stages=[('01_score_fidelity/cell_level_fidelity.parquet','01_fidelity.py'),('02_gene_evidence/gene_by_section_niche.tsv.gz','02_spatial_genes.py'),('03_attribution_calibration/known_mixture_scores.parquet','03_panel.py'),('03_attribution_calibration/frozen_replay_scores.parquet','05_replay.py'),('04_fixed_niche_contrasts/paired_dataset_effects.tsv','06_contrasts.py')]
for artifact,script in stages:
 if not (OUT/artifact).exists():run(script)
 else:log('skip completed '+script)
P=OUT/'03_attribution_calibration/synthetic_refits/model1'
if not (P/'fit_status.json').exists():
 if (P/'model/model.pt').exists():run('04b_resume_synthetic.py')
 else:raise RuntimeError('No recoverable synthetic checkpoint; two fresh attempts limit must not be exceeded')
if not (OUT/'03_attribution_calibration/end_to_end_scores.parquet').exists() or a.refresh_reports:run('08_calibration_summary.py')
if not (OUT/'reports/review.zip').exists() or a.refresh_reports:run('10_finish.py')
run('11_verify.py')
