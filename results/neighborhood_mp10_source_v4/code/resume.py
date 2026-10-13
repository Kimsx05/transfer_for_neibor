from base import *
import subprocess,argparse
p=argparse.ArgumentParser(description='v4 only. Old labels/models read-only; completed stages are skipped. No clustering/training.');p.add_argument('--stage',choices=['all','scores','compare','source','reports','verify'],default='all');a=p.parse_args()
h=pd.read_csv(OUT/'00_audit_and_freeze/niche_input_hashes.tsv',sep='\t');assert all(sha(r.path)==r.sha256 for r in h.itertuples()),'Frozen input changed: stop for reconciliation'
def run(name):subprocess.run(['/usr/bin/python3.12','-u',str(OUT/'code'/name)],cwd=ROOT,check=True)
def scores_ok():
 p=OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv';q=OUT/'04_niche_comparison/anchor_MP10_all_layers.parquet'
 return p.exists() and q.exists() and pd.read_csv(p,sep='\t').status.isin(['COMPLETE_FULL20','COMPLETE_PARTIAL_SIGNATURE','NO_COMMON_ANCHORS']).all()
if a.stage=='verify':
 c=OUT/'00_audit_and_freeze/delivery_checksums.tsv';assert c.exists(),'Final checksum manifest absent'
 for r in pd.read_csv(c,sep='\t').itertuples():assert sha(OUT/r.relative_path)==r.sha256,r.relative_path
 assert scores_ok();log('Verified frozen inputs and completed delivery checksums. No computation launched.');exit()
if a.stage in ['all','scores']:
 if not (OUT/'01_pathology_search/spot_pathology_overlay.tsv').exists():run('01_pathology.py')
 if not scores_ok():run('02_scores.py')
 else:log('scores complete: skip')
if a.stage in ['all','compare']:
 assert scores_ok()
 if a.stage=='compare' or not (OUT/'04_niche_comparison/technical_mixture_adjustment.tsv').exists() or not (OUT/'06_figures/signature_count_shares.pdf').exists():run('03_compare.py')
 else:log('comparison complete: skip')
if a.stage in ['all','source']:
 assert scores_ok()
 if a.stage=='source' or not (OUT/'05_source_validation/source_validation_manifest.json').exists():run('04_source.py')
 else:log('source controls complete: skip')
if a.stage in ['all','reports']:
 if not (OUT/'00_audit_and_freeze/model_checkpoint_hashes.tsv').exists():run('05_additional_audits.py')
 if a.stage=='reports' or not (OUT/'reports/review.zip').exists():run('06_finalize.py')
 else:log('reports and package complete: skip')
run('07_verify_delivery.py')
