from base import *
import re,zipfile
P=OUT/'04_niche_comparison';d=pd.read_parquet(P/'anchor_MP10_all_layers.parquet');e=pd.read_csv(OUT/'05_source_validation/source_evidence_all9.tsv',sep='\t');assert e.niche_name.tolist()==[f'niche{i}' for i in range(9)];assert len(d)==129632 and d.spot_uid.is_unique
orig=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet').set_index('spot_uid');assert np.array_equal(d.accepted_label.to_numpy(),orig.loc[d.spot_uid,'accepted_label'].to_numpy());assert d.spot_MP10_frozen.notna().all();v=d.epi_score_valid;assert d.loc[v,['spot_MP10_model_universe','epi_MP10_model_universe']].notna().all().all();assert d.loc[d.n_signature_available_model.lt(20),'epi_MP10_model_universe'].isna().all()
r=OUT/'01_pathology_search/user_review';review=pd.read_csv(r/'ROI_review_to_fill.tsv',sep='\t');assert len(review)==133
for c in ['evidence_figure','raw_HE']:
 for rel in review[c]:assert not str(rel).startswith('/');assert (r/rel).exists(),rel
for link in re.findall(r'(?:href|src)="([^"]+)"',(r/'review_index.html').read_text()):assert (r/link).exists(),link
with zipfile.ZipFile(OUT/'reports/review.zip') as z:
 assert z.testzip() is None
 for link in re.findall(r'(?:href|src)="([^"]+)"',z.read('01_pathology_search/user_review/review_index.html').decode()):assert '01_pathology_search/user_review/'+link in z.namelist(),link
 assert not any('counts_cache/' in n or n.endswith('.h5ad') or n.endswith('model.pt') for n in z.namelist())
paths=[P/'anchor_MP10_all_layers.parquet',P/'niche_MP10_summary_all9.tsv',OUT/'05_source_validation/source_evidence_all9.tsv',OUT/'reports/RUN_SUMMARY.md',OUT/'reports/review.zip',OUT/'03_epi_attribution/model_and_gene_universe_manifest.tsv',OUT/'03_epi_attribution/attribution_technical_validation.tsv',OUT/'02_spot_MP10/frozen_reproduction.tsv']
tab(OUT/'00_audit_and_freeze/delivery_checksums.tsv',[dict(relative_path=str(p.relative_to(OUT)),bytes=p.stat().st_size,sha256=sha(p)) for p in paths]);log('DELIVERY VALIDATED: all nine, frozen labels, score coverage, portable ROI links, ZIP integrity')
