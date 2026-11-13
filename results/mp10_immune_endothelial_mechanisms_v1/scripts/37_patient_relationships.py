"""Keep original patient-ID provenance and separately record the user-confirmed repeat relationship."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
import pathlib,json,csv,collections
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');T=B/'tables'
with (T/'specimen_exposure_eligibility.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
groups=collections.defaultdict(list)
for r in rows: groups[r['patient_uid_recorded']].append(r['sample_uid'])
out=[]
for r in rows:
 confirmed=r['sample_uid'] in ['InhouseData::T4463209a','InhouseData::T4463209b']
 out.append(dict(sample_uid=r['sample_uid'],dataset=r['dataset'],patient_uid_recorded=r['patient_uid_recorded'],source_patient_id_status=r['patient_id_status'],source_patient_id_source=r['patient_id_source'],n_specimens_with_same_recorded_id=len(groups[r['patient_uid_recorded']]),related_specimens=';'.join(groups[r['patient_uid_recorded']]),repeat_relationship_status='confirmed_same_patient' if confirmed else 'no_additional_confirmed_relationship',relationship_evidence='explicit_user_instruction' if confirmed else 'source_metadata_only',independence_from_other_recorded_ids='not_fully_verified',main_specimen=r['main_patient'],replacement_specimen=r['replacement_patient'],interpretation='Original provisional ID field is retained; a/b repeat relationship is independently confirmed' if confirmed else 'Do not merge short specimen names across datasets or assume verified independent patients'))
with (T/'patient_relationship_audit.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
assert len(rows)==40 and len(groups)==39
assert [v for v in groups.values() if len(v)>1]==[['InhouseData::T4463209a','InhouseData::T4463209b']]
import gzip
with gzip.open(T/'model_membership.tsv.gz','rt') as f: members=list(csv.DictReader(f,delimiter='\t'))
models=collections.defaultdict(set)
for r in members: models[r['model_id']].add(r['sample_uid'])
checks={k:not {'InhouseData::T4463209a','InhouseData::T4463209b'}.issubset(v) for k,v in models.items()}
assert all(checks.values())
(B/'qc/patient_relationships.json').write_text(json.dumps(dict(original_specimens=40,recorded_units=39,confirmed_repeat_relationship='T4463209a/b',verified_independent_patients='not established for all records',no_joint_a_b_in_any_model=all(checks.values()),n_models_checked=len(checks)),indent=2))
print('PATIENT RELATIONSHIP AUDIT COMPLETE')
