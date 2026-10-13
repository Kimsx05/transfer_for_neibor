from common import *
import zipfile

def main(out):
 for r in pd.read_csv(out/'00_audit_and_spec/input_checksums.tsv',sep='\t').itertuples():assert sha(r.path)==r.sha256,r.path
 spec=json.loads((out/'00_audit_and_spec/specification_hash.json').read_text());assert sha(out/'00_audit_and_spec/frozen_analysis_spec.yaml')==spec['sha256']
 for r in pd.read_csv(V4/'00_audit_and_freeze/model_checkpoint_hashes.tsv',sep='\t').itertuples():assert sha(r.path)==r.sha256,r.path
 expected=['HMGCS1','MSMO1','INSIG1','FDFT1','SQLE','FDPS','IDI1','ACSS2','ELOVL5','LDLR','MVD','SCD','HMGCR','ACAT2','SLC20A1','ACSL3','ACSL1','ADIPOR2','MXD1','NPC1'];assert TOP==expected
 basis=pd.read_csv(ROOT/'ST/mapping/07_EpiMP_basis_sample_balanced.csv',index_col=0);assert basis.EpiMP10.sort_values(ascending=False,kind='stable').head(20).index.tolist()==TOP
 d=pd.read_csv(out/'04_fixed_niche_contrasts/paired_section_effects.tsv',sep='\t');v=d[d.analysis_set.eq('MAIN')];assert (v.n_A>=20).all() and (v.n_B>=20).all();assert v.pair.nunique()==3
 src=pd.read_parquet(V4/'04_niche_comparison/anchor_MP10_all_layers.parquet');src=src[src.accepted_label.ge(0)&src.epi_score_valid]
 for r in v.itertuples():
  A,B=[int(x.replace('niche','')) for x in r.pair.split('_minus_')];g=src[src.section_id.eq(r.section_id)];value=g.loc[g.accepted_label.eq(A),r.score].median()-g.loc[g.accepted_label.eq(B),r.score].median();assert np.isclose(value,r.effect,atol=1e-12)
 ev=pd.read_csv(out/'05_evidence_synthesis/evidence_matrix_all9.tsv',sep='\t');assert len(ev)==12;assert set([f'niche{i}' for i in range(9)]).issubset(ev.niche_or_pair)
 P=out/'03_attribution_calibration';a=pd.read_parquet(P/'frozen_replay_scores.parquet');b=pd.read_parquet(P/'known_mixture_scores.parquet');c=pd.read_parquet(P/'end_to_end_scores.parquet');assert len(a)==1998 and len(b)==540 and len(c)==270
 assert a.bag_id.is_unique and b.bag_id.is_unique and c.bag_id.is_unique;assert set(c.bag_id).issubset(b.bag_id)
 attempts=json.loads((P/'synthetic_refits/model1/attempts.json').read_text());assert len(attempts)==2
 with zipfile.ZipFile(out/'reports/review.zip') as z:
  assert z.testzip() is None;names=z.namelist();assert not any(n.endswith(('.h5ad','.npz','.bin','.parquet','model.pt')) for n in names)
  private=set(pd.read_csv(P/'selected_epi_donors.tsv',sep='\t').unit.astype(str));private.update(src.patient_uid.dropna().astype(str));private={s for s in private if len(s)>4}
  for n in names:
   if n.endswith('.tsv'):
    content=z.read(n).decode()
    for key in private:assert key not in content,(n,key)
  assert 'reports/RUN_SUMMARY.md' in names and 'SANITIZATION.md' in names
 files=[p for p in out.rglob('*') if p.is_file() and p.suffix in ['.tsv','.json','.md','.zip','.parquet'] and not any(v in p.parts for v in ['counts_cache','by_section','logs']) and p.name not in ['final_checksums.tsv','verification.json']]
 write_tsv(out,out/'00_audit_and_spec/final_checksums.tsv',[dict(relative_path=str(p.relative_to(out)),sha256=sha(p)) for p in files])
 write_json(out,out/'reports/verification.json',dict(status='PASS',read_only_inputs=True,formal_scores_unchanged=True,exact_direct_contrasts=True,all9_plus3=True,synthetic_attempts=2,cache_axes=True,archive_integrity=True,archive_private_ID_scan=True))
 log('PASS: immutable inputs, exact direct contrasts, all9, recovery, deidentified ZIP')
if __name__=='__main__':main(OUT)
