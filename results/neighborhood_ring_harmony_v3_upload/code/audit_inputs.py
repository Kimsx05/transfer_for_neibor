from base import *
from scipy.spatial.distance import pdist
cov=pd.read_csv(OLD/'00_preflight/coverage63.tsv',sep='\t');rows=[]
for r in cov.itertuples():
 for layer in ['q05','mean']:
  path=Path(getattr(r,layer+'_path'));actual=sha(path);expected=getattr(r,layer+'_sha256');rows.append(dict(section_id=r.section_id,layer=layer,path=str(path),expected_sha256=expected,actual_sha256=actual,match=actual==expected))
assert all(r['match'] for r in rows);tab(OUT/'00_preflight/accepted_input_reverification.tsv',rows)
# Verify existing masks and features by independent algebra on stored q rather than reconstruct accepted estimates.
m,f,ref=inputs();checks=[]
for key,value in [('X1',np.sqrt(f['q1'])),('X2',np.sqrt(f['q2'])),('X12',np.concatenate([np.sqrt(f['q1']),np.sqrt(f['q2'])],1)/np.sqrt(2))]:
 err=float(np.max(np.abs(f[key]-value)));assert err<1e-6;checks.append(dict(check=key,max_error=err,status='PASS'))
assert np.array_equal(m.iloc[ref].spot_uid.to_numpy(),pd.read_csv(OLD/'03_features/reference_anchor_ids.tsv',sep='\t').spot_uid.to_numpy());checks.append(dict(check='reference UID order identical',max_error=0,status='PASS'))
for ring in ['R1','R2','R12']:
 p=joblib.load(OUT/'01_features'/f'{ring}_PCA.joblib');ev=p.explained_variance_ratio_;assert ev.sum()>=.95 and ev[:-1].sum()<.95;checks.append(dict(check=ring+' minimum PCA95 dimensions',max_error=0,status='PASS'))
tab(OUT/'00_preflight/feature_reuse_validation.tsv',checks);js(OUT/'00_preflight/composition_strata_definition.json',dict(definition='argmax((q1+q2)/2) in frozen30 subtype order; deterministic ties first index',purpose='descriptive conditional mixing; not pathology ground truth',frozen_input_hash=sha(OLD/'03_features/features.npz')));log('Accepted input hashes and feature algebra PASS')
