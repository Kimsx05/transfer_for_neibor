from common import *
# Bound the requested order statistic to the number of scores after CellOracle's own filter.
# No expression/motif/effect/perturbation gate is relaxed.
s=Path('cache/co310/lib/python3.10/site-packages/celloracle/trajectory/modified_VelocytoLoom_class.py')
before='nth_score = np.sort(score)[::-1][N]';after='nth_score = np.sort(score)[::-1][min(N, len(score) - 1)]'
src=s.read_text()
assert src.count(before)+src.count(after)==1
orig=src.replace(after,before);patched=src.replace(before,after)
if patched!=src:
 new=s.with_suffix('.tmp');new.write_text(patched);os.replace(new,s)
Path('qc/oracle_qc_compat.json').write_text(json.dumps(dict(package='celloracle==0.18.0',original_sha256=hashlib.sha256(orig.encode()).hexdigest(),patched_sha256=sha(s),original=before,replacement=after,reason='empirical CV score vector has fewer entries after original expression filters; cap requested order statistic; no change when original call succeeds'),indent=2))

