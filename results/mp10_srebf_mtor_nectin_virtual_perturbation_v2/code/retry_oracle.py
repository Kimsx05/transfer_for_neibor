from common import *
from concurrent.futures import ThreadPoolExecutor
from run_oracle import run
m=rd('cache/metadata.tsv.gz')
todo=[(j,s) for j in ['main','seed456']+['drop_'+x for x in sorted(m.dataset.unique())] for s in sorted(m.state.unique()) if not (Path('cache/oracle')/j/s/'complete.json').is_file()]
for j,s in todo:
 old=Path(f'logs/oracle_{j}_{s}.log')
 if old.exists():old.rename(old.with_suffix('.first_attempt'))
with ThreadPoolExecutor(max_workers=6) as pool:
 list(pool.map(run,todo))

