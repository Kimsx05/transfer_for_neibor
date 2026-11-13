"""Fail closed on changed original inputs before resuming a task."""
import pathlib,hashlib,csv,json
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(8*1024*1024),b''):h.update(z)
 return h.hexdigest()
rows=list(csv.DictReader(open(B/'tables/input_manifest.tsv'),delimiter='\t'))
for r in rows:
 p=pathlib.Path(r['path']);assert p.exists(),f'Missing original input: {p}'
 assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],f'Changed input; preserve this task and use a new result directory: {p}'
 print('Verified original',p.name,flush=True)
manifest=B/'config/delivered_cache_hashes.json'
if manifest.exists():
 for rel,h in json.loads(manifest.read_text()).items():
  p=B/rel;assert p.exists() and sha(p)==h,f'Changed frozen cache: {p}'
print('INPUT AND CACHE VERIFICATION COMPLETE')
