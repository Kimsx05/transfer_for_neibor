import pathlib,requests,hashlib,json,datetime
P=pathlib.Path.cwd();R=P/'resources/mp10_regulatory_priors';p=R/'hgnc_complete_set_20261118.txt';url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
if not p.exists():
 session=requests.Session();session.trust_env=False;r=session.get(url,timeout=90);r.raise_for_status();assert r.text.startswith('hgnc_id\t');p.write_bytes(r.content)
(R/'hgnc_complete_set_20261118.metadata.json').write_text(json.dumps(dict(url=url,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2));print(p,p.stat().st_size)
