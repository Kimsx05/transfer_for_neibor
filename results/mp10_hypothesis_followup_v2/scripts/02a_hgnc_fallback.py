import requests,pathlib,json,hashlib,datetime
P=pathlib.Path.cwd();R=P/'resources/mp10_regulatory_priors';p=R/'hgnc_complete_set_20261118.txt';attempts=[]
urls=['https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt','https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt']
for url in urls:
 for proxy in [True,False]:
  if p.exists():break
  try:
   s=requests.Session();s.trust_env=proxy;r=s.get(url,timeout=(12,35));r.raise_for_status();assert r.text.startswith('hgnc_id\t');p.write_bytes(r.content);attempts.append(dict(url=url,trust_env=proxy,status='completed'));break
  except Exception as e:attempts.append(dict(url=url,trust_env=proxy,status='failed',error=str(e)[:300]))
 if p.exists():break
(R/'hgnc_complete_set_20261118.metadata.json').write_text(json.dumps(dict(attempts=attempts,retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),path=str(p),bytes=p.stat().st_size if p.exists() else 0,sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None),indent=2));print(json.dumps(attempts),flush=True)
