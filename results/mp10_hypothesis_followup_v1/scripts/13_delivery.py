import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,hashlib,shutil,datetime,sys,importlib.metadata
import pandas as pd
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v1';S=P/'mp10_hypothesis_followup_v1/scripts'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def inventory():
 files=[p for p in B.rglob('*') if p.is_file()];v=dict(task='mp10_hypothesis_followup_v1',result_directory=str(B),code_directory=str(S.parent),result_file_bytes=sum(p.stat().st_size for p in files if p.name!='delivery_inventory.json'),code_file_bytes=sum(p.stat().st_size for p in S.rglob('*') if p.is_file()),result_file_count=len(files),raw_expression_copied=False,remote_published=False)
 (B/'config/delivery_inventory.json').write_text(json.dumps(v,indent=2));print(json.dumps(v))
if '--inventory-only' in sys.argv:inventory();sys.exit()
# Task-owned, reproducible literature-download intermediates are not part of a light delivery.
p=B/'logs/literature_source_records.json';records=json.loads(p.read_text());manifest=B/'config/literature_retrieval_manifest.json'
if not manifest.exists():
 manifest.write_text(json.dumps(dict(original_record_sha256=sha(p),retained='bibliographic identifiers, titles, MeSH and source metadata; abstracts not redistributed',removed_task_owned_fulltext_intermediates=['PMC6345596.xml','PMC8168716.xml']),indent=2))
if any('abstractText' in r for r in records):
 for r in records:r.pop('abstractText',None)
 p.write_text(json.dumps(records,indent=2))
for n in ['PMC6345596.xml','PMC8168716.xml']:
 f=B/'logs'/n
 if f.exists():f.unlink()
(B/'scripts').mkdir(exist_ok=True)
for p in S.iterdir():
 if p.is_file() and p.suffix in ['.R','.py','.sh']:shutil.copy2(p,B/'scripts'/p.name)
pd.DataFrame([dict(script=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in sorted(S.iterdir()) if p.is_file()]).to_csv(B/'config/code_manifest.tsv',sep='\t',index=False)
pd.DataFrame(sorted([dict(package=d.metadata['Name'],version=d.version) for d in importlib.metadata.distributions()],key=lambda x:x['package'].lower())).to_csv(B/'config/python_environment.tsv',sep='\t',index=False)
if (B/'logs/thread_audit.tsv').exists():shutil.copy2(B/'logs/thread_audit.tsv',B/'qc/thread_audit_delivery_snapshot.tsv')
rows=[]
for p in sorted(B.rglob('*')):
 if not p.is_file() or p.name in ['output_manifest.tsv','resume_state.json','delivery_inventory.json','thread_audit.tsv'] or p.name.startswith('resume_'):continue
 rows.append(dict(path=str(p.relative_to(B)),bytes=p.stat().st_size,sha256=sha(p)))
pd.DataFrame(rows).to_csv(B/'tables/output_manifest.tsv',sep='\t',index=False);inventory();print('Light delivery snapshot complete:',len(rows),'hashed files; no raw expression included')
