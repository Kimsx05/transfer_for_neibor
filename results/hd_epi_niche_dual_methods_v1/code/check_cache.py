from pathlib import Path
import json,hashlib,sys
R=Path(__file__).resolve().parents[1];A=R.parent;paths=[]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 base=A/('hd_bin8_comparison_v1' if s in ['Cervilla','D1_Johnson'] else 'hd_bin8_banksy_extension_v1')/s
 paths.extend([base/'bin_metadata.tsv.gz',base/'analysis_metadata.tsv.gz',A/'hd_bin8_banksy_extension_v1'/s/'comparison_metadata.tsv.gz',A/'hd_bin8_banksy_extension_v1'/s/'supplement_gene_coverage.tsv'])
state={'small_definition_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'config.json',R/'type_mapping.json',R/'code/prepare.py',R/'code/model.py',R/'code/zero_distance_review.py']},'input_metadata':{str(p.relative_to(A)):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in paths}}
f=R/'cache_contract.json'
if '--record' in sys.argv:
 assert not f.exists(),'Existing contract will not be silently replaced'
 f.write_text(json.dumps(state,indent=2));print('CACHE_CONTRACT_RECORDED')
else:
 assert json.loads(f.read_text())==state,'Input/definition changed: do not reuse old model cache. Reconcile explicitly in a new version.'
 print('CACHE_CONTRACT_OK')
