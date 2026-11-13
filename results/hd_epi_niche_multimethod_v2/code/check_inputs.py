from common import *
m=pd.read_csv(R/'input_manifest.tsv',sep='\t')
for _,v in m.iterrows():p=Path(v.path);assert p.stat().st_size==v.bytes and p.stat().st_mtime_ns==v.mtime_ns,str(p)
a=anchors();T=train_index(a);assert len(T)==80000 and a.key.is_unique
print('FROZEN_INPUTS_AND_T_CONFIRMED')

import hashlib
files=['config.yaml','type_mapping.json','code/common.py','code/graph_core.py','code/fit_composition.py','code/rings.py','code/expression.py','code/cellcharter_models.py','code/mender_ring_sensitivity.py','tables/T_v2.tsv']
sig={f:hashlib.sha256((R/f).read_bytes()).hexdigest() for f in files}
p=R/'cache_contract.json'
if p.exists():assert json.loads(p.read_text())==sig,'Definitions changed: explicitly reconcile existing task checkpoints before reuse'
else:p.write_text(json.dumps(sig,indent=2))
print('DEFINITION_CONTRACT_OK')

cm=pd.read_csv(R/'tables/counts_input_audit.tsv',sep='\t')
for _,v in cm.iterrows():
 p=Path(v.counts_path);assert p.stat().st_size==v.bytes and p.stat().st_mtime_ns==v.mtime_ns,str(p)
print('RAW_COUNT_PATH_IDENTITIES_CONFIRMED')
