from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];P=R.parents[2];v2=P/'neighborhood_ring_comparison_v2';v3=P/'neighborhood_ring_harmony_v3'
types=json.loads((v2/'03_features/feature_columns.json').read_text())['q1'];mapping={t:('Epithelial' if t=='Epithelial' else 'Endothelial' if t.endswith('EC') else 'Mesenchymal' if 'CAF' in t or t in ['Pericyte','SMC'] else 'Immune') for t in types}
pd.DataFrame([{'ordinary_type':t,'common_coarse':mapping[t]} for t in types]).to_csv(R/'tables/ordinary_visium_type_correspondence.tsv',sep='\t',index=False)
m=pd.read_parquet(v2/'03_features/common_anchors.parquet',columns=['spot_uid']);f=np.load(v2/'03_features/features.npz');lab=pd.read_parquet(v3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet');assert m.spot_uid.is_unique and lab.spot_uid.is_unique;lab=lab.set_index('spot_uid').reindex(m.spot_uid);assert lab.accepted_label.notna().all();q=(f['q1']+f['q2'])/2;assert len(q)==len(m) and q.shape[1]==len(types)
rows=[]
for label in sorted(lab.accepted_label.unique()):
 if label<0:continue
 z=q[lab.accepted_label.to_numpy()==label];rows.append(dict(ordinary_niche=int(label),n=len(z),**{co:float(z[:,[i for i,t in enumerate(types) if mapping[t]==co]].sum(1).mean()) for co in ['Epithelial','Mesenchymal','Endothelial','Immune']},definition='Descriptive equal ring mean (R1+R2)/2; frozen accepted R12 labels; c2l q05 per-spot composition, not HD votes'))
pd.DataFrame(rows).to_csv(R/'tables/ordinary_visium_frozen_patterns.tsv',sep='\t',index=False)
print('ORDINARY_FROZEN_READ',len(q),'anchors; no model or MP10 loaded')
