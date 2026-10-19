from core import *
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon,rankdata
extract(V8/'code/common.py',['signed','nearest','induced_components','bh'],globals())
extract(V8/'code/relations.py',['stratify','metric_record'],globals())
KEY=['layer','k','track','target_kind','target_id','metric','analysis_domain']
def family_id(branch,track,kind,metric,subset):
 end='O_D' if kind=='NICHE' and metric in ['O','D'] else metric
 return '|'.join([branch,track,kind,end,subset])
def summarize_spotted(s,K,branch):
 units=[];tests=[]
 for scope in ['ALL_DATASETS']+sorted(s.dataset_id.unique()):
  sub=s if scope=='ALL_DATASETS' else s[s.dataset_id.eq(scope)]
  for key,g in sub.groupby(KEY,sort=False):
   base=dict(zip(KEY,key));eligible=g[g.eligible & g.effect.notna() & g.analysis_unit_uid.fillna('').ne('')]
   ur=[]
   for uid,a in eligible.groupby('analysis_unit_uid',sort=False):
    r=dict(scope=scope,**base,analysis_unit_uid=uid,effect=a.effect.mean(),n_sections=len(a),known_patient=bool(a.known_patient.all()),datasets='|'.join(sorted(a.dataset_id.unique())),sections='|'.join(sorted(a.section_id)),stratified_effect=a.stratified_effect.mean(),component_equal_effect=a.component_equal_effect.mean());ur.append(r);units.append(r)
   for subset in ['ALL_POLICY_UNITS','VERIFIED_PATIENTS_ONLY']:
    u=[r for r in ur if subset=='ALL_POLICY_UNITS' or r['known_patient']];v=signed([r['effect'] for r in u]);f=family_id(branch,base['track'],base['target_kind'],base['metric'],subset)
    tests.append(dict(branch=branch,scope=scope,subset=subset,**base,family_id=f,planned_n=(4*K if base['metric'] in ['O','D'] else 2*K) if base['target_kind']=='NICHE' else 58,n_sections=sum(r['n_sections'] for r in u),n_sections_descriptive=int(g.effect.notna().sum()),n_sections_target_present=int(g.n_target_full.gt(0).sum()),n_sections_hot=int(g.n_source_hot.gt(0).sum()),**v))
 u=pd.DataFrame(units);r=pd.DataFrame(tests);r['q_within_BH']=np.nan
 for key,ids in r.groupby(['scope','family_id','k','analysis_domain']).groups.items():
  expected=int(r.loc[ids,'planned_n'].iloc[0]);assert len(ids)==expected,(key,len(ids),expected)
  r.loc[ids,'q_within_BH']=correction(r.loc[ids,'p'])
 return u,r
