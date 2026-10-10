from pathlib import Path
import os,json,hashlib,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ['PYTHONDONTWRITEBYTECODE']='1'
import numpy as np,pandas as pd
P=Path('/home/data/t070721/codex_workspace/Bladder Metabolism');V=P/'mp10_k32_reference_identifiability_v5';R=V/'scope_frozen_resume_v1'
P07='CNP0000460|CNP0000460_P07T'
def say(*x):print(datetime.datetime.now().isoformat(),*x,flush=True)
def json_default(x):
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,np.ndarray):return x.tolist()
 return str(x)
def js(p,x):
 p=Path(p);t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(x,indent=2,default=json_default)+'\n');t.replace(p)
def wt(x,p):x.to_csv(p,sep='\t',index=False)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def key(x):return hashlib.sha256(('\n'.join(map(str,x))+'\n').encode()).hexdigest()
def seed(purpose):return int.from_bytes(hashlib.sha256(('123|'+purpose).encode()).digest()[:4],'little')%(2**31-1)
def annotate_MP_sources(g):
 mp=pd.read_csv(V.parent/'mp10_k32_validation_v4/00_audit/top50_signatures.tsv',sep='\t');lookup=mp.groupby('gene').MP.agg(lambda x:set(x)).to_dict();membership=[];matched=[]
 for r in g.itertuples():
  source=set(str(r.source_features).split(';'))|{r.gene};hit=sorted(source&set(lookup));programs=set().union(*(lookup[x] for x in hit)) if hit else set();membership.append(';'.join(sorted(programs)));matched.append(';'.join(hit))
 g=g.copy();g['MP_memberships']=membership;g['MP10_input']=['EpiMP10' in x.split(';') for x in membership];g['MP_input_source_features']=matched;g['MP_annotation_rule']='union over exact original source_features and final canonical ID; no fuzzy gene matching';return g
