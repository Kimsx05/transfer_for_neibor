from pathlib import Path
import json, hashlib, time, os, resource, contextlib
import pandas as pd
import numpy as np
D=Path(__file__).resolve().parents[1]
os.chdir(D)
P=D.parents[1]
cfg=json.loads(Path('config/design.json').read_text())
base=P/'results/mp10_immune_endothelial_mechanisms_v1/tables'
v1=P/'results/mp10_hypothesis_followup_v1/tables'
v2=P/'results/mp10_hypothesis_followup_v2/tables'
raw=P/'subtype/v2/reference_training_v2/03_derived_reference/reference_v2_raw_counts_pre_harmonization.h5ad'
top=P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt'
old=P/'mp10_subtype_screen_round1/tables/epithelial_existing_metadata.tsv.gz'
targets=['SREBF1','SREBF2','NECTIN2','NECTIN1','NECTIN3','NECTIN4','PVR','HMGCR','SQLE','LDLR']
def rd(p,**kwargs):return pd.read_csv(p,sep='\t',**kwargs)
def save(x,p):x.to_csv(p,sep='\t',index=False,na_rep='')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(4194304),b''):h.update(block)
 return h.hexdigest()
def disk():
 return sum(p.stat().st_size for base in ['cache','tables','plot_data','figures','qc','logs','config','code'] for p in Path(base).rglob('*') if p.is_file())
@contextlib.contextmanager
def measured(stage):
 t=time.monotonic(); c=time.process_time(); d=disk()
 try: yield
 finally:
  f=Path('logs/stages.tsv')
  with f.open('a') as out:
   if out.tell()==0:out.write('stage\twall_s\tcpu_s\tprocess_highwater_RSS_KiB\toutput_byte_change\n')
   out.write(f'{stage}\t{time.monotonic()-t:.4f}\t{time.process_time()-c:.4f}\t{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}\t{disk()-d}\n')
def programs():
 g=rd(v1/'regulatory_gene_membership.tsv.gz')
 ans={k:set(x.gene) for k,x in g.groupby('signature') if k!='mTORC1'}
 for k,x in rd(v2/'scenic_corrected_regulon_edges.tsv').groupby('TF'):ans['SCENIC_'+k]=set(x.target)
 g=rd(v2/'official_gene_sets.tsv.gz');ans['CHOL_SYNTH']=set(g[g.signature=='CHOL_SYNTH'].gene)
 assert ans['MP10']==set(top.read_text().split())
 return ans

