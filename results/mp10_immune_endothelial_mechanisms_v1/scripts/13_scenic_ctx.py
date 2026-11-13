"""Motif pruning of a frozen lineage GRN; only this task's generated files are written."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import sys,pathlib,hashlib,json,pickle,subprocess
import pandas as pd
import pyarrow as pa
pa.set_cpu_count(1);pa.set_io_thread_count(1)
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df,df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase
P=pathlib.Path.cwd();B=P/'results/mp10_immune_endothelial_mechanisms_v1';R=P/'resources/mp10_regulatory_priors';lin=sys.argv[1];O=B/'objects/scenic';T=B/'tables';seed=123;workers=int(os.environ.get('SCENIC_CTX_WORKERS','1'));assert workers in (1,2)
exprfile=O/f'{lin}_training.tsv.gz';adjfile=O/f'{lin}_seed123_adjacencies.tsv';dbfile=R/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather';motif=R/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for z in iter(lambda:f.read(1024*1024),b''):h.update(z)
 return h.hexdigest()
assert dbfile.exists() and motif.exists(),'Complete matching motif references required'
assert adjfile.exists() and adjfile.stat().st_size>0,'GRNBoost2 unfinished'
params=dict(lineage=lin,seed=123,rho_threshold=.03,rho_mask_dropouts=False,min_genes=20,rank_threshold=5000,auc_threshold=.05,nes_threshold=3.0,num_workers=workers,mode='custom_multiprocessing',module_thresholds=[.75,.90],top_n_targets=[50],top_n_regulators=[5,10,50],rho_dichotomize=True,keep_only_activating=True,motif_similarity_fdr=.001,orthologuous_identity_threshold=0.0,training_sha256=sha(exprfile),adjacency_sha256=sha(adjfile),db_sha256=sha(dbfile),motif_sha256=sha(motif),method='pyscenic0.12.1_GRNBoost2_ctx')
key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest();done=O/f'{lin}_regulons.metadata.json'
if done.exists() and (O/f'{lin}_regulons.pkl').exists():
 old=json.loads(done.read_text());ignored={'num_workers','input_hash','n_modules','n_regulons'}
 if {k:v for k,v in old.items() if k not in ignored}=={k:v for k,v in params.items() if k not in ignored}:print('Verified ctx scientific-input cache',lin);sys.exit()

ex=pd.read_csv(exprfile,sep='\t',index_col=0).T;adj=pd.read_csv(adjfile,sep='\t');assert set(adj.columns)>={'TF','target','importance'}
modules=list(modules_from_adjacencies(adj,ex,min_genes=20,rho_threshold=.03,rho_mask_dropouts=False));print(lin,'modules',len(modules),flush=True)
db=FeatherRankingDatabase(str(dbfile),name='hg38_refseq-r80_mc9nr_10kb');df=prune2df([db],modules,str(motif),rank_threshold=5000,auc_threshold=.05,nes_threshold=3.0,num_workers=workers,client_or_address='custom_multiprocessing');df.to_csv(O/f'{lin}_motif_enrichment.csv');reg=df2regulons(df);print(lin,'regulons',len(reg),flush=True)
with open(O/f'{lin}_regulons.pkl','wb') as f:pickle.dump(reg,f)
rows=[];coverage=[]
for r in reg:
 for g,w in r.gene2weight.items():rows.append(dict(lineage=lin,regulon=r.name,TF=r.transcription_factor,target=g,weight=w))
 coverage.append(dict(lineage=lin,regulon=r.name,TF=r.transcription_factor,n_targets=len(r.genes),n_covered_targets=len(set(r.genes)&set(ex.columns)),small_regulon=len(r.genes)<20))
pd.DataFrame(rows).to_csv(T/f'scenic_{lin}_regulon_edges.tsv.gz',sep='\t',index=False);pd.DataFrame(coverage).to_csv(T/f'scenic_{lin}_regulon_coverage.tsv',sep='\t',index=False)
params.update(input_hash=key,n_modules=len(modules),n_regulons=len(reg));done.write_text(json.dumps(params,indent=2));print('CTX COMPLETE',lin,flush=True)
