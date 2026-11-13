import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import sys,pickle,pathlib,pandas as pd
from pyscenic.aucell import aucell
workers=min(2,len(os.sched_getaffinity(0)))
lin,inp,out=sys.argv[1:];B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1')
with open(B/f'objects/scenic/{lin}_regulons.pkl','rb') as f:reg=pickle.load(f)
x=pd.read_csv(inp,sep='\t',index_col=0).T
# Same gene universe/order and fixed tie-shuffle seed in every chunk.
a=aucell(x,reg,num_workers=workers,auc_threshold=.05,seed=123);assert set(a.index)==set(x.index);a=a.reindex(x.index);a.index.name='barcode';a.to_csv(out,sep='\t',compression='gzip');print(lin,len(x),'AUCell scored; workers',workers,flush=True)

if lin=='Epi':
 m20=set(pathlib.Path('neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().splitlines())
 masked=[]
 for r in reg:
  weights={g:w for g,w in r.gene2weight.items() if g not in m20}
  if weights:masked.append(r.copy(gene2weight=weights))
 xm=x.loc[:,~x.columns.isin(m20)]
 am=aucell(xm,masked,num_workers=workers,auc_threshold=.05,seed=123);assert set(am.index)==set(xm.index);am=am.reindex(xm.index);am.index.name='barcode';am.to_csv(out.replace('.tsv.gz','_M20-masked.tsv.gz'),sep='\t',compression='gzip')
 print('Epi M20-masked AUCell scored',len(xm),flush=True)
