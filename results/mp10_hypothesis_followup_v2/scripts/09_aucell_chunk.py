import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,pickle,pandas as pd,numpy as np
from pyscenic.aucell import aucell
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';O=P/'results/mp10_immune_endothelial_mechanisms_v1';inp,idsfile,out=sys.argv[1:];genes=pd.read_csv(O/'tables/scenic_genes_Epi.tsv',sep='\t').gene.tolist();ids=pd.read_csv(idsfile,sep='\t').barcode.tolist();x=np.fromfile(inp,dtype='<f8').reshape((len(genes),len(ids)),order='F').T;df=pd.DataFrame(x,index=ids,columns=genes)
regs=pickle.load(open(B/'objects/SREBF_corrected_regulons.pkl','rb'));orig=[r for r in pickle.load(open(O/'objects/scenic/Epi_regulons.pkl','rb')) if r.transcription_factor=='SREBF1'][0];regs.append(orig.copy(name='SREBF1_original(+)'));M20=set((P/'neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().split())
for mask in ['Full','M20']:
 rs=regs if mask=='Full' else [r.copy(gene2weight={g:w for g,w in r.gene2weight.items() if g not in M20}) for r in regs];xx=df if mask=='Full' else df.loc[:,~df.columns.isin(M20)];a=aucell(xx,rs,num_workers=1,auc_threshold=.05,seed=123);a=a.reindex(ids);a.index.name='barcode';a.to_csv(out.replace('.tsv.gz','_'+mask+'.tsv.gz'),sep='\t',compression='gzip')
print('AUCell',len(ids),'cells',flush=True)
