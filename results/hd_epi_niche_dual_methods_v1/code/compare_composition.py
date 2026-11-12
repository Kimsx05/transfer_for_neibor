from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];rows=[]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 for r in [50,200]:
  d=pd.read_parquet(R/s/f'r{r}.parquet');out=d[['bin_id']].copy();coarse={}
  for ref in ['blca','hpca']:
   dic={'Epithelial':['Epithelial'],'Mesenchymal':['Fibroblast','Mural'],'Endothelial':['Endothelial'],'Immune':['T/NK','Myeloid','B','Plasma','Mast']};v=[]
   for k,ts in dic.items():
    cols=[ref+'_P_'+t for t in ts if ref+'_P_'+t in d];x=d[cols].sum(axis=1,min_count=1);out[ref+'_'+k]=x;v.append(x)
   coarse[ref]=np.array(v).T
  out['coarse_L1']=np.abs(coarse['blca']-coarse['hpca']).sum(1);out['known_denominator_difference']=d.blca_n_known-d.hpca_n_known
  out.to_parquet(R/s/f'coarse_reference_r{r}.parquet',index=False)
  for band,mask in [('all',np.ones(len(d),bool)),('<20',d.nFeature_RNA<20),('20-99',d.nFeature_RNA.between(20,99)),('100-499',d.nFeature_RNA.between(100,499)),('>=500',d.nFeature_RNA>=500)]:
   q=out[mask];rows.append(dict(section=s,radius=r,center_information=band,n=len(q),coarse_L1_mean=q.coarse_L1.mean(),coarse_L1_q90=q.coarse_L1.quantile(.9),absolute_known_count_difference=q.known_denominator_difference.abs().mean()))
pd.DataFrame(rows).to_csv(R/'tables/reference_coarse_composition.tsv',sep='\t',index=False)
