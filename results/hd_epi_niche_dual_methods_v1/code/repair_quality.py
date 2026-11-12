from pathlib import Path
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 d=pd.read_parquet(R/s/'nodes.parquet');a=d[d.epi_anchor].set_index('bin_id',drop=False)
 for r in [50,200]:
  f=R/s/f'r{r}.parquet';p=pd.read_parquet(f).set_index('bin_id',drop=False)
  for comp,g in a.groupby('component'):
   pool=d[(d.component==comp)&d.analysis_qc];n=p.loc[g.index,'n_analysis_neighbors'].to_numpy()
   for field,mask,selfmask in [('low_information_fraction',pool.nFeature_RNA<20,g.nFeature_RNA<20),('reference_conflict_fraction',pool.reference_conflict,g.reference_conflict)]:
    q=pool[mask];v=cKDTree(q[['x_um','y_um']]).query_ball_point(g[['x_um','y_um']],r+1e-7,return_length=True,workers=4) if len(q) else np.zeros(len(g));p.loc[g.index,field]=np.divide(v-selfmask.to_numpy(),n,out=np.full(len(g),np.nan),where=n>0)
  assert p[['low_information_fraction','reference_conflict_fraction']].notna().all().all();assert p[['low_information_fraction','reference_conflict_fraction']].ge(0).all().all();assert p[['low_information_fraction','reference_conflict_fraction']].le(1).all().all();p.reset_index(drop=True).to_parquet(f,index=False)
 print('QUALITY_FIXED',s,flush=True)
