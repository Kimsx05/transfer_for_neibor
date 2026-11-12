from pathlib import Path
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1]
for s in ['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']:
 d=pd.read_parquet(R/s/'nodes.parquet');a=d[d.epi_anchor];tree=cKDTree(d[['x_um','y_um']]);out=a[['bin_id']].copy()
 for r in [50,200]:
  out[f'r{r}_supplied_all']=tree.query_ball_point(a[['x_um','y_um']],r+1e-7,return_length=True,workers=4)-1
  if not s.startswith('MPBC'):
   t=d[d.in_tissue_original==1];out[f'r{r}_original_in_tissue']=cKDTree(t[['x_um','y_um']]).query_ball_point(a[['x_um','y_um']],r+1e-7,return_length=True,workers=4)-1
  else:out[f'r{r}_original_in_tissue']=np.nan
 out.to_parquet(R/s/'geometric_coverage.parquet',index=False)
