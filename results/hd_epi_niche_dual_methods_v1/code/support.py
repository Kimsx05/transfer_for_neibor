from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parents[1];E=R.parent/'hd_bin8_banksy_extension_v1';sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2']
for s in sections:
 o=R/s
 if (o/'neighbor_marker_r200.parquet').exists():continue
 d=pd.read_parquet(o/'nodes.parquet');base=R.parent/'hd_bin8_comparison_v1' if s in sections[:2] else E
 m=pd.read_csv(base/s/'analysis_metadata.tsv.gz',sep='\t').set_index('bin_id');d=d.set_index('bin_id',drop=False)
 panels={'epithelial':d.epi_marker_supported,'collagen':d.ext_Collagen_stroma_detected>=2}
 for p in ['T/NK','Myeloid','Endothelial','B','Plasma','Mast','Mural']:
  if 'detected_'+p in m:panels[p]=m['detected_'+p].reindex(d.index).fillna(0)>=2
 for radius in [50,200]:
  a=d[d.epi_anchor];out=a[['bin_id']].copy()
  for comp,g in a.groupby('component'):
   pool=d[(d.component==comp)&d.analysis_qc];xy=g[['x_um','y_um']];n=cKDTree(pool[['x_um','y_um']]).query_ball_point(xy,radius+1e-7,return_length=True,workers=4)-1
   for p,flag in panels.items():
    z=pool[flag.reindex(pool.index)];count=cKDTree(z[['x_um','y_um']]).query_ball_point(xy,radius+1e-7,return_length=True,workers=4) if len(z) else np.zeros(len(g));out.loc[g.index,'neighbor_marker_'+p]=np.divide(count-flag.loc[g.index],n,out=np.full(len(g),np.nan),where=n>0)
  out.reset_index(drop=True).to_parquet(o/f'neighbor_marker_r{radius}.parquet',index=False)
 print('MARKER_SUPPORT',s,flush=True)
