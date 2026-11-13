from common import *
from scipy.spatial import cKDTree

def compute():
 for s in SECTIONS:
  n=pd.read_parquet(V/s/'nodes.parquet');a=n[n.epi_anchor].copy()
  for r in [50,200]:
   dest=R/'tables'/f'{s}_raw_panel_neighbors_r{r}.parquet'
   if dest.exists():continue
   out=a[['bin_id','section']].copy().set_index('bin_id');den=pd.read_parquet(V/s/f'r{r}.parquet').set_index('bin_id').loc[out.index,'n_analysis_neighbors']
   for name,col in [('basal_squamous','ext_Keratin_extended_detected'),('urothelial','ext_Urothelial_controls_detected')]:
    for comp,q in a.groupby('component'):
     pool=n[n.analysis_qc&(n.component==comp)&(n[col]>=2)];c=cKDTree(pool[['x_um','y_um']]).query_ball_point(q[['x_um','y_um']],r+1e-7,return_length=True,workers=4) if len(pool) else np.zeros(len(q),int);c=c-(q[col]>=2).to_numpy();out.loc[q.bin_id,'neighbor_marker_'+name]=c/den.loc[q.bin_id].to_numpy()
   out.reset_index().to_parquet(dest,index=False)
 print('RAW_PANEL_NEIGHBOR_EVIDENCE_COMPLETE',flush=True)
if __name__=='__main__':gate();compute()
