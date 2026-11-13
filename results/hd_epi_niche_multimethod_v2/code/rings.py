from common import *
from scipy.spatial import cKDTree
if __name__=='__main__':
 gate()
 for s in SECTIONS:
  o=R/'models'/f'{s}_rings.npz'
  if o.exists():continue
  d=pd.read_parquet(V/s/'nodes.parquet');a=d[d.epi_anchor].copy();cols=['supplied','analysis_qc']+[ref+'::'+t for ref in TYPES for t in TYPES[ref]]
  counts_ring=np.zeros((len(a),8,len(cols)),np.int32);a['ordinal']=np.arange(len(a))
  for comp,ac in a.groupby('component'):
   pool=d[d.component==comp];queries=ac[['x_um','y_um']].to_numpy()
   for j,c in enumerate(cols):
    if c=='supplied':z=pool;selfn=np.ones(len(ac),int)
    elif c=='analysis_qc':z=pool[pool.analysis_qc];selfn=np.ones(len(ac),int)
    else:
     ref,t=c.split('::');z=pool[pool.analysis_qc&(pool['type_'+ref]==t)];selfn=(ac['type_'+ref]==t).to_numpy().astype(int)
    tr=cKDTree(z[['x_um','y_um']]);prev=np.zeros(len(ac),int)
    for b,r in enumerate(range(25,201,25)):
     cur=tr.query_ball_point(queries,r+1e-7,return_length=True,workers=4)-selfn if len(z) else np.zeros(len(ac),int)
     counts_ring[ac.ordinal,b,j]=cur-prev;prev=cur
  for ref in TYPES:
   for r,nring in [(50,2),(200,8)]:
    old=pd.read_parquet(V/s/f'r{r}.parquet');ix=[cols.index(ref+'::'+t) for t in TYPES[ref]];assert np.array_equal(counts_ring[:,:nring,ix].sum(1),old[[ref+'_C_'+t for t in TYPES[ref]]].to_numpy())
  known=np.stack([counts_ring[:,:,[cols.index(ref+'::'+t) for t in TYPES[ref]]].sum(2) for ref in TYPES],axis=2);unknown=counts_ring[:,:,1,None]-known
  bounds=np.arange(25,201,25);offsets=np.arange(-25,26)*8;nominal=np.diff([0]+[int(((offsets[:,None]**2+offsets[None,:]**2)<=b*b+1e-7).sum())-1 for b in bounds])
  np.savez_compressed(o,counts=counts_ring,columns=np.array(cols),known=known,unknown=unknown,reference_order=np.array(list(TYPES)),empty_QC=(counts_ring[:,:,1]==0),empty_known=(known==0),nominal_grid_positions=nominal,bin_id=a.bin_id.to_numpy().astype(str));print('RINGS_COMPLETE',s,flush=True)

 from neighbor_evidence import compute
 compute()
