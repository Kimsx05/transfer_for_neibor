from common import *
from geometry import pixel,edges,clear_paths
from numba import njit
from scipy.ndimage import label,distance_transform_edt
from scipy.spatial import cKDTree
import shapely
@njit
def area_obs(rc,comp,r,cc,dist):
 out=np.zeros(len(rc))
 for i,(a,b) in enumerate(rc):
  rr=int(np.rint(a));ss=int(np.rint(b));c=comp[i]
  if c<=0 or not(0<=rr<cc.shape[0] and 0<=ss<cc.shape[1]) or cc[rr,ss]!=c:continue
  if dist[rr,ss]*8>r+12:out[i]=np.pi*r*r;continue
  for x in range(max(0,int(np.ceil(a-r/8))),min(cc.shape[0],int(np.floor(a+r/8))+1)):
   for y in range(max(0,int(np.ceil(b-r/8))),min(cc.shape[1],int(np.floor(b+r/8))+1)):
    if (x-a)**2+(y-b)**2>(r/8)**2 or cc[x,y]!=c:continue
    n=max(1,int(np.ceil(2*max(abs(x-a),abs(y-b)))));ok=True
    for k in range(1,n):
     if cc[int(np.floor(a+.5+(x-a)*k/n)),int(np.floor(b+.5+(y-b)*k/n))]!=c:ok=False;break
    if ok:out[i]+=64
  out[i]=min(np.pi*r*r,out[i])
 return out
def polygons(s,ids):
 with gzip.open(next((V/'xenium/geo').glob('*_'+s+'_cell_boundaries*')),'rb') as f:b=pd.read_parquet(f)
 b=b.sort_values('cell_id',kind='stable');ix,order=pd.factorize(b.cell_id,sort=False)
 arr=shapely.polygons(shapely.linearrings(b[['vertex_x','vertex_y']].to_numpy(),indices=ix))
 ix=pd.Index(order).get_indexer(ids);assert min(ix)>=0;return arr[ix]
def graph(s):
 hd=s in HD;v,d=folders(s);z=labels(s);gn=pd.read_parquet(v/'graph_nodes.parquet').set_index('id').loc[z.index]
 xy=z[['x_um','y_um']].values;w=np.load(v/'window.npz');mask=w['mask'];aff=w['affine'];rc=pixel(xy,aff);cp=gn.component.to_numpy('i4');cc=label(mask,structure=np.ones((3,3)))[0]
 if hd:
  b=pd.read_parquet(ROOT/'highres_VisiumHD/analysis/hd_epi_niche_dual_methods_v1'/s/'nodes.parquet',columns=['in_tissue_original','array_row','array_col','component']);b=b[b.in_tissue_original==1];cc[:]=0
  # Original component identifiers can start at zero; offset both maps and IDs.
  cc[b.array_row,b.array_col]=b.component+1;cp=cp+1
 dist=distance_transform_edt(mask)
 for res in ([16,32] if hd else ['cell']):
  cs=pd.read_parquet(d/f'states_{res}.parquet');ee=cs[cs.eligible].copy();ex=ee[['x_um','y_um']].values;ec=ee.component.to_numpy('i4')+(1 if hd else 0);er=pixel(ex,aff)
  if res==32:parents=pd.read_parquet(d/'parent32.parquet').loc[z.index].region_id32
  ei=z.index.get_indexer(ee.index)
  for r in [50,200]:
   if res!=32 and (v/f'graph_{r}.npz').exists():
    full=sparse.load_npz(v/f'graph_{r}.npz');g=full[ei];del full
   else:g=edges(ex,xy,er,rc,mask,ec,cp,r)
   i,j=g.nonzero();k=(parents.to_numpy()[j]==ee.index.to_numpy()[i]) if res==32 else (ei[i]==j)
   shared=int(k.sum());g=sparse.csr_matrix((np.ones((~k).sum(),np.int8),(i[~k],j[~k])),shape=g.shape)
   assert (ec[g.nonzero()[0]]==cp[g.indices]).all()
   area=area_obs(er,ec,r,cc,dist)
   if hd:area=np.maximum(0,area-64*ee.n_bins.to_numpy())
   ee[f'area_{r}']=area;ee[f'coverage_{r}']=area/(np.pi*r*r);ee[f'neighbors_{r}']=np.diff(g.indptr)
   sparse.save_npz(d/f'epi_{res}_{r}.npz',g)
   print(s,res,r,len(ee),g.nnz,shared,flush=True)
  ee.to_parquet(d/f'centers_{res}.parquet')
 gr=pd.read_parquet(d/'identity.parquet');rr=gr.grade.isin(['supported','possible']).to_numpy()
 tq=pd.DataFrame(index=z.index[rr]);tq['component']=cp[rr]
 for r in [50,200]:
  ar=area_obs(rc[rr],cp[rr],r,cc,dist)
  if hd:ar=np.maximum(0,ar-z.n_bins.to_numpy()[rr]*64)
  tq[f'coverage_{r}']=ar/(np.pi*r*r)
 tq.to_parquet(d/'target_cover.parquet')
 if not hd:
  po=polygons(s,z.index);valid=~(pd.read_parquet(v/'boundary_geometry.parquet').loc[z.index].review_boundary|z.area_outlier|z.boundary_degenerate).to_numpy()
  bi=shapely.boundary(po);ei=z.index.get_indexer(ee.index);ix,j=shapely.STRtree(bi).query(bi[ei],predicate='dwithin',distance=5)
  ok=valid[ei[ix]]&valid[j]&(ei[ix]!=j)&clear_paths(ei[ix],j,rc,rc,mask,cp,cp)
  ix,j=ix[ok],j[ok];gap=shapely.distance(bi[ei[ix]],bi[j]);pd.DataFrame({'row':ix,'neighbor':j,'gap_um':gap}).to_parquet(d/'boundary_pairs.parquet',index=False)
if __name__=='__main__':
 for s in sys.argv[1:] or HD+XE:graph(s)
