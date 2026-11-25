from base import *
from PIL import Image,ImageDraw
import gzip
from scipy import ndimage
from scipy.spatial import cKDTree
from numba import njit
@njit
def clear_paths(i,j,a,b,mask,comp_a,comp_b):
    ok=np.ones(len(i),np.bool_)
    for k in range(len(i)):
        u,v=i[k],j[k]
        if comp_a[u]!=comp_b[v]:
            ok[k]=False; continue
        dx,dy=b[v]-a[u];n=max(1,int(np.ceil(2*max(abs(dx),abs(dy)))))
        for step in range(n+1):
            p=int(np.floor(a[u,0]+dx*step/n+.5)); q=int(np.floor(a[u,1]+dy*step/n+.5))
            if p<0 or q<0 or p>=mask.shape[0] or q>=mask.shape[1] or not mask[p,q]:
                ok[k]=False;break
    return ok

def mask_for(s):
    if s in HD:
        b=pd.read_parquet(A/'hd_epi_niche_dual_methods_v1'/s/'nodes.parquet',columns=['bin_id','array_row','array_col','x_um','y_um','component','in_tissue_original'])
        b=b[b.in_tissue_original==1];rc=b[['array_row','array_col']].to_numpy(int);mask=np.zeros(rc.max(0)+2,bool);mask[tuple(rc.T)]=True
        xy=b[['x_um','y_um']].to_numpy();aff=np.linalg.lstsq(np.c_[xy,np.ones(len(xy))],rc,rcond=None)[0];err=np.max(np.abs(np.c_[xy,np.ones(len(xy))]@aff-rc));assert err<.1,err
        return mask,aff,b
    d=O/'xenium'/s;b=pd.read_parquet(d/'cells.parquet');p=next((O/'xenium/geo').glob('*_'+s+'_cell_boundaries*'))
    with gzip.open(p,'rb') as f: poly=pd.read_parquet(f)
    xy=poly[['vertex_x','vertex_y']].to_numpy();origin=xy.min(0)-40;rc=(xy-origin)/8
    # Using axis 0 for x and axis 1 for y in both PIL and the raster.
    dims=np.ceil(rc.max(0)).astype(int)+6; im=Image.new('1',(int(dims[1]),int(dims[0])))
    draw=ImageDraw.Draw(im)
    for idx in poly.groupby('cell_id',sort=False).indices.values():
        if len(idx)>=3: draw.polygon(list(map(tuple,rc[idx,::-1].round().astype(int))),fill=1)
    mask=ndimage.binary_closing(np.asarray(im),structure=ndimage.iterate_structure(ndimage.generate_binary_structure(2,2),2))
    aff=np.array([[.125,0],[0,.125],-origin/8])
    return mask,aff,b

def pixel(xy,aff):return np.ascontiguousarray(np.c_[xy,np.ones(len(xy))]@aff)
def edges(xa,xb,pa,pb,mask,ca,cb,r,self_graph=False):
    tree=cKDTree(xb);blocks=[]
    for start in range(0,len(xa),256):
        near=tree.query_ball_point(xa[start:start+256],r)
        lens=np.fromiter(map(len,near),int,count=len(near));i=np.repeat(np.arange(start,start+len(near)),lens)
        j=np.concatenate(near).astype(np.int32) if np.sum(lens)>0 else np.empty(0,np.int32)
        if self_graph:
            k=i<j;i=i[k];j=j[k]
        k=clear_paths(i,j,pa,pb,mask,ca,cb);blocks.append((i[k].astype('i4'),j[k]))
    i=np.concatenate([v[0] for v in blocks]);j=np.concatenate([v[1] for v in blocks])
    a=sparse.csr_matrix((np.ones(len(i),np.uint8),(i,j)),shape=(len(xa),len(xb)))
    return a.maximum(a.T).tocsr() if self_graph else a

def build(s):
    mask,aff,b=mask_for(s); is_hd=s in HD; d=O/('hd' if is_hd else 'xenium')/s
    z=pd.read_csv(d/'units.tsv.gz',sep='\t').set_index('region_id') if is_hd else b
    xy=z[['x_um','y_um'] if is_hd else ['x_centroid','y_centroid']].to_numpy()
    rc=pixel(xy,aff); cc=ndimage.label(mask,structure=np.ones((3,3)))[0];ri=np.rint(rc).astype(int)
    comp=z.component.to_numpy('i4') if is_hd else cc[tuple(ri.T)]
    noded=pd.DataFrame({'id':z.index,'x_um':xy[:,0],'y_um':xy[:,1],'component':comp})
    for r in ([50] if is_hd else [50,200]):
        g=edges(xy,xy,rc,rc,mask,comp,comp,r,True);sparse.save_npz(d/f'graph_{r}.npz',g)
        noded[f'degree_{r}']=np.diff(g.indptr)
        if not is_hd:
            ix=np.arange(-int(r/8),int(r/8)+1);disk=(ix[:,None]**2+ix[None,:]**2)*64<=r*r
            from scipy.signal import fftconvolve
            cov=fftconvolve(mask.astype(float),disk.astype(float),mode='same')/disk.sum()
            noded[f'coverage_{r}']=np.clip(cov[tuple(ri.T)],0,1)
        print(s,r,g.nnz,flush=True)
    noded.to_parquet(d/'graph_nodes.parquet',index=False)
    np.savez_compressed(d/'window.npz',mask=mask,affine=aff)
    if not is_hd:return
    e=pd.read_parquet(F/s/'Epi64_regions.parquet',columns=['region_id','x_um','y_um','component','measurement_eligible','H05','MP10','n_bins','nFeature'])
    e.to_parquet(d/'epi64_frozen.parquet',index=False);e=e[e.measurement_eligible].copy()
    high=e.MP10>e.MP10.quantile(.95);assert np.array_equal(high,e.H05)
    em=pd.read_parquet(F/s/'Epi64_members.parquet')
    tm=pd.read_parquet(d/'members.parquet').rename(columns={'region_id':'neighbor_id'})
    sh=em.merge(tm,on='bin_id',validate='one_to_one');sh=sh[sh.region_id.isin(e.region_id)].drop_duplicates(['region_id','neighbor_id'])
    row=pd.Index(e.region_id).get_indexer(sh.region_id);col=z.index.get_indexer(sh.neighbor_id)
    shared=sparse.csr_matrix((np.ones(len(sh),np.uint8),(row,col)),shape=(len(e),len(z)))
    exy=e[['x_um','y_um']].to_numpy();er=pixel(exy,aff);ec=e.component.to_numpy('i4')
    for r in [50,200]:
        g=edges(exy,xy,er,rc,mask,ec,comp,r);over=g.multiply(shared).tocsr();g=(g-over).tocsr();g.eliminate_zeros()
        sparse.save_npz(d/f'epi64_independent_{r}.npz',g)
        cov=np.zeros(len(e))
        for c in np.unique(ec):
            jj=np.flatnonzero(ec==c); points=b.loc[b.component==c,['x_um','y_um']].to_numpy()
            cov[jj]=np.minimum(1,cKDTree(points).query_ball_point(exy[jj],r,return_length=True)*64/(np.pi*r*r))
        e[f'coverage_{r}']=cov;e[f'shared_removed_{r}']=np.diff(over.indptr);e[f'independent_{r}']=np.diff(g.indptr)
    e.to_parquet(d/'epi64_centers.parquet',index=False)
if __name__=='__main__':
    for sample in sys.argv[1:] or HD+XE:build(sample)

