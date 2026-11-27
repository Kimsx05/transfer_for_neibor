# Frozen v1 geometry primitives; source hashes in audit/software_inputs.tsv.
import hashlib
from itertools import product
import numpy as np
import pandas as pd
import shapely
from scipy import sparse
from scipy.sparse.csgraph import connected_components
CATS=[''.join(x) for x in product('EN','IN','SN')]
NAMES=['niche_low_immune','niche_all_epi']
OFFSETS=[(0,2),(1,1),(1,-1),(0,-2),(-1,-1),(-1,1)]
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for part in iter(lambda:f.read(2**20),b''):h.update(part)
    return h.hexdigest()

def tissue_graph(d,valid):
    loc={(r,c,b):i for i,(r,c,b) in enumerate(zip(d.array_row,d.array_col,d.tissue_component_id)) if valid[i]}
    edges=[(i,j) for (r,c,b),i in loc.items() for dr,dc in OFFSETS if (j:=loc.get((r+dr,c+dc,b))) is not None]
    i,j=np.array(edges,dtype=int).reshape(-1,2).T
    return sparse.csr_matrix((np.ones(len(i),np.float32),(i,j)),shape=(len(d),len(d)))

def band(dist,w):
    a=np.asarray(dist,float)
    return np.select([~np.isfinite(a),a< -w,a<=0,a<=w],['NA','B2','B1','G1'],default='G2')

def shared(xy,step):
    mid=xy+step/2
    off=step[:,[1,0]]*np.array([-1,1])/(2*np.sqrt(3))
    return np.stack([mid-off,mid+off],axis=1)

def nearest(xy,lines):
    tree=shapely.STRtree(shapely.linestrings(lines))
    q=shapely.points(xy)
    idx=tree.nearest(q)
    return shapely.distance(q,tree.geometries.take(idx)),idx

def geometry(d,labels,name):
    xy=d[['x_um','y_um']].to_numpy()
    ok=labels!='Unknown'
    inside=np.isin(labels,['ENN','ENS'] if name==NAMES[0] else ['ENN','ENS','EIN','EIS'])
    adj=tissue_graph(d,ok)
    coo=sparse.triu(adj,k=1).tocoo()
    i,j=coo.row,coo.col
    cross=inside[i]!=inside[j]
    same=~cross
    graph=sparse.coo_matrix((np.ones(sum(same)),(i[same],j[same])),shape=adj.shape)
    _,comp=connected_components(graph,directed=False)
    i,j=i[cross],j[cross]
    i,j=np.where(inside[i],i,j),np.where(inside[i],j,i)
    lines=shared(xy[i],xy[j]-xy[i])
    seg=pd.DataFrame(lines.reshape(-1,4),columns=['x1','y1','x2','y2'])
    seg['inside_barcode']=d.barcode.to_numpy()[i];seg['outside_barcode']=d.barcode.to_numpy()[j]
    seg['block']=d.tissue_component_id.to_numpy()[i]
    seg['inside_component']=comp[i];seg['outside_component']=comp[j]
    dist=np.full(len(d),np.nan);unknown=dist.copy()
    nearest_seg=np.full(len(d),-1);nearest_comp=np.full(len(d),-1)
    for block,ids in d[ok].groupby('tissue_component_id').groups.items():
        ids=np.asarray(ids)
        loc=set(zip(d.array_row[ids],d.array_col[ids]))
        pieces=[]
        for dr,dc in OFFSETS:
            miss=[idx for idx,r,c in zip(ids,d.array_row[ids],d.array_col[ids]) if (r+dr,c+dc) not in loc]
            if miss:
                step=np.tile([50*dc,50*np.sqrt(3)*dr],(len(miss),1))
                pieces.append(shared(xy[miss],step))
        if pieces:unknown[ids]=nearest(xy[ids],np.concatenate(pieces))[0]
        for side,key in [(True,'inside_component'),(False,'outside_component')]:
            for group,ss in seg[seg.block==block].groupby(key):
                use=ids[(comp[ids]==group)&(inside[ids]==side)]
                if not len(use):continue
                g,k=nearest(xy[use],lines[ss.index])
                segids=ss.index.to_numpy()[k]
                dist[use]=g*(-1 if side else 1)
                nearest_seg[use]=segids
                nearest_comp[use]=comp[i[segids]]
    strict=np.where(np.abs(dist)+1e-6<unknown,dist,np.nan)
    return dict(d=dist,strict=strict,unknown_edge=unknown,nearest_segment=nearest_seg,nearest_region=nearest_comp,inside=pd.array(np.where(ok,inside.astype(float),np.nan),dtype='Int8')),seg
