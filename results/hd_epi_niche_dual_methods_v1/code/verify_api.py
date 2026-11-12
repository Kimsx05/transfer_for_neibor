import os
os.environ['NUMBA_NUM_THREADS']='4';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,squidpy as sq,json,importlib.metadata as md
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.metrics import adjusted_rand_score
R=Path(__file__).resolve().parents[1]
xy=np.array([[0,0],[8,0],[16,0],[0,100],[8,100],[16,100]],float);labels=pd.Categorical(['Epi','Epi',None,'Other',None,None],categories=['Epi','Other'])
g=cKDTree(xy).sparse_distance_matrix(cKDTree(xy),8.001,output_type='coo_matrix').tocsr();g.data[:]=1;g.setdiag(0);g.eliminate_zeros();a=ad.AnnData(sparse.csr_matrix((6,1)),obs=pd.DataFrame({'type':labels}));a.obsp['spatial_connectivities']=g
from squidpy.gr._niche import _nhood_profile_embedding
p=_nhood_profile_embedding(a,groups='type',spatial_connectivities_key='spatial_connectivities',scale=False,distance=1,abs_nhood=False,n_hop_weights=None)
assert np.allclose(p,[[1,0],[1,0],[1,0],[0,0],[0,1],[0,0]])
# Real ROI: full-node graph vs query-row graph, same selected centers and same random seed.
d=pd.read_parquet(R/'Cervilla/nodes.parquet');c=d[d.epi_anchor].iloc[len(d[d.epi_anchor])//2];z=d[(abs(d.x_um-c.x_um)<256)&(abs(d.y_um-c.y_um)<256)].reset_index(drop=True)
xy=z[['x_um','y_um']].to_numpy();g=cKDTree(xy).sparse_distance_matrix(cKDTree(xy),50.000001,output_type='coo_matrix').tocsr();g.data[:]=1;g.setdiag(0);g.eliminate_zeros();mask=z.epi_anchor.to_numpy();a=ad.AnnData(sparse.csr_matrix((len(z),1)),obs=pd.DataFrame({'type':pd.Categorical(z.type_blca,categories=['Epithelial','T/NK','Myeloid','Fibroblast','Endothelial','Mast','B','Plasma','Mural'])}));a.obsp['spatial_connectivities']=g;b=a.copy();b.obsp['spatial_connectivities']=sparse.diags(mask.astype(float))@g
for obj in [a,b]:sq.gr.calculate_niche_neighborhood(obj,groups='type',resolutions=[.2],n_neighbors=15,scale=False,distance=1,abs_nhood=False,cluster_mask=pd.Series(mask,index=obj.obs_names),rng=123,min_niche_size=None)
assert np.allclose(a.obsm['niche_embedding'][mask],b.obsm['niche_embedding'][mask]);assert adjusted_rand_score(a.obs.loc[mask,'nhood_niche_res=0.2'],b.obs.loc[mask,'nhood_niche_res=0.2'])==1
