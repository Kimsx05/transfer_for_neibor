import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
import scanpy as sc,anndata as ad
import igraph as ig,leidenalg as la

def source_graph(x,k=15):
 a=ad.AnnData(np.asarray(x,dtype=np.float64));sc.pp.neighbors(a,n_neighbors=min(k,len(x)-1),use_rep='X',random_state=123);g=a.obsp['connectivities'].tocsr();g.setdiag(0);g.eliminate_zeros();return g

def graph(x,k=15):
 # Exact weighted quotient of the unchanged Scanpy full-training composition graph.
 # Membership constraint: all rows of one complete identical input profile together.
 # Q=G.T A G: off-diagonal weights count each undirected cross-group edge once;
 # Q_ii/2 self-loop is the internal edge mass. No invented clique or n_i*n_j edges.
 x=np.asarray(x,dtype=np.float64);u,inv,n=np.unique(x,axis=0,return_inverse=True,return_counts=True);full=source_graph(x,k);membership=sparse.csr_matrix((np.ones(len(x)),(np.arange(len(x)),inv)),shape=(len(x),len(u)));q=(membership.T@full@membership).tocsr();co=sparse.triu(q,k=1).tocoo();edges=list(zip(co.row.tolist(),co.col.tolist()));weights=co.data.tolist()
 for i,w in enumerate(q.diagonal()):
  if w>0:edges.append((i,i));weights.append(float(w/2))
 g=ig.Graph(n=len(u),edges=edges,directed=False);g.es['weight']=weights;g.vs['frequency']=n.tolist();g['source_components']=int(connected_components(full,directed=False,return_labels=False));g['source_edge_mass']=float(full.sum()/2);g['source_degree_median']=float(np.median(full.getnnz(1)));assert np.isclose(sum(weights),g['source_edge_mass'],rtol=1e-6);return g,u,inv,n

def partition(g,inv,res,seed):
 p=la.find_partition(g,la.RBConfigurationVertexPartition,weights='weight',resolution_parameter=float(res),seed=int(seed),n_iterations=-1);return np.asarray(p.membership)[inv]

def diagnostics(g,u,inv,n,y):
 deg=np.array(g.degree());w=np.array(g.es['weight']);lo=np.full(len(u),int(y.max())+1);hi=np.full(len(u),-1);np.minimum.at(lo,inv,y);np.maximum.at(hi,inv,y);split=int(((lo!=hi)&(n>1)).sum())
 return dict(n=len(inv),unique_profiles=len(u),duplicate_rows=int(n[n>1].sum()),identical_profiles_split=split,components=len(g.connected_components()),isolates=int((deg==0).sum()),degree_median=float(np.median(deg)),degree_max=int(deg.max()),edge_weight_median=float(np.median(w)) if len(w) else 0,edge_weight_max=float(w.max()) if len(w) else 0,source_components=g['source_components'],source_degree_median=g['source_degree_median'],source_edge_mass=g['source_edge_mass'],quotient_edge_mass=float(w.sum()),graph='exact weighted quotient; identical-profile constraint')
