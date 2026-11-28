from data import *
from scipy.spatial import cKDTree
import igraph as ig,leidenalg as la
sys.modules["tensorflow"]=None  # classical UMAP does not use the optional parametric backend
from umap.umap_ import fuzzy_simplicial_set,find_ab_params,simplicial_set_embedding
from sklearn.utils import check_random_state
import importlib.metadata as md
p=sys.argv[1];d=O/p;lat=pd.read_parquet(d/'latent.parquet');z=lat.to_numpy(dtype='f4',copy=True);keys=lat.index
q=pd.read_parquet(d/'qc_eligible.parquet').set_index('cell_key').loc[keys];assert q.index.equals(keys)
n=len(q);k=30;Gpath=d/'expression_graph.npz'
if not Gpath.exists():
 tree=cKDTree(z);neighbors=np.empty((n,k),'i4');dist=np.empty((n,k),'f4')
 for start in range(0,n,4096):
  end=min(n,start+4096);vv,ii=tree.query(z[start:end],k=k,workers=1)
  # UMAP counts the observation itself in its 30-entry list. Make identity explicit in case of exact ties.
  for j in np.flatnonzero(ii[:,0]!=np.arange(start,end)):
   me=start+j;other=ii[j]!=me
   ii[j]=np.r_[me,ii[j,other][:k-1]];vv[j]=np.r_[0.,vv[j,other][:k-1]]
  neighbors[start:end]=ii;dist[start:end]=vv
 G,_,_=fuzzy_simplicial_set(z,k,check_random_state(1729),'euclidean',knn_indices=neighbors,knn_dists=dist,set_op_mix_ratio=1.,local_connectivity=1.)
 G=G.tocsr();G.setdiag(0);G.eliminate_zeros();assert (G!=G.T).nnz==0
 sparse.save_npz(Gpath,G);write(pd.DataFrame({'cell_key':keys}),d/'graph_rows.tsv.gz')
else:
 assert pd.read_csv(d/'graph_rows.tsv.gz',sep='\t').cell_key.tolist()==list(keys);G=sparse.load_npz(Gpath)
tri=sparse.triu(G,k=1).tocoo();graph=ig.Graph(n=n,edges=np.column_stack((tri.row,tri.col)),directed=False);graph.es['weight']=tri.data
for r in R:
 dr=d/('r'+r);dr.mkdir(exist_ok=True);path=dr/'labels.parquet'
 if path.exists():
  labels=pd.read_parquet(path);assert labels.index.equals(keys);q[cname(r)]=labels[cname(r)];continue
 if p=='hd':
  part=la.find_partition(graph,la.RBConfigurationVertexPartition,weights='weight',resolution_parameter=float(r),n_iterations=-1,seed=1729)
 else:
  # Same repeated C++ optimiser and stop test as n_iterations=-1, with an auditable timing/gain trace.
  import time
  part=la.RBConfigurationVertexPartition(graph,weights='weight',resolution_parameter=float(r))
  optimiser=la.Optimiser();optimiser.set_rng_seed(1729)
  with (dr/'optimisation_trace.tsv').open('w') as trace:
   print('iteration','gain','quality','elapsed_s',sep='\t',file=trace)
   start=time.monotonic();it=0
   while True:
    gain=optimiser.optimise_partition(part,n_iterations=1);it+=1
    print(it,gain,part.quality(),round(time.monotonic()-start,3),sep='\t',file=trace,flush=True)
    if gain<=0:break
 m=np.asarray(part.membership);original,counts=np.unique(m,return_counts=True);order=np.lexsort((original,-counts));mapping={int(original[o]):j for j,o in enumerate(order)}
 labels=pd.DataFrame({'original_leiden_id':m,cname(r):np.array([mapping[v] for v in m],dtype='i4')},index=keys);labels.to_parquet(path);q[cname(r)]=labels[cname(r)]
 write(pd.DataFrame({'original_leiden_id':original,'cluster':[mapping[i] for i in original],'n_objects':counts}),dr/'cluster_key.tsv')
 print('saved Leiden',p,r,len(counts),flush=True)
a,b=find_ab_params(1.,.3)
um,_=simplicial_set_embedding(z,G.copy(),2,1.,a,b,1.,5,200,'pca',check_random_state(1729),'euclidean',{},False,{},False,parallel=False)
assert um.shape==(n,2) and np.isfinite(um).all();q['UMAP1']=um[:,0];q['UMAP2']=um[:,1];q.reset_index().to_parquet(d/'objects.parquet',index=False)
(d/'graph_config.json').write_text(json.dumps(dict(neighbors=30,knn='scipy cKDTree exact, including self at input',metric='euclidean',embedding='umap-learn simplicial_set_embedding',umap_version=md.version('umap-learn'),initialization='pca then 200 UMAP SGD epochs',min_dist=.3,spread=1.,a=a,b=b,negative_sample_rate=5,seed=1729,n_edges=G.nnz//2,partition='RBConfigurationVertexPartition',weighted=True,directed=False,iterations=-1),indent=2))

