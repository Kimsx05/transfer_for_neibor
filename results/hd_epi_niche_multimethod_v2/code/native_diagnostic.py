from common import *
import scanpy as sc,anndata as ad
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import NearestNeighbors
if __name__=='__main__':
 gate()
 if (R/'tables/native_graph_diagnostics.tsv').exists():print('NATIVE_DIAGNOSTIC_CACHE',flush=True);raise SystemExit(0)
 rows=[]
 for ref in TYPES:
  for r in [50,200]:
   z=np.load(V/f'{ref}_r{r}/B_training.npz');p=z['P'];y=z['labels'];a=ad.AnnData(p);sc.pp.neighbors(a,n_neighbors=15,use_rep='X',random_state=123);g=a.obsp['connectivities'];u,inv,n=np.unique(p,axis=0,return_inverse=True,return_counts=True);codes=pd.Categorical(y).codes;lo=np.full(len(u),codes.max()+1);hi=np.full(len(u),-1);np.minimum.at(lo,inv,codes);np.maximum.at(hi,inv,codes)
   dd,_=NearestNeighbors(n_neighbors=2,n_jobs=4).fit(u).kneighbors(u);rows.append(dict(reference=ref,radius=r,training=len(p),unique_profiles=len(u),zero_tie_rows=int(n[n>1].sum()),max_exact_multiplicity=int(n.max()),near_unique_profiles_lt_1e_6=int((dd[:,1]<1e-6).sum()),components=int(connected_components(g,directed=False,return_labels=False)),isolates=int((g.getnnz(1)==0).sum()),degree_median=float(np.median(g.getnnz(1))),edge_weight_q01=float(np.quantile(g.data,.01)),edge_weight_median=float(np.median(g.data)),native_clusters=len(set(y)),pure_Epi_training_rows=int((p[:,0]==1).sum()),pure_Epi_native_labels=len(set(y[p[:,0]==1])),identical_profiles_split=int((lo!=hi).sum()),split_rows=int(n[lo!=hi].sum())))
 pd.DataFrame(rows).to_csv(R/'tables/native_graph_diagnostics.tsv',sep='\t',index=False);print('NATIVE_DIAGNOSTIC_COMPLETE',flush=True)
