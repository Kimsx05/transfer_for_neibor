from common import *
import inspect,importlib.metadata as md,subprocess
import harmonypy,cellcharter as cc,scvi,MENDER,scanpy as sc,squidpy as sq,anndata as ad
from scipy import sparse
from scipy.spatial import distance_matrix
from graph_core import *
versions={p:md.version(p) for p in ['numpy','scipy','pandas','anndata','scanpy','squidpy','scikit-learn','igraph','leidenalg','harmonypy','MENDER','cellcharter','scvi-tools','torch','torchgmm']};versions['python']=__import__('sys').version;versions['MENDER_commit']=subprocess.check_output(['git','-C',str(R/'environment/MENDER'),'rev-parse','HEAD'],text=True).strip();(R/'versions.json').write_text(json.dumps(versions,indent=2))
# Small full geometry: exact-ring implementation matches pinned official MENDER.
xy=np.array([(i*8,j*8) for i in range(12) for j in range(12)],float);ct=np.array(['Epithelial','Fibroblast','Unknown']*48);a=ad.AnnData(sparse.csr_matrix((len(xy),1)),obs=pd.DataFrame({'ct':pd.Categorical(ct),'section':pd.Categorical(['toy']*len(xy))}));a.obsm['spatial']=xy
m=MENDER.MENDER_single(a,ct_obs='ct',random_seed=123);m.set_MENDER_para(n_scales=2,nn_mode='radius',nn_para=25,count_rep='s',include_self=False);m.run_representation()
d=distance_matrix(xy,xy);prev=np.zeros((len(xy),3));types=list(m.ct_unique)
for b,r in enumerate([25,50]):
 c=np.column_stack([((d<=r+1e-7)&(d>0)&(ct[None,:]==t)).sum(1) for t in types]);assert np.array_equal(c-prev,a.obsm[f'scale{b}']);prev=c
# Unknown column removed before normalize_total, not a feature; geometric nodes retained.
# Official CellCharter n_layers=[1] excludes layer zero.
g=sparse.csr_matrix(((d<=50)&(d>0)).astype(np.float32));a.obsp['spatial_connectivities']=g;a.obsm['latent']=np.arange(len(a)*3).reshape(len(a),3).astype(float)
cc.gr.aggregate_neighbors(a,n_layers=[1],aggregations='mean',connectivity_key='spatial_connectivities',use_rep='latent',sample_key='section');expected=g@a.obsm['latent']/np.asarray(g.sum(1));assert np.allclose(a.obsm['X_cellcharter'],expected)
x=np.repeat(np.array([[1.,0],[0,1],[.5,.5]]),[5,3,2],axis=0);gg,u,inv,n=graph(x,k=2);y=partition(gg,inv,.1,123);assert diagnostics(gg,u,inv,n,y)['identical_profiles_split']==0
# Exact quotient preserves the original graph modularity under profile constraints.
full=source_graph(x,k=2);co=sparse.triu(full,k=1).tocoo();ge=ig.Graph(n=len(x),edges=list(zip(co.row.tolist(),co.col.tolist())));yy=np.array([0,1,1]);assert np.isclose(gg.modularity(yy,weights='weight'),ge.modularity(yy[inv],weights=co.data.tolist()))
(R/'tables/implementation_checks.json').write_text(json.dumps(dict(MENDER_exact_ring_counts=True,CellCharter_layer1_mean=True,weighted_contraction_modularity=True,identical_input_split=False,Unknown='counted toy column removed before any preprocessing',source='pinned installed public API and source'),indent=2));print('INTERFACE_CHECKS_PASSED',flush=True)
# Normalization/PCA equivalence to the installed MENDER preprocessing on toy training.
from sklearn.decomposition import PCA
count=np.asarray(m.adata_MENDER.X);names=np.asarray(m.adata_MENDER.var_names).astype(str)
print('MENDER_FEATURE_NAMES',names[:6],flush=True)
known=np.repeat(np.asarray(m.ct_unique)!='Unknown',m.n_scales);raw=count[:,known].copy();v=ad.AnnData(raw.copy());sc.pp.normalize_total(v);sc.pp.log1p(v);target=np.median(raw.sum(1));direct=np.log1p(raw*target/raw.sum(1,keepdims=True));assert np.allclose(v.X,direct,atol=1e-5)
from fit_composition import transfer
xt=np.array([[1.,0.],[1.,0.],[1.,0.],[0.,1.],[0.,1.]]);yy=np.array([0,0,1,2,2]);p,vote,dist=transfer(xt,yy,xt[:1]);assert p[0]==0 and np.isclose(vote[0],2/3) and dist[0]==0
print('PREPROCESS_AND_EXACT_TRANSFER_PASSED',flush=True)
