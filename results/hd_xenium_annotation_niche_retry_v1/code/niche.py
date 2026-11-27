from base_data import *
from scipy.spatial import cKDTree
from scipy.ndimage import label
from scipy.sparse.csgraph import breadth_first_order
from sklearn.metrics import adjusted_rand_score
import itertools,torch,cellcharter as cc
from geometry import clear_paths,pixel
import pickle
torch.set_num_threads(2);torch.set_num_interop_threads(1)
platform=sys.argv[1];samples=HD if platform=='hd' else XE
zs=[];training=[];order=[];rng=np.random.default_rng(1729);rows=[]
for s in samples:
 d=paths(s);a=dataset(s);latent=np.load(d/'latent.npy');idx=np.flatnonzero(np.isfinite(latent[:,0])&eligible(a,s).to_numpy());xy=a.obs[['x_um','y_um']].values
 w=np.load(V/platform/s/'window.npz');mask=w['mask'];rc=pixel(xy,w['affine']);pix=np.floor(rc+.5).astype(int);inside=((pix>=0)&(pix<mask.shape)).all(1);comp=np.zeros(len(xy),'i4');comp[inside]=label(mask,np.ones((3,3)))[0][tuple(pix[inside].T)]
 # Complete observed section before eligibility filtering. Subsetting removes edges without creating new neighbors.
 dist,nb=cKDTree(xy).query(xy,k=7,workers=1);ii=np.repeat(np.arange(len(xy)),7);jj=nb.ravel();ok=(ii!=jj)&(dist.ravel()<=30);ii=ii[ok];jj=jj[ok];keep=clear_paths(ii,jj,rc,rc,mask,comp,comp)
 G=sparse.csr_matrix((np.ones(keep.sum(),'f4'),(ii[keep],jj[keep])),shape=(len(xy),len(xy)));G=G.maximum(G.T);G.setdiag(0);G.eliminate_zeros();G=G[idx][:,idx];deg=np.diff(G.indptr)
 geo=xy[idx];shells=[[] for _ in range(3)]
 for node in rng.choice(len(idx),min(2000,len(idx)),replace=False):
  front={node};seen=front.copy()
  for h in range(3):
   front={j for k in front for j in G.indices[G.indptr[k]:G.indptr[k+1]]}-seen;seen.update(front)
   if front:shells[h].extend(np.linalg.norm(geo[list(front)]-geo[node],axis=1))
 for h,r in enumerate(shells):
  rows.append((s,len(xy),len(idx),G.nnz//2,int((deg==0).sum()),h+1,len(r),*np.quantile(r,[.05,.5,.95])))
 sparse.save_npz(d/'niche_graph.npz',G);pd.Series(a.obs_names[idx]).to_csv(d/'niche_graph_ids.tsv.gz',header=False,index=False)
 obs=pd.DataFrame({'sample':pd.Categorical([s]*len(idx))},index=a.obs_names[idx]);q=ad.AnnData(sparse.csr_matrix((len(idx),1)),obs=obs);q.obsm['latent']=latent[idx];q.obsp['spatial_connectivities']=G
 cc.gr.aggregate_neighbors(q,n_layers=[0,1,2,3],sample_key='sample',connectivity_key='spatial_connectivities',use_rep='latent')
 feat=q.obsm['X_cellcharter'].astype('f4');zs.append(feat);order.append((s,a.obs_names,idx,deg));training.append(rng.choice(len(idx),min(3200,min([len(idx),3200])),replace=False)+sum(map(len,zs[:-1])));del a,q,G
write(pd.DataFrame(rows,columns=['sample','all_objects','eligible','edges','isolates','hop','pairs_in_distance_audit','p05_um','median_um','p95_um']),platform+'_graph.tsv')
X=np.vstack(zs);train=np.concatenate(training);np.save(O/platform/'gmm_training_indices.npy',train)
stability=[];models={}
for k in [4,6,8]:
 pred=[];conv=[]
 for seed in [1729,1733,1741]:
  m=cc.tl.GaussianMixture(n_clusters=k,random_state=seed,batch_size=1024,trainer_params={'accelerator':'cpu','devices':1,'max_epochs':80,'enable_progress_bar':False})
  m.fit(X[train]);pp=np.asarray(m.predict(X[train]));pred.append(pp);conv.append(bool(m.converged_))
  if seed==1729:models[k]=m
 ari=np.median([adjusted_rand_score(a,b) for a,b in itertools.combinations(pred,2)])
 stability.append((k,ari,min(np.bincount(pred[0],minlength=k))/len(train),sum(conv)))
tab=pd.DataFrame(stability,columns=['K','median_seed_ARI','min_balanced_fraction','n_converged']);write(tab,platform+'_niche_stability.tsv');cand=tab[tab.min_balanced_fraction>=.015];cand=tab if cand.empty else cand
k=int(cand.sort_values(['median_seed_ARI','K'],ascending=[False,True]).iloc[0].K);m=models[k]
with open(O/platform/'cellcharter_model.pkl','wb') as f:pickle.dump(m,f)
off=0
for s,allids,idx,deg in order:
 feat=X[off:off+len(idx)];off+=len(idx);p=np.concatenate([np.asarray(m.predict(feat[i:i+2048])) for i in range(0,len(idx),2048)])
 out=pd.DataFrame(index=allids);out['niche_id']='unassigned_QC';out.iloc[idx,0]=[f'{platform}_N{j+1}' for j in p];out['graph_degree']=pd.array([pd.NA]*len(out),dtype='Int32');out.iloc[idx,1]=deg;out.to_parquet(paths(s)/'niche.parquet')

