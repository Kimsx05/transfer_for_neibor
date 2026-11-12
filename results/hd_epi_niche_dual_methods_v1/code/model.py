from pathlib import Path
import os
os.environ['TF_NUM_INTRAOP_THREADS']='4';os.environ['TF_NUM_INTEROP_THREADS']='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]))
for k in ['OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='4'
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd,json,joblib,gc
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score as ari,normalized_mutual_info_score as nmi,silhouette_score
from sklearn.neighbors import NearestNeighbors
import anndata as ad,squidpy as sq
R=Path(__file__).resolve().parents[1];sections=['Cervilla','D1_Johnson','UC_SCC','MPBC1','MPBC2'];types=json.loads((R/'type_mapping.json').read_text());seeds=[123,456,789]
def sample(d,seed=123):
 rng=np.random.default_rng(seed);out=[]
 # round-robin through physical 400um blocks, no label-based sampling
 for s in sections[:4]:
  z=d[(d.section==s)&d.evaluable];blocks=(z.x_um//400).astype(str)+'_'+(z.y_um//400).astype(str)
  groups=[rng.permutation(g.index) for _,g in z.groupby(blocks)];rng.shuffle(groups)
  chosen=[]
  for level in range(max(map(len,groups))):
   chosen.extend(g[level] for g in groups if len(g)>level)
   if len(chosen)>=20000:break
  out.extend(chosen[:20000])
 return np.sort(out)
def km(x,k,seed):return MiniBatchKMeans(k,random_state=seed,n_init=10,batch_size=4096,max_iter=200,max_no_improvement=20,tol=1e-4).fit(x)
def native(P,res,seed):
 # Version-pinned Squidpy's actual native Leiden stage; same representation as public API.
 from squidpy.gr._niche import _leiden_clusterers, calculate_niche_custom
 a=ad.AnnData(P)
 # Public custom pipeline uses identical frozen P rather than recomputing query graph per seed.
 return calculate_niche_custom(a,embedder=lambda _:P,clusterers=_leiden_clusterers(base_colname='nhood_niche',resolutions=res,n_neighbors=15,flavor='igraph',n_iterations=-1),rng=seed,copy=True)
for ref,cats in types.items():
 for radius in [50,200]:
  o=R/f'{ref}_r{radius}';o.mkdir(exist_ok=True)
  if (o/'complete.json').exists():continue
  d=pd.concat([pd.read_parquet(R/s/f'r{radius}.parquet') for s in sections],ignore_index=True);d['evaluable']=d[ref+'_n_known']>=10;P=d[[ref+'_P_'+t for t in cats]].to_numpy();T=sample(d);S=np.flatnonzero(d.evaluable);assert np.isfinite(P[S]).all() and np.allclose(P[S].sum(1),1)
  d.loc[T,['section','bin_id']].to_csv(o/'training_anchor_ids.tsv',sep='\t',index=False)
  # All nodes retained; non-training graph rows empty. Each T row has its complete true radius neighborhood.
  graph=[];obs=[];offset=0;api_indices=[]
  for s in sections:
   nodes=pd.read_parquet(R/s/'nodes.parquet');ns=len(nodes);index=pd.Series(np.arange(ns),index=nodes.bin_id);ts=d.loc[T];ts=ts[ts.section==s];ti=index.loc[ts.bin_id].to_numpy();api_indices.extend(offset+ti);offset+=ns
   rows=[];cols=[]
   for comp,g in ts.groupby('component'):
    pool=nodes[(nodes.component==comp)&nodes.analysis_qc];js=index.loc[pool.bin_id].to_numpy();tree=cKDTree(pool[['x_um','y_um']]);ci=index.loc[g.bin_id].to_numpy()
    for a in range(0,len(g),256):
     lists=tree.query_ball_point(g.iloc[a:a+256][['x_um','y_um']],radius+1e-7,workers=4)
     for ii,ls in zip(ci[a:a+256],lists):
      jj=js[ls];jj=jj[jj!=ii];rows.append(np.full(len(jj),ii,np.int32));cols.append(jj.astype(np.int32))
   rr=np.concatenate(rows) if rows else np.array([],np.int32);cc=np.concatenate(cols) if cols else np.array([],np.int32)
   graph.append(sparse.csr_matrix((np.ones(len(rr),np.float32),(rr,cc)),shape=(ns,ns)));obs.append(nodes[['bin_id','type_'+ref]].assign(section=s));del nodes,rows,cols;gc.collect()
  ob=pd.concat(obs,ignore_index=True);ob.index=ob.section+'::'+ob.bin_id;ob['type']=pd.Categorical(ob['type_'+ref],categories=cats)
  a=ad.AnnData(sparse.csr_matrix((len(ob),1)),obs=ob);a.obsp['spatial_connectivities']=sparse.block_diag(graph,format='csr');del graph,obs;gc.collect()
  mask=pd.Series(False,index=a.obs_names);mask.iloc[api_indices]=True
  # API fits masked rows in observation order, recover by unique identity.
  resolutions=[.1,.2,.4,.6,.8,1.0]
  sq.gr.calculate_niche_neighborhood(a,groups='type',resolutions=resolutions,n_neighbors=15,scale=False,distance=1,abs_nhood=False,cluster_mask=mask,min_niche_size=None,flavor='igraph',n_iterations=-1,rng=123)
  train_names=(d.loc[T,'section']+'::'+d.loc[T,'bin_id']).tolist();ix=a.obs_names.get_indexer(train_names);assert (ix>=0).all();assert np.allclose(a.obsm['niche_embedding'][ix],P[T],atol=1e-6)
  B={res:a.obs.loc[train_names,f'nhood_niche_res={res}'].astype(str).to_numpy() for res in resolutions};del a;gc.collect()
  rng=np.random.default_rng(123);evalix=rng.choice(len(T),min(2500,len(T)),replace=False);pt=P[T];zt=np.sqrt(pt);metrics=[];models={};labels={}
  for k in range(2,9):
   fit=km(zt,k,123);y=fit.predict(zt);stab=[ari(y,km(zt,k,seed).predict(zt)) for seed in seeds[1:]]
   # Fixed alternate physical-block holdout diagnostic; no new training weights.
   block=((d.loc[T,'x_um']//400+d.loc[T,'y_um']//400).astype(int)%5!=0).to_numpy();yst=km(zt[block],k,456).predict(zt)
   sil=silhouette_score(zt[evalix],y[evalix]);metrics.append(dict(method='A',parameter=k,n_clusters=len(np.unique(y)),seed_ARI=np.mean(stab),block_ARI=ari(y,yst),silhouette=sil));models[k]=fit;labels[k]=y
  # Native stage repeated on exactly P_T (public API above is the baseline).
  bs={seed:native(pt,resolutions,seed) for seed in seeds[1:]};block=((d.loc[T,'x_um']//400+d.loc[T,'y_um']//400).astype(int)%5!=0).to_numpy();bb=native(pt[block],resolutions,456)
  for res,y in B.items():
   col=f'nhood_niche_res={res}';stab=[ari(y,bs[seed].obs[col].astype(str)) for seed in seeds[1:]];sil=silhouette_score(pt[evalix],y[evalix]) if len(np.unique(y[evalix]))>1 else -1
   metrics.append(dict(method='B',parameter=res,n_clusters=len(np.unique(y)),seed_ARI=np.mean(stab),block_ARI=ari(y[block],bb.obs[col].astype(str)),silhouette=sil))
  met=pd.DataFrame(metrics);met['selection_score']=.5*met.seed_ARI+.3*met.block_ARI+.2*met.silhouette;met.to_csv(o/'candidate_metrics.tsv',sep='\t',index=False)
  ak=int(met[met.method=='A'].sort_values(['selection_score','parameter'],ascending=[False,True]).iloc[0].parameter);br=float(met[met.method=='B'].sort_values(['selection_score','parameter'],ascending=[False,True]).iloc[0].parameter)
  joblib.dump(models[ak],o/'A_model.joblib');by=B[br];np.savez_compressed(o/'B_training.npz',P=pt,labels=by.astype(str),T=T)
  result=d[['section','bin_id','evaluable']].copy();coverage=[]
  for method in ['A','B']:
   candidate=np.full(len(d),'',object);final=np.full(len(d),'unresolved',object);source=np.full(len(d),'unresolved',object);reason=np.full(len(d),'n_known_below_10',object);support=np.full(len(d),np.nan);distance=np.full(len(d),np.nan)
   if method=='A':
    fit=models[ak];ds=fit.transform(np.sqrt(P[S]));pred=fit.predict(np.sqrt(P[S])).astype(str);near=np.sort(ds,axis=1)[:,:2];vote=(near[:,1]-near[:,0])/np.maximum(near[:,1],1e-12);dist=near[:,0];td=fit.transform(zt).min(1);threshold=np.quantile(td,.995)*1.5+1e-8;ok=dist<=threshold;ty=labels[ak].astype(str)
   else:
    nn=NearestNeighbors(n_neighbors=15,n_jobs=4).fit(pt);pred=[];votes=[];dist=[]
    for st in range(0,len(S),10000):
     dd,jj=nn.kneighbors(P[S[st:st+10000]])
     for ds,js in zip(dd,jj):
      sel=js[ds<1e-10] if (ds<1e-10).any() else js;u,c=np.unique(by[sel],return_counts=True);chosen=sorted(zip(-c,u))[0][1];pred.append(chosen);votes.append(c.max()/len(sel));dist.append(ds[0])
    pred=np.asarray(pred);vote=np.asarray(votes);dist=np.asarray(dist);held=~block;hdd,_=NearestNeighbors(n_neighbors=1).fit(pt[block]).kneighbors(pt[held]);threshold=np.quantile(hdd,.995)*1.5+1e-8;ok=(dist<=threshold)&(vote>.5);ty=by
   candidate[S]=pred;final[S[ok]]=pred[ok];source[S[ok]]='model_predict' if method=='A' else 'external_projection';reason[S]='outside_training_support_or_vote_tie';reason[S[ok]]='assigned';support[S]=vote;distance[S]=dist
   final[T]=ty;candidate[T]=ty;source[T]='native_train';reason[T]='native_fit';result[method+'_candidate']=candidate;result[method+'_niche']=final;result[method+'_assignment_source']=source;result[method+'_reason']=reason;result[method+'_support']=support;result[method+'_distance']=distance
   for sec,g in result.groupby('section'):
    for state,n in g[method+'_assignment_source'].value_counts().items():coverage.append(dict(section=sec,method=method,state=state,n=n,anchors=len(g),evaluable=int(g.evaluable.sum())))
  result.to_parquet(o/'labels.parquet',index=False);pd.DataFrame(coverage).to_csv(o/'coverage.tsv',sep='\t',index=False)
  (o/'complete.json').write_text(json.dumps(dict(ref=ref,radius=radius,K=ak,resolution=br,n_train=len(T),n_evaluable=len(S),native_ARI=ari(labels[ak],by),native_NMI=nmi(labels[ak],by),API_composition_verified=True),indent=2));print('COMPLETE',ref,radius,flush=True)
