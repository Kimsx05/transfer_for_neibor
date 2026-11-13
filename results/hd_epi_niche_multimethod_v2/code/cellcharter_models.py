from common import *
from scipy.spatial import cKDTree
from scipy import sparse
import anndata as ad,cellcharter as cc,torch,joblib,gc
from sklearn.metrics import adjusted_rand_score as ari,silhouette_score
from sklearn.neighbors import NearestNeighbors
import concurrent.futures
_Executor=concurrent.futures.ThreadPoolExecutor
class LimitedExecutor(_Executor):
 def __init__(self,max_workers=None,*args,**kwargs):super().__init__(min(max_workers or 4,4),*args,**kwargs)
concurrent.futures.ThreadPoolExecutor=LimitedExecutor
from fit_composition import select

def aggregate():
 for s in SECTIONS:
  if all((R/'models'/f'{s}_expression_neighbors_r{r}.npz').exists() for r in [50,200]):continue
  n=pd.read_parquet(V/s/'nodes.parquet');a=n[n.epi_anchor].copy();a['ordinal']=np.arange(len(a));pool=n[n.analysis_qc].copy();pool['ordinal']=np.arange(len(pool));zs=[]
  for rep in ['scVI','scVI_without_Top20']:
   f=np.load(R/'models/expression'/rep/f'{s}_latent.npz');ix=pd.Index(f['bin_id']).get_indexer(pool.bin_id);assert (ix>=0).all();zs.append(f['latent'][ix])
  zs=np.column_stack(zs)
  for r in [50,200]:
   dest=R/'models'/f'{s}_expression_neighbors_r{r}.npz'
   if dest.exists():continue
   X=np.zeros((len(a),20),np.float32);num=np.zeros(len(a),int)
   for comp,aa in a.groupby('component'):
    pp=pool[pool.component==comp];tr=cKDTree(pp[['x_um','y_um']]);pi=pp.ordinal.to_numpy();self_index=pd.Index(pp.bin_id)
    for st in range(0,len(aa),256):
     ac=aa.iloc[st:st+256];ls=tr.query_ball_point(ac[['x_um','y_um']],r+1e-7,workers=4);rr=[];ccs=[]
     for j,(li,own) in enumerate(zip(ls,self_index.get_indexer(ac.bin_id))):
      q=np.asarray(li);q=q[q!=own];rr.extend([j]*len(q));ccs.extend(q.tolist())
     g=sparse.csr_matrix((np.ones(len(rr),np.float32),(rr,ccs)),shape=(len(ac),len(pp)));nc=np.asarray(g.sum(1)).ravel();X[ac.ordinal]=np.divide(g@zs[pi],nc[:,None],out=np.zeros((len(ac),20),np.float32),where=nc[:,None]>0);num[ac.ordinal]=nc
   old=pd.read_parquet(V/s/f'r{r}.parquet');assert np.array_equal(num,old.n_analysis_neighbors);np.savez_compressed(dest,X=X,n=num,bin_id=a.bin_id.to_numpy().astype(str));print('EXPRESSION_NEIGHBORS',s,r,flush=True)

def array(v):return v.detach().cpu().numpy() if hasattr(v,'detach') else np.asarray(v)
def fit(r,rep):
 o=R/'models'/f'CellCharter_both_r{r}_{rep}';o.mkdir(exist_ok=True)
 if (o/'complete.json').exists():return
 if rep!='scVI' and json.loads((R/'models/expression'/rep/'complete.json').read_text()).get('identical_to')=='scVI':
  main=R/'models'/f'CellCharter_both_r{r}_scVI';d=pd.read_parquet(main/'assignments.parquet');d['representation']=rep;d.to_parquet(o/'assignments.parquet',index=False);j=json.loads((main/'complete.json').read_text());j.update(representation=rep,identity_sensitivity=True,reason='Top20 absent from main model features; no duplicate training');(o/'complete.json').write_text(json.dumps(j,indent=2));
  for name in ['model','train_labels.npy','buffer_validation.json','candidate_metrics.tsv']:
   f=o/name
   if not f.exists():f.symlink_to(Path('..')/main.name/name,target_is_directory=name=='model')
  return
 a=anchors();T=train_index(a);x=[];nc=[]
 for s in SECTIONS:
  z=np.load(R/'models'/f'{s}_expression_neighbors_r{r}.npz');x.append(z['X'][:,:10] if rep=='scVI' else z['X'][:,10:]);nc.append(z['n'])
 x=np.concatenate(x);nc=np.concatenate(nc);ok=nc>=10;assert ok[T].all();S=np.flatnonzero(ok);xt=x[T];adata=ad.AnnData(np.zeros((len(T),1)));adata.obsm['X_cellcharter']=xt
 device='gpu' if torch.cuda.is_available() else 'cpu';params=dict(batch_size=8192,trainer_params=dict(accelerator=device,devices=1,max_epochs=100,enable_progress_bar=False),random_state=123);torch.set_num_threads(4)
 if rep=='scVI':
  aut=cc.tl.ClusterAutoK(n_clusters=list(range(2,11)),max_runs=3,convergence_tol=.01,model_params=params);aut.fit(adata)
  (o/'autok_run_info.json').write_text(json.dumps(dict(candidate_K=list(range(2,11)),native_auxiliary_K=[1,11],derived_run_seeds=[int(v.generate_state(1)[0]) for v in np.random.SeedSequence(123).spawn(3)],actual_runs=len(aut.labels[2]),max_runs=3,thread_pool_limit=4),indent=2))
  rows=[];ix=np.random.default_rng(123).choice(len(T),2500,replace=False)
  for j,k in enumerate(aut.n_clusters):
   ys=[array(v) for v in aut.labels[k]];f=aut.best_models[k];y=array(f.predict(xt));rows.append(dict(parameter=k,n_clusters=len(set(y)),seed_ARI=np.mean([ari(ys[0],v) for v in ys[1:]]),silhouette=silhouette_score(xt[ix],y[ix]),AutoK_stability=float(np.mean(aut.stability[j])),actual_runs=len(ys),converged=bool(f.converged_),iterations=int(f.num_iter_)))
  m=pd.DataFrame(rows);best=m.AutoK_stability.max();cand=m[(m.AutoK_stability>=best-.05)&(m.seed_ARI>=.8)];cand=cand if len(cand) else m[m.AutoK_stability>=best-.05];chosen=cand.sort_values('parameter').iloc[0];k=int(chosen.parameter);f=aut.best_models[k];aut.best_models[2].save(o/'K2_model');np.savez_compressed(o/'K2.npz',S=S,labels=array(aut.best_models[2].predict(x[S])));m['selected']=m.parameter==k;m.to_csv(o/'candidate_metrics.tsv',sep='\t',index=False)
 else:
  main=json.loads((R/'models'/f'CellCharter_both_r{r}_scVI/complete.json').read_text());k=int(main['selected_parameter']);f=cc.tl.GaussianMixture(n_clusters=k,**params);f.fit(xt);pd.DataFrame([dict(parameter=k,selected=True,reason='fixed main K, refitted in new expression space',converged=bool(f.converged_))]).to_csv(o/'candidate_metrics.tsv',sep='\t',index=False)
 y=array(f.predict(xt));pr=array(f.predict_proba(x[S]));pred=pr.argmax(1);confidence=pr.max(1)
 keep,test,bd=buffer_split(a,T,r);near=NearestNeighbors(n_neighbors=1,n_jobs=4).fit(xt[keep]);hd=near.kneighbors(xt[test])[0][:,0];threshold=float(np.quantile(hd,.995)*1.5+1e-8);dist=NearestNeighbors(n_neighbors=1,n_jobs=4).fit(xt).kneighbors(x[S])[0][:,0];valid=(dist<=threshold)&(confidence>.5)
 bf=cc.tl.GaussianMixture(n_clusters=k,**{**params,'random_state':456});bf.fit(xt[keep]);buf=dict(retained_train=int(keep.sum()),validation=int(test.sum()),minimum_cross_distance_um=float(bd[keep].min()),required_gt_um=2*r,ARI_to_full_on_buffer_test=ari(y[test],array(bf.predict(xt[test]))),scope='fixed expression/reference-query model conditional validation')
 z=a[['section','bin_id']].copy();z['candidate']=-1;z['label']=-1;z['source']='unresolved';z['reason']='insufficient_QC_expression_neighbors';z['confidence']=np.nan;z['distance']=np.nan;z.loc[S,'candidate']=pred;z.loc[S,'confidence']=confidence;z.loc[S,'distance']=dist;z.loc[S,'reason']='GMM_membership_or_extrapolation';z.loc[S[valid],'label']=pred[valid];z.loc[S[valid],'source']='native_model_predict';z.loc[S[valid],'reason']='assigned';z.loc[T,'label']=y;z.loc[T,'candidate']=y;z.loc[T,'source']='native_train';z.loc[T,'reason']='native_fit';z['evaluable']=ok;z['method']='CellCharter';z['reference']='both';z['radius']=r;z['representation']=rep;z.to_parquet(o/'assignments.parquet',index=False);np.save(o/'train_labels.npy',y);f.save(o/'model');(o/'buffer_validation.json').write_text(json.dumps(buf,indent=2));(o/'complete.json').write_text(json.dumps(dict(method='CellCharter',reference='both',radius=r,representation=rep,selected_parameter=k,n_clusters=len(set(y)),n_train=len(T),evaluable=int(ok.sum()),unresolved=int((z.label<0).sum()),distance_limit=threshold,converged=bool(f.converged_),**buf),indent=2));print('CELLCHARTER_COMPLETE',r,rep,flush=True)
if __name__=='__main__':
 gate();aggregate()
 for rep in ['scVI','scVI_without_Top20']:
  for r in [50,200]:fit(r,rep)
