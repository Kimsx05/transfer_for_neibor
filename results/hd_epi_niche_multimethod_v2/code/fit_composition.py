from common import *
from graph_core import graph,partition,diagnostics
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score as ari,silhouette_score
from sklearn.neighbors import NearestNeighbors
import joblib,gc

def km(x,k,seed):return MiniBatchKMeans(n_clusters=k,random_state=seed,n_init=10,batch_size=4096,max_iter=200,tol=1e-4,max_no_improvement=20).fit(x)
def sil(x,y,ix):
 return float(silhouette_score(x[ix],y[ix])) if 1<len(np.unique(y[ix]))<len(ix) else np.nan

def select(m):
 # A transparent parsimony band, not decimal-score optimization across feature spaces.
 good=m[(m.n_clusters>=2)&(m.seed_ARI>=.8)]
 if len(good)==0:good=m[m.n_clusters>=2];reason='exploratory_no_seed_stability_ge_0.8'
 else:reason='stable_parsimony_band'
 if len(good)==0:return m.iloc[0], 'no_nontrivial_partition'
 band=good[good.silhouette>=good.silhouette.max()-.05]
 if not len(band):band=good
 return band.sort_values(['n_clusters','parameter']).iloc[0],reason

def transfer(x,y,q):
 nn=NearestNeighbors(n_neighbors=min(15,len(x)),algorithm='kd_tree',n_jobs=4).fit(x);u,inv=np.unique(x,axis=0,return_inverse=True);exact={};codes=np.unique(y)
 hist=np.zeros((len(u),int(y.max())+1),np.int32);np.add.at(hist,(inv,y),1)
 for v,c in zip(u,hist):exact[v.tobytes()]=(int(c.argmax()),float(c.max()/c.sum()))
 out=np.empty(len(q),int);vote=np.empty(len(q));dist=np.empty(len(q))
 for st in range(0,len(q),10000):
  dd,jj=nn.kneighbors(q[st:st+10000]);ys=y[jj]
  for j,(row,d,z) in enumerate(zip(q[st:st+10000],dd,ys)):
   pair=exact.get(row.tobytes())
   if pair is None:c=np.bincount(z,minlength=hist.shape[1]);pair=(int(c.argmax()),float(c.max()/len(z)))
   out[st+j],vote[st+j]=pair;dist[st+j]=d[0]
 return out,vote,dist

def run(x,ok,method,ref,r,rep,preprocessor=None):
 a=anchors();T=train_index(a);assert ok[T].all();S=np.flatnonzero(ok);xt=np.asarray(x[T],np.float64);tag=f'{method}_{ref}_r{r}_{rep}';o=R/'models'/tag;o.mkdir(exist_ok=True)
 if (o/'complete.json').exists():print('CACHE',tag,flush=True);return
 ix=np.random.default_rng(123).choice(len(T),2500,replace=False);metrics=[];fits={};labels={};diags=[]
 if method=='hellinger_kmeans':
  for k in CFG['K']:
   f=km(xt,k,123);y=f.predict(xt);st=[ari(y,km(xt,k,seed).predict(xt)) for seed in [456,789]];metrics.append(dict(parameter=k,n_clusters=len(set(y)),seed_ARI=np.mean(st),seed_ARI_456=st[0],seed_ARI_789=st[1],silhouette=sil(xt,y,ix),converged_iterations=int(f.n_iter_)));fits[k]=f;labels[k]=y;print("CANDIDATE",tag,k,metrics[-1],flush=True)
 else:
  g,u,inv,n=graph(xt,15)
  for res in CFG['resolutions']:
   y=partition(g,inv,res,123);st=[ari(y,partition(g,inv,res,seed)) for seed in [456,789]];metrics.append(dict(parameter=res,n_clusters=len(set(y)),seed_ARI=np.mean(st),seed_ARI_456=st[0],seed_ARI_789=st[1],silhouette=sil(xt,y,ix),converged_iterations=-1));labels[res]=y;print("CANDIDATE",tag,res,metrics[-1],flush=True)
  d=pd.DataFrame(metrics);chosen,_=select(d);rr=float(chosen.parameter)
  for k in [15,30,50] if method=='adapted_leiden' else [15]:
   gg,uu,ii,nn=(g,u,inv,n) if k==15 else graph(xt,k);yy=labels[rr] if k==15 else partition(gg,ii,rr,123);z=diagnostics(gg,uu,ii,nn,yy);z.update(k=k,ARI_to_k15=ari(labels[rr],yy),resolution=rr);diags.append(z)
  # Fixed primary k15; k30/50 are finite construction sensitivity, not an unrestricted grid.
  pd.DataFrame(diags).to_csv(o/'graph_diagnostics.tsv',sep='\t',index=False)
 m=pd.DataFrame(metrics);chosen,why=select(m);p=float(chosen.parameter);p=int(p) if method=='hellinger_kmeans' else p;y=labels[p];m['selected']=m.parameter==p;m.to_csv(o/'candidate_metrics.tsv',sep='\t',index=False)
 keep,test,distbuffer=buffer_split(a,T,r);buf={'retained_train':int(keep.sum()),'validation':int(test.sum()),'minimum_cross_distance_um':float(distbuffer[keep].min()),'required_gt_um':2*r,'scope':'fixed-representation buffered validation; not independent test'}
 if method=='hellinger_kmeans':
  fit=fits[p];ds=fit.transform(x[S]);pred=fit.predict(x[S]);dd=np.sort(ds,axis=1);vote=(dd[:,1]-dd[:,0])/np.maximum(dd[:,1],1e-12);dist=dd[:,0];threshold=float(np.quantile(fit.transform(xt).min(1),.995)*1.5+1e-8);valid=dist<=threshold
  bf=km(xt[keep],p,456);buf['ARI_to_full_on_buffer_test']=float(ari(y[test],bf.predict(xt[test])));joblib.dump(fit,o/'model.joblib')
  coarse=fits[2];coarse_pred=coarse.predict(x[S]);np.savez_compressed(o/'K2.npz',S=S,labels=coarse_pred);joblib.dump(coarse,o/'K2_model.joblib')
 else:
  pred,vote,dist=transfer(xt,y,np.asarray(x[S],np.float64));nn=NearestNeighbors(n_neighbors=1,n_jobs=4).fit(xt[keep]);hd=nn.kneighbors(xt[test])[0][:,0];threshold=float(np.quantile(hd,.995)*1.5+1e-8);valid=(dist<=threshold)&(vote>.5)
  bg,bu,bi,bn=graph(xt[keep],15);by=partition(bg,bi,p,456);bp,bv,bd=transfer(xt[keep],by,xt[test]);buf['ARI_to_full_on_buffer_test']=float(ari(y[test],bp));np.savez_compressed(o/'training.npz',X=xt,labels=y,T=T)
 z=a[['section','bin_id']].copy();z['candidate']=-1;z['label']=-1;z['source']='unresolved';z['reason']='n_known_below_10';z['confidence']=np.nan;z['distance']=np.nan;z.loc[S,'candidate']=pred;z.loc[S,'confidence']=vote;z.loc[S,'distance']=dist;z.loc[S,'reason']='outside_training_support_or_vote_ambiguity';z.loc[S[valid],'label']=pred[valid];z.loc[S[valid],'source']='model_predict' if method=='hellinger_kmeans' else 'external_projection';z.loc[S[valid],'reason']='assigned';z.loc[T,'label']=y;z.loc[T,'candidate']=y;z.loc[T,'source']='native_train';z.loc[T,'reason']='native_fit';z['method']=method;z['reference']=ref;z['radius']=r;z['representation']=rep;z['evaluable']=ok
 z.to_parquet(o/'assignments.parquet',index=False);np.save(o/'train_labels.npy',y);joblib.dump(preprocessor,o/'preprocessor.joblib');(o/'buffer_validation.json').write_text(json.dumps(buf,indent=2));(o/'complete.json').write_text(json.dumps(dict(method=method,reference=ref,radius=r,representation=rep,selected_parameter=p,n_clusters=int(chosen.n_clusters),selection_reason=why,n_train=len(T),evaluable=int(ok.sum()),unresolved=int((z.label<0).sum()),distance_limit=threshold,seed_ARI=float(chosen.seed_ARI),silhouette=float(chosen.silhouette),**buf),indent=2));print('COMPLETE',tag,flush=True)

def main(ref_filter=None,radius_filter=None):
 gate();a=anchors();T=train_index(a)
 for ref in ([ref_filter] if ref_filter else TYPES):
  for r in ([radius_filter] if radius_filter else [50,200]):
   d,P,ok=composition(ref,r);H=np.sqrt(np.nan_to_num(P,nan=0));pca=PCA(svd_solver='full').fit(H[T]);rank=int((pca.singular_values_>pca.singular_values_[0]*1e-10).sum());x=pca.transform(H)[:,:rank];pre={'PCA':pca,'rank':rank,'transform':'sqrt(P); centering; all nonzero PCs'}
   assert np.allclose(np.linalg.norm(H[T[:100]]-H[T[0]],axis=1),np.linalg.norm(x[T[:100]]-x[T[0]],axis=1))
   for method in ['hellinger_kmeans','adapted_leiden']:run(x,ok,method,ref,r,'NONE',pre)
   import harmonypy
   hp=R/'models'/f'harmony_{ref}_r{r}.npz'
   if hp.exists():hx=np.load(hp)['X']
   else:
    ho=harmonypy.run_harmony(x[ok].T,pd.DataFrame({'section':a.loc[ok,'section'].to_numpy()}),'section',theta=1,lamb=5,nclust=50,max_iter_harmony=10,random_state=123,ncores=4);z=np.asarray(ho.Z_corr);z=z.T if z.shape[1]==ok.sum() else z;assert z.shape==x[ok].shape;hx=np.full(x.shape,np.nan);hx[ok]=z;np.savez_compressed(hp,X=hx)
    (R/'tables'/f'harmony_{ref}_r{r}.json').write_text(json.dumps(dict(scope='all_five_sections_transductive',theta=1,lamb=5,rank=rank,MPBC2_in_representation=True,MPBC2_in_niche_fit=False,objective=[float(v) for v in getattr(ho,'objective_harmony',[])]),indent=2))
   for method in ['hellinger_kmeans','adapted_leiden']:run(hx,ok,method,ref,r,'WEAK',pre)
 for ref in ([ref_filter] if ref_filter else TYPES):
  for r in ([radius_filter] if radius_filter else [50,200]):
   d,P,ok=composition(ref,r);c=[]
   for s in SECTIONS:
    z=np.load(R/'models'/f'{s}_rings.npz');ix=[list(z['columns']).index(ref+'::'+t) for t in TYPES[ref]];c.append(z['counts'][:,:r//25,ix].reshape(len(z['bin_id']),-1))
   c=np.concatenate(c).astype(float);target=float(np.median(c[T].sum(1)));v=np.log1p(np.divide(c*target,c.sum(1,keepdims=True),out=np.zeros_like(c),where=c.sum(1,keepdims=True)>0));pca=PCA(n_components=min(50,v.shape[1]-1),svd_solver='full').fit(v[T]);x=pca.transform(v);run(x,ok,'MENDER_adapted',ref,r,'ring_counts',{'PCA':pca,'target_total':target,'normalization':'native median T total then log1p','feature_order':'ring-major known types','include_self':False});del c,v,x;gc.collect()
if __name__=='__main__':main()
