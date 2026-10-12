from base import *
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score,adjusted_rand_score,normalized_mutual_info_score
from scipy.spatial.distance import pdist,squareform
from scipy.optimize import linear_sum_assignment
import igraph as ig,leidenalg
import subprocess

def km(z,k,seed=123):
 model=KMeans(n_clusters=int(k),n_init=50,max_iter=300,algorithm='lloyd',random_state=seed).fit(z)
 if model.n_iter_>=300:model=KMeans(n_clusters=int(k),n_init=1,init=model.cluster_centers_,max_iter=1000,algorithm='lloyd',random_state=seed).fit(z)
 return model

def graph(z):
 n=len(z);ids=NearestNeighbors(n_neighbors=30,n_jobs=1).fit(z).kneighbors(z,return_distance=False)
 a=sparse.csr_matrix((np.ones(n*30,np.float32),(np.repeat(np.arange(n),30),ids.ravel())),shape=(n,n));inter=(a@a.T).tocoo();w=inter.data/(60-inter.data);ok=(inter.row<inter.col)&(w>1/15)
 g=ig.Graph(n=n,edges=list(zip(inter.row[ok].tolist(),inter.col[ok].tolist())),directed=False);g.es['weight']=w[ok].astype(float).tolist();return g

def leiden(g,r,seed=123):return np.asarray(leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,weights='weight',resolution_parameter=float(r),seed=seed,n_iterations=10).membership,dtype=np.int32)

def neighbors(z,ref):
 nn=NearestNeighbors(n_neighbors=31,n_jobs=1).fit(z[ref]);ds=[];ids=[];local=np.full(len(z),-1,int);local[ref]=np.arange(len(ref))
 for start in range(0,len(z),5000):
  d,i=nn.kneighbors(z[start:start+5000]);keep=i!=local[start:start+len(d),None];order=np.argsort(~keep,axis=1,kind='stable')[:,:30];ds.append(np.take_along_axis(d,order,axis=1));ids.append(np.take_along_axis(i,order,axis=1))
 return np.concatenate(ds),np.concatenate(ids)

def votes(lab,ids,ds,candidate=None):
 ls=lab[ids[:,:15]];dist=ds[:,:15];k=int(lab.max())+1;counts=np.stack([(ls==c).sum(1) for c in range(k)],1);best=counts.max(1);cost=np.stack([np.where(counts[:,c]>0,(dist*(ls==c)).sum(1)/np.maximum(counts[:,c],1),np.inf) for c in range(k)],1);cost[counts!=best[:,None]]=np.inf;pred=cost.argmin(1)
 if candidate is None:candidate=pred
 support=(ls==candidate[:,None]).mean(1);return pred,support

def profiles(out,lab,m,f):
 rows=[];by=[];sup=[];q=np.concatenate([f['q1'],f['q2']],1);glob=q.mean(0)
 for c in np.unique(lab):
  ix=np.flatnonzero(lab==c);p=q[ix].mean(0);mm=m.iloc[ix];ds=mm.dataset_id.value_counts(normalize=True);ss=mm.sample_uid.value_counts(normalize=True)
  sup.append(dict(cluster=int(c),n=len(ix),fraction=len(ix)/len(m),n_dataset=len(ds),n_sample=len(ss),n_source_known_patients=mm.loc[mm.source_patient_identity_known,'patient_uid'].nunique(),n_analysis_units=mm.analysis_unit_uid.nunique(),max_dataset_fraction=ds.max(),max_sample_fraction=ss.max(),dataset_dominated=ds.max()>.8,sample_dominated=ss.max()>.5,small=len(ix)<.005*len(m),dispersion=float(np.mean(np.sum((q[ix]-p)**2,1)))))
  for j,v in enumerate(p):rows.append(dict(cluster=int(c),ring=['R1','R2'][j//30],cell_type=TYPES[j%30],mean_relative_composition=float(v),global_mean=float(glob[j]),log2_enrichment=float(np.log2((v+1e-6)/(glob[j]+1e-6)))))
  for field in ['dataset_id','sample_uid']:
   for unit,gg in mm.groupby(field):
    pp=q[gg.index].mean(0);disp=float(np.mean(np.sum((q[gg.index]-pp)**2,1)))
    for j,v in enumerate(pp):by.append(dict(cluster=int(c),level=field,unit=unit,n=len(gg),ring=['R1','R2'][j//30],cell_type=TYPES[j%30],mean_relative_composition=float(v),dispersion=disp))
 tab(out/'profiles_R1_R2.tsv',rows);tab(out/'profile_by_dataset_sample.tsv',by);tab(out/'cluster_support.tsv',sup)
 return pd.DataFrame(sup)

def integration(out,z,m,f,ref,ids,g):
 mr=m.iloc[ref];ri=ids[ref];raw=f['X12'].astype(float);edges=np.asarray(g.get_edgelist());w=np.asarray(g.es['weight']);rows=[]
 for field in ['dataset_id','sample_uid']:
  a=mr[field].to_numpy();same=a[ri]==a[:,None];ent=[]
  for r in a[ri]:
   _,counts=np.unique(r,return_counts=True);p=counts/len(r);ent.append(float(-(p*np.log(p)).sum()))
  es=a[edges[:,0]]==a[edges[:,1]];hd=np.linalg.norm(raw[ref,None,:]-raw[ref[ri]],axis=2)
  rng=np.random.default_rng(123);rnd=rng.integers(0,len(ref),size=ri.shape);rhd=np.linalg.norm(raw[ref,None,:]-raw[ref[rnd]],axis=2)
  rows.append(dict(source_key=field,n_eval=len(ref),k=30,same_source_neighbor_fraction=same.mean(),neighbor_entropy=np.mean(ent),same_source_edge_fraction=es.mean(),same_source_edge_weight_fraction=w[es].sum()/w.sum(),raw_Hellinger_neighbor_mean=hd.mean(),raw_Hellinger_within=hd[same].mean(),raw_Hellinger_cross=hd[~same].mean() if (~same).any() else np.nan,random_pair_Hellinger=rhd.mean(),components=len(g.connected_components()),isolates=sum(np.asarray(g.degree())==0)))
 # Composition strata are frozen dominant major-like subtype from uncorrected two-ring q, descriptive only.
 strata=(f['q1'][ref]+f['q2'][ref]).argmax(1);cond=[]
 for c in np.unique(strata):
  ix=np.flatnonzero(strata==c)
  for field in ['dataset_id','sample_uid']:
   a=mr[field].to_numpy();cond.append(dict(stratum=TYPES[c],source_key=field,n=len(ix),same_source_neighbor_fraction=float((a[ri[ix]]==a[ix,None]).mean()),raw_Hellinger_neighbor_mean=float(np.linalg.norm(raw[ref[ix],None,:]-raw[ref[ri[ix]]],axis=2).mean())))
 tab(out/'conditional_composition_metrics.tsv',cond);tab(out/'integration_metrics.tsv',rows)
 # One representative per non-overlapping 8x8 array-coordinate block in each section.
 block=mr.sample_uid+'::'+(mr.array_row//8).astype(str)+'::'+(mr.array_col//8).astype(str);sub=np.flatnonzero(~block.duplicated());np.save(out/'spatial_block_reference_indices.npy',sub)
 if len(sub)>31:
  jj=NearestNeighbors(n_neighbors=31,n_jobs=1).fit(z[ref[sub]]).kneighbors(z[ref[sub]],return_distance=False)[:,1:];br=[]
  for field in ['dataset_id','sample_uid']:
   a=mr.iloc[sub][field].to_numpy();br.append(dict(source_key=field,n=len(sub),k=30,same_source_fraction=float((a[jj]==a[:,None]).mean()),definition='first frozen reference anchor per section 8x8 array block; descriptive sparse diagnostic'))
  tab(out/'spatial_block_mixing.tsv',br)
 return rows

def _run_rep(ring,mode,strength):
 rid=f'{mode}_{strength}_{ring}' if mode!='NONE' else f'NONE_{ring}';out=OUT/'02_R12_pilot'/rid if ring=='R12' else OUT/'01_features'/rid;out.mkdir(exist_ok=True)
 if (out/'COMPLETE.json').exists():return out
 m,f,ref=inputs();m=pd.read_parquet(OUT/'01_features/metadata.parquet');z=np.load(OUT/'01_features'/f'{ring}_PCA.npy')
 theta,lam=(1,5) if strength=='WEAK' else (2,1)
 if mode!='NONE':
  if not (out/'embedding.npy').exists():
   if not (out/'harmony_manifest.json').exists():
    log('Harmony '+rid)
    with open(out/'harmony.log','w') as logf:subprocess.run(['/usr/local/bin/Rscript',str(OUT/'code/harmony.R'),str(OUT/'01_features'/f'{ring}_PCA.tsv'),str(OUT/'01_features/harmony_metadata.tsv'),'dataset_id' if mode=='DATASET' else 'sample_uid',str(theta),str(lam),str(out)],stdout=logf,stderr=subprocess.STDOUT,check=True)
   zz=np.loadtxt(out/'embedding.tsv',delimiter='\t');assert zz.shape==z.shape and np.isfinite(zz).all();np.save(out/'embedding.npy',zz)
  z=np.load(out/'embedding.npy')
 else:np.save(out/'embedding.npy',z)
 zr=z[ref];evaluid=pd.read_csv(OLD/'05_comparison/silhouette_evaluation_UIDs.tsv',sep='\t').spot_uid;ev=np.flatnonzero(m.iloc[ref].spot_uid.isin(evaluid));D=squareform(pdist(zr[ev]));rawD=squareform(pdist(f[{'R1':'X1','R2':'X2','R12':'X12'}[ring]][ref[ev]].astype(float)))
 if (out/'neighbors.npz').exists():a=np.load(out/'neighbors.npz');ds=a['distance'];ids=a['ids']
 else:log('Neighbors '+rid);ds,ids=neighbors(z,ref);np.savez_compressed(out/'neighbors.npz',distance=ds,ids=ids)
 if (out/'graph.pickle').exists():g=ig.Graph.Read_Pickle(str(out/'graph.pickle'))
 else:g=graph(zr);g.write_pickle(str(out/'graph.pickle'))
 if not (out/'integration_metrics.tsv').exists():integration(out,z,m,f,ref,ids,g)
 cache=out/'candidates.joblib';models,labs,rows=joblib.load(cache) if cache.exists() else ({},{},[])
 def candidate(alg,p):
  tag=f'{alg}_{p}'
  if tag in labs:return
  log(rid+' '+tag)
  if alg=='KM':models[tag]=km(zr,p);lab=models[tag].labels_;conv=models[tag].n_iter_<models[tag].max_iter
  else:lab=leiden(g,p);conv=True
  labs[tag]=lab;k=len(set(lab));sizes=np.bincount(lab);q=np.concatenate([f['q1'],f['q2']],1);cent=np.array([q[ref][lab==c].mean(0) for c in range(k)]);dist=pdist(cent)
  pred=models[tag].predict(z) if alg=='KM' else votes(lab,ids,ds)[0];pred[ref]=lab;tr,support=votes(lab,ids,ds,pred);threshold=np.quantile(ds[ref,0],.99);accepted=(support>=.6)&(ds[:,0]<=threshold)
  sil=lambda dd:float(silhouette_score(dd,lab[ev],metric='precomputed')) if 1<len(set(lab[ev]))<len(ev) else np.nan
  dom=max(m.iloc[ref].assign(cluster=lab).groupby('cluster').dataset_id.value_counts(normalize=True))
  row=dict(algorithm=alg,parameter=p,k=k,converged=conv,silhouette_embedding=sil(D),silhouette_raw_Hellinger=sil(rawD),min_cluster=int(sizes.min()),tiny_clusters=int((sizes<max(20,.001*len(ref))).sum()),min_profile_distance=float(dist.min()) if len(dist) else np.nan,raw_composition_dispersion=float(np.mean(np.sum((q[ref]-cent[lab])**2,1))),max_dataset_fraction=float(dom),low_support_fraction=float((support<.6).mean()),distance_outlier_fraction=float((ds[:,0]>threshold).mean()),accepted_fraction=float(accepted.mean()))
  rows.append(row);cd=out/'candidates'/tag;cd.mkdir(parents=True,exist_ok=True);np.save(cd/'reference_labels.npy',lab);tab(cd/'reference_sizes.tsv',pd.Series(lab).value_counts().rename_axis('cluster').reset_index(name='n'));tab(cd/'full_sizes.tsv',pd.Series(pred).value_counts().rename_axis('cluster').reset_index(name='n'));tab(cd/'profiles.tsv',[dict(cluster=c,ring=['R1','R2'][j//30],cell_type=TYPES[j%30],mean=float(v)) for c,pr in enumerate(cent) for j,v in enumerate(pr)]);tab(cd/'source_support.tsv',m.iloc[ref].assign(cluster=lab).groupby(['cluster','dataset_id','sample_uid']).size().reset_index(name='n'));joblib.dump((models,labs,rows),cache);tab(out/'candidate_parameters.tsv',rows)
 for p in [2,4,6,8,10,12]:candidate('KM',p)
 for p in [.005,.01,.02,.05,.1]:candidate('LEIDEN',p)
 lc=[r['k'] for r in rows if r['algorithm']=='LEIDEN'];extra=[]
 if max(lc)<6:extra=[.2,.4]
 elif min(lc)>12:extra=[.0005,.001,.002]
 js(out/'resolution_extension.json',dict(initial_K=lc,extra=extra,trigger='all below 6' if max(lc)<6 else 'all above 12' if min(lc)>12 else 'none'))
 for p in extra:candidate('LEIDEN',p)
 df=pd.DataFrame(rows)
 # Predefined conservative granularity: K6 unless extra composition separation and support justify increasing.
 valid=df[(df.algorithm=='KM')&df.converged&(df.k>=6)&(df.tiny_clusters==0)].sort_values('k');best=valid.iloc[0]
 for _,r in valid.iloc[1:].iterrows():
  gain=1-r.raw_composition_dispersion/best.raw_composition_dispersion
  if gain>=.15 and r.min_profile_distance>=.10 and r.accepted_fraction>=best.accepted_fraction-.03 and r.silhouette_embedding>=best.silhouette_embedding:best=r
 choices={'KM':best};lv=df[(df.algorithm=='LEIDEN')&df.converged&(df.k>1)&(df.tiny_clusters==0)];separated=lv[(lv.silhouette_embedding>0)&(lv.silhouette_raw_Hellinger>0)];lv=separated if len(separated) else lv;moderate=lv[lv.k.between(6,12)];lv=moderate if len(moderate) else lv
 if not len(lv):lv=df[(df.algorithm=='LEIDEN')&(df.k>1)]
 choices['LEIDEN']=lv.sort_values(['silhouette_embedding','accepted_fraction'],ascending=False).iloc[0]
 for alg,best in choices.items():
  od=out/alg;od.mkdir(exist_ok=True)
  if (od/'working_config.json').exists():continue
  param=int(best.parameter) if alg=='KM' else float(best.parameter);tag=f'{alg}_{param}';lab=labs[tag];st=[]
  for seed in range(124,129):
   uid=pd.read_csv(OLD/f'05_comparison/stability_UIDs_{seed}.tsv',sep='\t').spot_uid;sub=np.flatnonzero(m.iloc[ref].spot_uid.isin(uid));sp=od/f'stability_{seed}.npz'
   if sp.exists():sl=np.load(sp)['labels']
   else:
    log(rid+' '+alg+' stability '+str(seed));sl=km(zr[sub],param,seed).labels_ if alg=='KM' else leiden(graph(zr[sub]),param,seed);np.savez_compressed(sp,indices=sub,labels=sl)
   cont=pd.crosstab(lab[sub],sl).to_numpy();aa,bb=linear_sum_assignment(-cont);q=np.concatenate([f['q1'],f['q2']],1)[ref[sub]];c1=np.array([q[lab[sub]==c].mean(0) for c in sorted(set(lab[sub]))]);c2=np.array([q[sl==c].mean(0) for c in sorted(set(sl))]);st.append(dict(seed=seed,n=len(sub),k=len(set(sl)),ARI=adjusted_rand_score(lab[sub],sl),matched_fraction=cont[aa,bb].sum()/len(sub),matched_profile_L2=float(np.linalg.norm(c1[aa]-c2[bb],axis=1).mean()),stability_type='conditional_on_fixed_representation'))
  tab(od/'stability.tsv',st)
  if ring=='R12':
   # Remove a deterministic 20% of complete 8x8 section blocks, refit cluster only.
   mr=m.iloc[ref];block=mr.sample_uid+'::'+(mr.array_row//8).astype(str)+'::'+(mr.array_col//8).astype(str);uu=np.sort(block.unique());rng=np.random.default_rng(321);keep=rng.choice(uu,int(.8*len(uu)),replace=False);sub=np.flatnonzero(block.isin(keep));sp=od/'block_stability.npz'
   if sp.exists():sl=np.load(sp)['labels']
   else:sl=km(zr[sub],param,321).labels_ if alg=='KM' else leiden(graph(zr[sub]),param,321);np.savez_compressed(sp,indices=sub,labels=sl)
   tab(od/'block_stability.tsv',[dict(n=len(sub),ARI=adjusted_rand_score(lab[sub],sl),k=len(set(sl)),stability_type='section_spatial_8x8_block_subsample_fixed_representation')])
  candidate_lab=models[tag].predict(z) if alg=='KM' else votes(lab,ids,ds)[0];candidate_lab[ref]=lab;tr,support=votes(lab,ids,ds,candidate_lab);threshold=float(np.quantile(ds[ref,0],.99));low=support<.6;far=ds[:,0]>threshold;accepted=np.where(low|far,-1,candidate_lab);native=np.full(len(z),-1,int);native[ref]=lab
  labdf=m[['spot_uid','dataset_id','section_id','sample_uid','patient_uid','analysis_unit_uid','reference','complete_rings']].copy();labdf['reference_native_label']=native;labdf['candidate_label']=candidate_lab;labdf['common_transfer_label']=tr;labdf['accepted_label']=accepted;labdf['candidate_vote_fraction']=support;labdf['nearest_reference_distance']=ds[:,0];labdf['low_support']=low;labdf['distance_outlier']=far;labdf['label_origin']=np.where(m.reference,'reference_native','nearest_centroid' if alg=='KM' else '15NN_vote');labdf['rejection_reason']=np.select([low&far,low,far],['low_support_and_distance','low_support','distance_outlier'],default='');labdf['assignment_status']=np.where(accepted<0,'UNASSIGNED','ASSIGNED');labdf.to_parquet(od/'labels.parquet',index=False);np.save(od/'reference_native_labels.npy',lab);joblib.dump(models[tag] if alg=='KM' else dict(reference_labels=lab,representation=rid),od/'working_model.joblib')
  sup=profiles(od,candidate_lab,m,f);coverage=labdf.groupby(['dataset_id','sample_uid']).agg(n=('spot_uid','size'),low_support=('low_support','mean'),distance_outlier=('distance_outlier','mean'),accepted_fraction=('accepted_label',lambda v:(v>=0).mean())).reset_index();tab(od/'support_and_assignment.tsv',coverage)
  conv=json.loads((out/'harmony_manifest.json').read_text())['convergence_status'] if mode!='NONE' else 'not_applicable';config=dict(config_id=f'{mode}_{ring}_{alg}',integration_mode=mode,batch_key='none' if mode=='NONE' else 'dataset_id' if mode=='DATASET' else 'sample_uid',strength=strength,theta=None if mode=='NONE' else theta,**{'lambda':None if mode=='NONE' else lam},ring_input=ring,representation_id=rid,n_pc=z.shape[1],algorithm=alg,k_or_resolution=param,n_clusters=int(best.k),n_reference=len(ref),n_all=len(m),fit_status='COMPLETE',convergence_status=conv,stability_type='conditional_on_fixed_representation',stability_ARI_median=float(pd.DataFrame(st).ARI.median()),low_support_fraction=float(low.mean()),distance_outlier_fraction=float(far.mean()),accepted_fraction=float((accepted>=0).mean()),dataset_dominance_summary=int(sup.dataset_dominated.sum()),sample_dominance_summary=int(sup.sample_dominated.sum()),composition_preservation_summary=float(best.raw_composition_dispersion),recommendation_status='provisional_pending_integration_comparison',silhouette_embedding=float(best.silhouette_embedding),silhouette_raw_Hellinger=float(best.silhouette_raw_Hellinger),distance_Q99=threshold,selection_rule='KM starts K6; increase only >=15% dispersion gain, min profile distance>=.10, no silhouette loss, coverage loss<=.03; Leiden requires positive embedding and raw silhouette when available; then non-tiny K6-12 if supported, else finite-grid best silhouette; no forced K range')
  js(od/'working_config.json',config);tab(od/'candidate_parameters.tsv',df[df.algorithm==alg]);js(od/'graph_or_model_manifest.json',dict(graph='exact 30NN including self, binary SNN Jaccard>1/15, undirected no selfloops',components=len(g.connected_components()),isolates=sum(np.asarray(g.degree())==0),seed=123,n_init=50,max_iter=300,leiden_n_iterations=10,partition='RBConfigurationVertexPartition',representation=rid))
 k=int(choices['LEIDEN'].k)
 if 2<=k<=20:
  kk=f'KM_{k}';mod=models.get(kk)
  if mod is None:mod=km(zr,k)
  ll=labs[f"LEIDEN_{float(choices['LEIDEN'].parameter)}"];tab(out/'sameK_diagnostics.tsv',[dict(representation=rid,k=k,ARI=adjusted_rand_score(mod.labels_,ll),NMI=normalized_mutual_info_score(mod.labels_,ll))]);joblib.dump(mod,out/'sameK_KM.joblib')
 if not 2<=k<=20:tab(out/'sameK_diagnostics.tsv',[dict(representation=rid,k=k,ARI=np.nan,NMI=np.nan,status='not_run_selected_K_outside_2_20')])
 js(out/'COMPLETE.json',dict(status='complete',mp10_values_loaded=False));return out

def run_rep(ring,mode,strength):
 import fcntl
 lock=open(OUT/'logs'/f'{mode}_{strength}_{ring}.lock','w')
 fcntl.flock(lock,fcntl.LOCK_EX)
 try:return _run_rep(ring,mode,strength)
 finally:fcntl.flock(lock,fcntl.LOCK_UN);lock.close()

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('ring');p.add_argument('mode');p.add_argument('strength',default='NONE',nargs='?');a=p.parse_args();run_rep(a.ring,a.mode,a.strength)
