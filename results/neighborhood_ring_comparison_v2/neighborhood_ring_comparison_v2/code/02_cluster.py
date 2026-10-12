from common import *
import joblib,itertools,platform,importlib.metadata
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score,adjusted_rand_score
from scipy.spatial.distance import pdist,squareform
from scipy.optimize import linear_sum_assignment
import igraph as ig,leidenalg
import matplotlib.pyplot as plt
KS=[2,4,6,8,10,12,16];RES=[.1,.2,.4,.6,.8,1.];SEEDS=[124,125,126,127,128]

def fitkm(x,k,seed=123):return KMeans(n_clusters=int(k),n_init=50,max_iter=300,random_state=seed,algorithm='lloyd').fit(x)
def pca95(x):
 p=PCA(n_components=.95,svd_solver='full',whiten=False);z=p.fit_transform(x);return p,z

def graph(z):
 n=len(z);k=min(30,n);nn=NearestNeighbors(n_neighbors=k,metric='euclidean',n_jobs=1).fit(z);ids=nn.kneighbors(z,return_distance=False)
 a=sparse.csr_matrix((np.ones(n*k,np.float32),(np.repeat(np.arange(n),k),ids.ravel())),shape=(n,n))
 inter=(a@a.T).tocoo();w=inter.data/(2*k-inter.data);ok=(inter.row<inter.col)&(w>1/15)
 g=ig.Graph(n=n,edges=list(zip(inter.row[ok].tolist(),inter.col[ok].tolist())),directed=False);g.es['weight']=w[ok].astype(float).tolist();return g

def leiden(g,r,seed=123):return np.asarray(leidenalg.find_partition(g,leidenalg.RBConfigurationVertexPartition,weights='weight',resolution_parameter=float(r),seed=seed,n_iterations=10).membership,dtype=np.int32)
def transfer(zref,z,labels,ref):
 n=len(z);nn=NearestNeighbors(n_neighbors=16,n_jobs=1).fit(zref);pred=np.empty(n,int);vote=np.empty(n);nearest=np.empty(n)
 local={int(v):i for i,v in enumerate(ref)}
 for start in range(0,n,4000):
  d,ids=nn.kneighbors(z[start:start+4000]);
  for a in range(len(d)):
   i=start+a;keep=ids[a]!=local[i] if i in local else np.ones(16,bool);ni=ids[a][keep][:15];di=d[a][keep][:15];ls=labels[ni];u,counts=np.unique(ls,return_counts=True);best=u[counts==counts.max()];chosen=min(best,key=lambda c:(float(di[ls==c].mean()),int(c)));pred[i]=chosen;vote[i]=counts.max()/len(ni);nearest[i]=di[0]
 return pred,vote,nearest

def groups(meta,labels):
 t=meta[['dataset_id','patient_uid','section_id','source_patient_identity_known']].copy();t['cluster']=labels;out=[]
 for c,g in t.groupby('cluster'):
  known=g[g.source_patient_identity_known & g.patient_uid.notna()];d=g.dataset_id.value_counts(normalize=True);p=known.patient_uid.value_counts(normalize=True)
  out.append(dict(cluster=int(c),n=len(g),n_dataset=g.dataset_id.nunique(),n_source_patients=known.patient_uid.nunique(),n_sections=g.section_id.nunique(),max_dataset_fraction=d.max(),max_known_patient_fraction=p.max() if len(p) else np.nan))
 return out

def run_input(name,key,meta,features,ref,evalix,subs):
 out=OUT/'04_clustering'/name;out.mkdir(exist_ok=True)
 if (out/'COMPLETE.json').exists():log(name+' checkpoint complete');return
 x=features[key].astype(float);xr=x[ref];q=np.concatenate([features['q1'],features['q2']],axis=1);candidate=[];models={};labels={}
 cache=out/'candidates.joblib'
 if cache.exists():models,labels,candidate=joblib.load(cache)
 p,zr=pca95(xr);z=p.transform(x);joblib.dump(p,out/'pca95.joblib');np.save(out/'pca_projection.npy',z.astype('float32'));js(out/'pca.json',dict(n_pc=p.n_components_,variance_ratio=p.explained_variance_ratio_.tolist(),center=True,scale=False,whiten=False))
 D=squareform(pdist(xr[evalix]));g=graph(zr)
 for alg,params in [('KM',KS),('LEIDEN',RES)]:
  for param in params:
   tag=f'{alg}_{param}'
   if tag in labels:continue
   log(f'{name} {tag}')
   if alg=='KM':model=fitkm(xr,param);lab=model.labels_;models[tag]=model;converged=model.n_iter_<300
   else:lab=leiden(g,param);converged=True
   labels[tag]=lab;k=len(np.unique(lab));sc=silhouette_score(D,lab[evalix],metric='precomputed') if 1<len(np.unique(lab[evalix]))<len(evalix) else np.nan
   sup=groups(meta.iloc[ref],lab);sizes=np.bincount(lab);cent=np.array([q[ref][lab==c].mean(0) for c in sorted(np.unique(lab))]);dist=pdist(cent)
   candidate.append(dict(input=name,algorithm=alg,parameter=param,k=k,silhouette=sc,min_reference_cluster=int(sizes.min()),tiny_cluster=bool(sizes.min()<max(20,.001*len(ref))),max_dataset_dependency=max(t['max_dataset_fraction'] for t in sup),max_known_patient_dependency=max(t['max_known_patient_fraction'] for t in sup),min_profile_euclidean=float(dist.min()) if len(dist) else np.nan,profile_redundancy_flag=bool(len(dist) and dist.min()<.03),converged=converged,n_pc=int(p.n_components_) if alg=='LEIDEN' else 0))
   joblib.dump((models,labels,candidate),cache);tab(out/'candidates.tsv',candidate)
 cand=pd.DataFrame(candidate);summaries=[]
 for alg in ['KM','LEIDEN']:
  candidates=cand[(cand.algorithm==alg)&(cand.k>1)&np.isfinite(cand.silhouette)&cand.converged];assert len(candidates)
  near=candidates[candidates.silhouette>=candidates.silhouette.max()-.02];best=near.sort_values(['k','parameter']).iloc[0];param=int(best.parameter) if alg=='KM' else float(best.parameter);tag=f'{alg}_{param}';lab=labels[tag];config=name+'_'+alg
  od=out/alg;od.mkdir(exist_ok=True);np.save(od/'reference_native_labels.npy',lab)
  stability=[]
  for seed,sub in zip(SEEDS,subs):
   sp=od/f'stability_{seed}.npz'
   if sp.exists():slab=np.load(sp)['labels']
   else:
    log(config+f' stability {seed}')
    if alg=='KM':slab=fitkm(xr[sub],param,seed).labels_
    else:pp,zz=pca95(xr[sub]);slab=leiden(graph(zz),param,seed)
    np.savez_compressed(sp,reference_local_indices=sub,labels=slab)
   ari=adjusted_rand_score(lab[sub],slab);cont=pd.crosstab(lab[sub],slab).to_numpy();aa,bb=linear_sum_assignment(-cont);match=cont[aa,bb].sum()/len(sub)
   stability.append(dict(config=config,seed=seed,n=len(sub),k=len(np.unique(slab)),ARI=ari,matched_fraction=match))
  tab(od/'stability.tsv',stability)
  space=x if alg=='KM' else z;rawtransfer,votes,dist=transfer(space[ref],space,lab,ref);quality_ref=rawtransfer[ref]==lab;common=rawtransfer.copy();common[ref]=lab
  native=models[tag].predict(x) if alg=='KM' else common.copy();native[ref]=lab;threshold=float(np.quantile(dist[ref],.99));
  labelsout=meta[['spot_uid','dataset_id','patient_uid','section_id','reference']].copy();labelsout['cluster_label_native_application']=native;labelsout['cluster_label_common_transfer']=common;labelsout['leave_self_out_transfer_label']=rawtransfer;labelsout['top_vote_fraction']=votes;labelsout['nearest_reference_distance']=dist;labelsout['distance_outlier']=dist>threshold;labelsout['low_vote']=votes<.6;labelsout['label_origin']=np.where(meta.reference,'reference_native','nearest_centroid' if alg=='KM' else 'PCA_15NN_vote');labelsout['centroid_distance']=np.linalg.norm(x-models[tag].cluster_centers_[native],axis=1) if alg=='KM' else np.nan
  labelsout.to_parquet(od/'labels.parquet',index=False);joblib.dump(models[tag] if alg=='KM' else dict(pca=p,reference_coordinates=zr,reference_labels=lab,space='PCA95'),od/'working_model.joblib')
  support=groups(meta,native);tab(od/'cluster_support.tsv',support)
  profiles=[]
  types=json.loads((OUT/'03_features/feature_columns.json').read_text())['q1']
  for c in np.unique(native):
   for ring,arr in [('R1',features['q1']),('R2',features['q2'])]:
    prof=arr[native==c].mean(0)
    for t,v in zip(types,prof):profiles.append(dict(cluster=int(c),ring=ring,cell_type=t,mean_relative_composition=float(v)))
  tab(od/'profiles.tsv',profiles)
  sizes=labelsout.groupby(['dataset_id','patient_uid','section_id','cluster_label_native_application'],dropna=False).size().reset_index(name='n');tab(od/'composition_by_section_patient.tsv',sizes)
  stability=pd.DataFrame(stability);summary=dict(config=config,status='BLIND_FIT_COMPLETE',dimensions=x.shape[1],n_pc=int(p.n_components_) if alg=='LEIDEN' else 0,k=int(best.k),parameter=param,n_reference=len(ref),n_all=len(x),silhouette=float(best.silhouette),stability_median=float(stability.ARI.median()),stability_min=float(stability.ARI.min()),stability_max=float(stability.ARI.max()),reference_leave_self_out_agreement=float(quality_ref.mean()),native_common_agreement=float((native==common).mean()),low_vote_fraction=float((votes<.6).mean()),distance_outlier_fraction=float((dist>threshold).mean()),distance_outlier_threshold=threshold,stable_warning=bool(stability.ARI.min()<.5))
  js(od/'working_config.json',summary);summaries.append(summary)
  # Every candidate full application sizes, saved without retaining redundant full label matrices.
 for tag,lab in labels.items():
  if tag.startswith('KM'):nl=models[tag].predict(x)
  else:
   zz,vv,dd=transfer(zr,z,lab,ref);nl=zz;nl[ref]=lab
  tab(out/(tag+'__full_sizes.tsv'),pd.Series(nl).value_counts().rename_axis('cluster').reset_index(name='n'))
 # same K diagnostic, maximum one additional K per input.
 kl=int([r['k'] for r in summaries if r['config'].endswith('LEIDEN')][0]);lpar=float([r['parameter'] for r in summaries if r['config'].endswith('LEIDEN')][0]);ll=labels[f'LEIDEN_{lpar}'];diagnostics=[]
 if 2<=kl<=20:
  kt=f'KM_{kl}';km=models[kt] if kt in models else fitkm(xr,kl);joblib.dump(km,out/'sameK_KM.joblib');np.save(out/'sameK_KM_labels.npy',km.labels_);diagnostics.append(dict(input=name,k=kl,comparison='sameK_full_Hellinger_KM_vs_Leiden',ARI=adjusted_rand_score(km.labels_,ll)))
  if name=='R12' and diagnostics[-1]['ARI']<.8:
   bridge=fitkm(zr,kl);joblib.dump(bridge,out/'bridge_PCA_KM.joblib');diagnostics.extend([dict(input=name,k=kl,comparison='full_KM_vs_PCA_KM',ARI=adjusted_rand_score(km.labels_,bridge.labels_)),dict(input=name,k=kl,comparison='PCA_KM_vs_Leiden',ARI=adjusted_rand_score(bridge.labels_,ll))]);js(out/'bridge_trigger.json',dict(trigger='sameK reference ARI <0.8, predefined blind divergence flag',MP10_used=False))
 else:diagnostics.append(dict(input=name,k=kl,comparison='not_run_K_outside_2_20',ARI=None))
 tab(out/'sameK_diagnostics.tsv',diagnostics)
 # One shared UMAP of the reference only, explicitly display subset.
 import umap
 log(name+' UMAP display')
 um=umap.UMAP(n_neighbors=30,min_dist=.3,random_state=123,n_jobs=1);u=um.fit_transform(zr);np.save(out/'umap_reference.npy',u);joblib.dump(um,out/'umap_reference_model.joblib')
 fig,axes=plt.subplots(1,2,figsize=(13,5))
 for ax,alg in zip(axes,['KM','LEIDEN']):
  labs=np.load(out/alg/'reference_native_labels.npy');sc=ax.scatter(u[:,0],u[:,1],c=labs,s=2,cmap='tab20');ax.set_title(name+' '+alg+' | shared reference UMAP');ax.set_xlabel('UMAP1');ax.set_ylabel('UMAP2');fig.colorbar(sc,ax=ax)
 fig.tight_layout();fig.savefig(out/'UMAP_comparison.png',dpi=180);plt.close(fig)
 tab(out/'summary.tsv',summaries);js(out/'COMPLETE.json',dict(status='complete',MP10_values_loaded=False))

def main():
 features=np.load(OUT/'03_features/features.npz');meta=pd.read_parquet(OUT/'03_features/common_anchors.parquet');ref=features['reference'];rng=np.random.default_rng(123);evalix=np.sort(rng.choice(len(ref),min(2000,len(ref)),replace=False));tab(OUT/'05_comparison/silhouette_evaluation_UIDs.tsv',meta.iloc[ref[evalix]][['spot_uid']]);subs=[]
 refmeta=meta.iloc[ref].reset_index(drop=True)
 for seed in SEEDS:
  sub,_=sample_hierarchy(refmeta,int(.8*len(ref)),seed);subs.append(sub);tab(OUT/f'05_comparison/stability_UIDs_{seed}.tsv',refmeta.iloc[sub][['spot_uid']])
 js(OUT/'00_preflight/software.json',{m:importlib.metadata.version(m) for m in ['numpy','pandas','scipy','scikit-learn','igraph','leidenalg','umap-learn','pyarrow']})
 for name,key in [('R1','X1'),('R2','X2'),('R12','X12')]:run_input(name,key,meta,features,ref,evalix,subs)
 summaries=pd.concat([pd.read_csv(OUT/'04_clustering'/n/'summary.tsv',sep='\t') for n in ['R1','R2','R12']],ignore_index=True);tab(OUT/'05_comparison/six_scheme_comparison.tsv',summaries)
 candidates=pd.concat([pd.read_csv(OUT/'04_clustering'/n/'candidates.tsv',sep='\t') for n in ['R1','R2','R12']],ignore_index=True);tab(OUT/'05_comparison/all_candidates.tsv',candidates)
 labs={};alllabs={}
 for c in summaries.config:
  n,a=c.rsplit('_',1);labs[c]=np.load(OUT/'04_clustering'/n/a/'reference_native_labels.npy');alllabs[c]=pd.read_parquet(OUT/'04_clustering'/n/a/'labels.parquet')['cluster_label_native_application'].to_numpy()
 aris=[];by=[]
 for a,b in itertools.combinations(labs,2):
  aris.append(dict(config1=a,config2=b,ARI=adjusted_rand_score(labs[a],labs[b]),k1=len(set(labs[a])),k2=len(set(labs[b]))));tab(OUT/'05_comparison'/f'overlap__{a}__{b}.tsv',pd.crosstab(labs[a],labs[b]).reset_index())
  for field in ['dataset_id','patient_uid','section_id']:
   for unit,g in refmeta.groupby(field):
    ix=g.index.to_numpy();by.append(dict(config1=a,config2=b,level=field,unit=unit,n=len(ix),ARI=adjusted_rand_score(labs[a][ix],labs[b][ix])))
 tab(OUT/'05_comparison/pairwise_ARI.tsv',aris);tab(OUT/'05_comparison/ARI_by_dataset_patient_section.tsv',by)
 # Full application maps retain each dataset and section, equal aspect and fixed cluster colors per config.
 for (ds,s),m in meta.groupby(['dataset_id','section_id']):
  od=OUT/'05_comparison'/ds;od.mkdir(exist_ok=True);fig,axes=plt.subplots(2,3,figsize=(15,10));ix=m.index.to_numpy()
  for ax,c in zip(axes.ravel(),['R1_KM','R2_KM','R12_KM','R1_LEIDEN','R2_LEIDEN','R12_LEIDEN']):
   ax.scatter(m.array_col,m.array_row*np.sqrt(3),c=alllabs[c][ix],s=4,cmap='tab20',vmin=0,vmax=max(alllabs[c]));ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(c);ax.set_xlabel('array_col');ax.set_ylabel('sqrt(3) array_row')
  fig.suptitle(f'{ds} / {s} | anchors={len(m)} | region cnv_candidate_v1');fig.tight_layout();fig.savefig(od/(s+'__six_schemes.png'),dpi=140);plt.close(fig)
 files=[OUT/'03_features/features.npz',OUT/'03_features/common_anchors.parquet',OUT/'03_features/reference_anchor_ids.tsv',OUT/'00_preflight/input_contract.json']+list((OUT/'code').glob('*.py'))
 manifest=dict(created=time.strftime('%F %T'),blind_labels_frozen=True,region_final_frozen=False,MP10_unblinding_allowed=False,region_reason='pending human ROI review',hashes={str(p):sha(p) for p in files},selection_rule='within 0.02 of best silhouette: smallest K then parameter',k_grid=KS,resolution_grid=RES,seeds=SEEDS,SNN='exact Euclidean 30NN incl self, Jaccard shared-neighbor sparse product; prune <=1/15; undirected no self loops',Leiden='leidenalg RBConfigurationVertexPartition, n_iterations=10',transfer='15 equal votes; exclude self in reference QC; tie mean distance then stable integer ID; distance outlier > reference leave-self-out nearest distance Q99',UMAP='reference-only shared display; not selection input',source_validation='required; original top20 AUCell; counts-first pseudospots; same labels; composition-only control; high Epi gate in input_contract')
 js(OUT/'05_comparison/frozen_manifest.json',manifest);log('Six schemes and blind comparison complete')
if __name__=='__main__':main()
