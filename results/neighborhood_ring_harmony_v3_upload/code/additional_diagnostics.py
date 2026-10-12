from base import *
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
from scipy.spatial.distance import pdist,squareform
import itertools,igraph as ig
from cluster import votes
C=OUT/'04_comparison';m,f,ref=inputs();m=pd.read_parquet(OUT/'01_features/metadata.parquet');configs=sorted((OUT/'03_main').iterdir());reps=sorted({p.resolve().parent for p in configs});rows=[];profiles=[];candidate_support=[]
for od in reps:
 g=ig.Graph.Read_Pickle(str(od/'graph.pickle'));degree=np.asarray(g.degree());rows.append(dict(representation_id=od.name,n=g.vcount(),edges=g.ecount(),components=len(g.connected_components()),isolates=int((degree==0).sum()),degree_min=int(degree.min()),degree_Q25=float(np.quantile(degree,.25)),degree_median=float(np.median(degree)),degree_Q75=float(np.quantile(degree,.75)),degree_max=int(degree.max())))
 for p in (od/'candidates').iterdir():
  labs=np.load(p/'reference_labels.npy');mm=m.iloc[ref].copy();mm['cluster']=labs
  for key in ['dataset_id','sample_uid','analysis_unit_uid','patient_uid']:
   aa=mm if key!='patient_uid' else mm[mm.source_patient_identity_known]
   counts=aa.groupby(['cluster',key]).size().reset_index(name='n');counts.insert(0,'source_key',key);counts=counts.rename(columns={key:'source_uid'});counts.insert(0,'candidate',p.name);counts.insert(0,'representation_id',od.name);candidate_support.append(counts)
for p in configs:
 pr=pd.read_csv(p/'profiles_R1_R2.tsv',sep='\t');cent=pr.pivot(index='cluster',columns=['ring','cell_type'],values='mean_relative_composition').to_numpy();d=squareform(pdist(cent));cor=np.corrcoef(cent)
 for a,b in itertools.combinations(range(len(cent)),2):profiles.append(dict(config_id=p.name,cluster1=a,cluster2=b,profile_L2=d[a,b],profile_Pearson=cor[a,b]))
tab(C/'graph_degree_diagnostics.tsv',rows);tab(C/'centroid_similarity.tsv',profiles);tab(C/'candidate_source_support.tsv',pd.concat(candidate_support,ignore_index=True))
# Pairwise reference native agreement within each dataset/sample; no biological alignment assumed.
labs={p.name:np.load(p/'reference_native_labels.npy') for p in configs};mr=m.iloc[ref].reset_index(drop=True);out=[]
for a,b in itertools.combinations(labs,2):
 for key in ['dataset_id','sample_uid']:
  for unit,gg in mr.groupby(key):
   ix=gg.index.to_numpy();out.append(dict(config1=a,config2=b,source_key=key,source_uid=unit,n=len(ix),ARI=adjusted_rand_score(labs[a][ix],labs[b][ix]),NMI=normalized_mutual_info_score(labs[a][ix],labs[b][ix])))
tab(C/'pairwise_ARI_NMI_by_dataset_sample.tsv',out);log('Additional diagnostics complete')

coverage=[]
for od in reps:
 models,labs,cr=joblib.load(od/'candidates.joblib');z=np.load(od/'embedding.npy',mmap_mode='r');n=np.load(od/'neighbors.npz');ids=n['ids'];dist=n['distance'];far=dist[:,0]>np.quantile(dist[ref,0],.99)
 for tag,lab in labs.items():
  pred=models[tag].predict(z) if tag.startswith('KM') else votes(lab,ids,dist)[0];pred[ref]=lab;tr,v= votes(lab,ids,dist,pred)
  t=m[['dataset_id','sample_uid']].copy();t['low_support']=v<.6;t['distance_outlier']=far;t['accepted']=(v>=.6)&~far
  a=t.groupby(['dataset_id','sample_uid']).agg(n=('accepted','size'),accepted_fraction=('accepted','mean'),low_support_fraction=('low_support','mean'),distance_outlier_fraction=('distance_outlier','mean')).reset_index();a.insert(0,'candidate',tag);a.insert(0,'representation_id',od.name);coverage.append(a)
tab(C/'candidate_coverage_by_dataset_sample.tsv',pd.concat(coverage,ignore_index=True))
