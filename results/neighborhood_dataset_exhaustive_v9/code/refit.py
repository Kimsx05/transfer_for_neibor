# Dedicated v3 numerical runtime, with no import-time writes into prior projects.
import sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path=[p for p in sys.path if 'dist-packages' not in p]
sys.path[:0]=[str(ROOT/'neighborhood_ring_comparison_v2/runtime'),str(ROOT/'spatial_ecology_top20_v2/runtime'),str(ROOT/'cell2location_setup_v1/cell2location_env/lib/python3.12/site-packages')]
sys.path.append(str(ROOT/'neighborhood_v1/.venv/lib/python3.12/site-packages'))
from core import *
from scipy import sparse
from scipy.spatial.distance import pdist,squareform,cdist
from scipy.optimize import linear_sum_assignment
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import silhouette_score,adjusted_rand_score
from sklearn.decomposition import PCA
import igraph as ig,leidenalg,joblib,importlib.metadata
extract(V2/'code/common.py',['allocate','sample_hierarchy'],globals())
extract(V3/'code/cluster.py',['graph','leiden','neighbors','votes'],globals())

def run(cid,rebuild=False):
 start=time.time();folder=OUT/'combinations'/cid;folder.mkdir(exist_ok=True)
 if not rebuild and valid(folder,'NICHE_FROZEN'):log(cid+' NICHE_FROZEN reused');return
 guard(cid+':representation',300*1024**2)
 c=pd.read_csv(OUT/'combinations.tsv',sep='\t',dtype={'combo_id':str}).set_index('combo_id').loc[cid];members=c.datasets.split('|');allm=pd.read_parquet(OUT/'shared/master_anchors.parquet');ix=np.flatnonzero(allm.dataset_id.isin(members));m=allm.iloc[ix].reset_index(drop=True);m['sample_uid']=m.dataset_id+'::'+m.section_id
 # Metadata whitelist from common anchors contains no MP10 columns.
 assert not any('MP10' in col for col in m.columns)
 unit=pd.read_parquet(V3/'01_features/metadata.parquet',columns=['spot_uid','analysis_unit_uid']);m=m.merge(unit,on='spot_uid',validate='one_to_one')
 if 'patient_stratum' not in m:
  m['patient_stratum']=np.where(m.source_patient_identity_known,m.patient_uid,'SECTION::'+m.section_id)
 ref,alloc=sample_hierarchy(m,20000,123);assert len(np.unique(ref))==len(ref) and len(ref)>=31
 f=np.load(V2/'03_features/features.npz');x=f['X12'][ix].astype('float64');q=np.concatenate([f['q1'][ix],f['q2'][ix]],1).astype(float)
 pca=PCA(n_components=.95,svd_solver='full',random_state=123).fit(x[ref]);z=pca.transform(x)
 inp=temp_path(cid,'pca.tsv');meta=temp_path(cid,'metadata.tsv');emb=temp_path(cid,'embedding.tsv');logpath=temp_path(cid,'harmony.log');obj=temp_path(cid,'objective.tsv');hm=temp_path(cid,'harmony_manifest.json');rs=temp_path(cid,'R_session.txt')
 np.savetxt(inp,z,delimiter='\t',fmt='%.17g');m[['spot_uid','sample_uid']].to_csv(meta,sep='\t',index=False)
 exchange=folder/'temporary_representation.json'
 reuse=(not rebuild) and valid(folder,'REPRESENTATION_DONE') and exchange.exists() and emb.exists() and hm.exists()
 if reuse:
  checkpoint=json.loads(exchange.read_text());reuse=checkpoint['embedding_hash']==sha(emb) and checkpoint['harmony_hash']==sha(hm) and checkpoint['reference_hash']==arrhash(ix[ref])
 if not reuse:
  with open(logpath,'w') as logf:subprocess.run(['/usr/local/bin/Rscript',str(OUT/'code/harmony.R'),str(inp),str(meta),'sample_uid','1','5',str(inp.parent)],check=True,stdout=logf,stderr=subprocess.STDOUT,env=dict(os.environ,TMPDIR='/tmp/neighborhood_v9_rtemp_t070721_20261018'))
  if not rebuild:js(exchange,dict(embedding_hash=sha(emb),harmony_hash=sha(hm),reference_hash=arrhash(ix[ref]),note='Exchange is regenerable after NICHE_FROZEN; absence is normal after cleanup'))
 else:log(cid+' exact temporary Harmony checkpoint reused')
 z=np.loadtxt(emb,delimiter='\t');assert z.shape==(len(m),pca.n_components_) and np.isfinite(z).all()
 harmony=json.loads(hm.read_text())
 if rebuild:
  frozen=np.load(folder/'reference.npz');assert np.array_equal(frozen['reference_index'],ix[ref]);model=np.load(folder/'pca_model.npz');assert np.array_equal(model['components'],pca.components_) and np.array_equal(model['mean'],pca.mean_)
  dest=OUT/'review'/cid;dest.mkdir(parents=True,exist_ok=True);npz(dest/'reference_embedding.npz',reference_coordinates=z[ref],reference_labels=np.load(folder/'labels.npz')['reference_native'],reference_index=ix[ref]);js(dest/'representation_rebuild.json',dict(harmony=harmony,code_hash=sha(__file__),purpose='display only; frozen labels unchanged'));clean(cid);return
 js(folder/'harmony.json',harmony);table(folder/'harmony_objective.tsv',pd.read_csv(obj,sep='\t'));(folder/'R_session.txt').write_text(rs.read_text())
 npz(folder/'reference.npz',anchor_index=ix.astype(np.int32),reference_index=ix[ref].astype(np.int32),reference_local=ref.astype(np.int32));table(folder/'reference_allocation.tsv',alloc)
 npz(folder/'pca_model.npz',components=pca.components_,mean=pca.mean_,explained_variance_ratio=pca.explained_variance_ratio_)
 js(folder/'representation.json',dict(n_anchors=len(m),n_reference=len(ref),n_pc=pca.n_components_,members=members,master_uid_hash=arrhash(allm.spot_uid.to_numpy(dtype=str)),MP10_used=False,source_feature_dtype=str(f['X12'].dtype),PCA_dtype='float64',software={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','scikit-learn','igraph','leidenalg']},code_hash=sha(__file__)))
 mark(folder,'REPRESENTATION_DONE',[folder/'reference.npz',folder/'pca_model.npz',folder/'harmony.json',folder/'representation.json'])
 log(cid+' Harmony complete; PC='+str(pca.n_components_))
 guard(cid+':graph',0);zr=z[ref];g=graph(zr);ds,ids=neighbors(z,ref);threshold=float(np.quantile(ds[ref,0],.99));ev=np.sort(np.random.default_rng(123).choice(len(ref),min(2000,len(ref)),replace=False));D=squareform(pdist(zr[ev]));rawD=squareform(pdist(x[ref[ev]]));rows=[];labs={}
 def candidate(res):
  lab=leiden(g,res);labs[res]=lab;k=len(set(lab));sizes=np.bincount(lab);cent=np.array([q[ref][lab==i].mean(0) for i in range(k)]);dist=pdist(cent);pred=votes(lab,ids,ds)[0];pred[ref]=lab;tr,support=votes(lab,ids,ds,pred);accepted=(support>=.6)&(ds[:,0]<=threshold)
  sil=lambda a:float(silhouette_score(a,lab[ev],metric='precomputed')) if 1<len(set(lab[ev]))<len(ev) else np.nan
  rows.append(dict(parameter=res,k=k,silhouette_embedding=sil(D),silhouette_raw_Hellinger=sil(rawD),min_cluster=int(sizes.min()),tiny_clusters=int((sizes<max(20,.001*len(ref))).sum()),min_profile_distance=float(dist.min()) if len(dist) else np.nan,raw_composition_dispersion=float(np.mean(np.sum((q[ref]-cent[lab])**2,1))),accepted_fraction=accepted.mean(),low_support_fraction=(support<.6).mean(),distance_outlier_fraction=(ds[:,0]>threshold).mean()))
 for res in [.005,.01,.02,.05,.1]:candidate(res)
 ks=[r['k'] for r in rows];extra=[.2,.4] if max(ks)<6 else [.0005,.001,.002] if min(ks)>12 else []
 for res in extra:candidate(res)
 df=pd.DataFrame(rows);lv=df[(df.k>1)&df.tiny_clusters.eq(0)];sep=lv[(lv.silhouette_embedding>0)&(lv.silhouette_raw_Hellinger>0)];lv=sep if len(sep) else lv;moderate=lv[lv.k.between(6,12)];lv=moderate if len(moderate) else lv
 fallback=False
 if not len(lv):lv=df[df.k>1];fallback=True
 if not len(lv):lv=df;fallback=True
 best=lv.sort_values(['silhouette_embedding','accepted_fraction'],ascending=False).iloc[0];res=float(best.parameter);lab=labs[res];K=int(best.k);candidate_lab=votes(lab,ids,ds)[0];candidate_lab[ref]=lab;tr,support=votes(lab,ids,ds,candidate_lab);accepted=np.where((support>=.6)&(ds[:,0]<=threshold),candidate_lab,-1)
 st=[]
 for seed in range(124,129):
  sub,_=sample_hierarchy(m.iloc[ref].reset_index(drop=True),int(.8*len(ref)),seed);sl=leiden(graph(zr[sub]),res,seed);ct=pd.crosstab(lab[sub],sl).to_numpy();aa,bb=linear_sum_assignment(-ct);qsub=q[ref[sub]];c1=np.array([qsub[lab[sub]==c].mean(0) for c in sorted(set(lab[sub]))]);c2=np.array([qsub[sl==c].mean(0) for c in sorted(set(sl))]);st.append(dict(seed=seed,n=len(sub),k=len(set(sl)),ARI=adjusted_rand_score(lab[sub],sl),matched_fraction=ct[aa,bb].sum()/len(sub),matched_profile_L2=float(np.linalg.norm(c1[aa]-c2[bb],axis=1).mean()),kind='balanced80_fixed_representation'))
 mr=m.iloc[ref];block=mr.sample_uid+'::'+(mr.array_row//8).astype(str)+'::'+(mr.array_col//8).astype(str);uu=np.sort(block.unique());keep=np.random.default_rng(321).choice(uu,int(.8*len(uu)),replace=False);sub=np.flatnonzero(block.isin(keep));sl=leiden(graph(zr[sub]),res,321);st.append(dict(seed=321,n=len(sub),k=len(set(sl)),ARI=adjusted_rand_score(lab[sub],sl),kind='spatial8x8_fixed_representation'))
 table(folder/'stability.tsv',st);table(folder/'candidates.tsv',df)
 npz(folder/'labels.npz',anchor_index=ix.astype(np.int32),candidate=candidate_lab.astype(np.int16),accepted=accepted.astype(np.int16),transfer=tr.astype(np.int16),support=support,nearest_distance=ds[:,0],reference_native=lab.astype(np.int16))
 profiles=[];flags=[];coverage=[]
 old=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet',columns=['spot_uid','accepted_label']).set_index('spot_uid').loc[m.spot_uid,'accepted_label'].to_numpy();mapping=[];oldfull=pd.read_parquet(V3/'03_main/SAMPLE_R12_LEIDEN/labels.parquet',columns=['accepted_label']).accepted_label.to_numpy();qfull=np.concatenate([f['q1'],f['q2']],1).astype(float);oldcent=np.array([qfull[oldfull==j].mean(0) for j in range(9)])
 for k in range(K):
  use=accepted==k;mm=m[use];pr=q[use].mean(0) if use.any() else np.full(60,np.nan);dc=mm.dataset_id.value_counts(normalize=True);sc=mm.section_id.value_counts(normalize=True)
  flags.append(dict(niche=k,n_accepted=int(use.sum()),n_candidate=int((candidate_lab==k).sum()),n_datasets=mm.dataset_id.nunique(),n_sections=mm.section_id.nunique(),n_units=mm.analysis_unit_uid.nunique(),max_dataset_fraction=dc.max() if len(dc) else np.nan,max_section_fraction=sc.max() if len(sc) else np.nan,dataset_dominated=bool(len(dc) and dc.max()>.8),sample_dominated=bool(len(sc) and sc.max()>.5),small=use.sum()<.005*len(m)))
  for j,v in enumerate(pr):profiles.append(dict(niche=k,ring=['R1','R2'][j//30],celltype=TYPES[j%30],mean=v,log2_enrichment=np.log2((v+1e-6)/(q[:,j].mean()+1e-6))))
  for sec,ss in m.groupby('section_id'):
   cc=use[ss.index];coverage.append(dict(niche=k,dataset_id=ss.dataset_id.iloc[0],section_id=sec,analysis_unit_uid=ss.analysis_unit_uid.iloc[0],n_accepted=int(cc.sum()),n_common=len(ss)))
  for j in range(9):mapping.append(dict(niche=k,old_niche=j,profile_L2=float(np.linalg.norm(pr-oldcent[j])),n_overlap=int((use&(old==j)).sum()),new_n=int(use.sum()),old_subset_n=int((old==j).sum()),overlap_fraction=float((use&(old==j)).sum()/use.sum()) if use.any() else np.nan))
 table(folder/'profiles.parquet',profiles);table(folder/'niche_flags.tsv',flags);table(folder/'coverage.parquet',coverage);table(folder/'old_niche_mapping.tsv',mapping)
 plan=[dict(niche=k,cell_index=j,celltype=TYPES[j],pair_id=f'n{k}_c{j:02d}') for k in range(K) for j in NON];table(folder/'ring_plan.tsv',plan)
 choice=dict(K=K,resolution=res,extension=extra,selection=best.to_dict(),fallback_quality_limited=fallback,threshold=threshold,stability_min=min(r['ARI'] for r in st),seconds=time.time()-start,accepted_n=int((accepted>=0).sum()),label_hash=sha(folder/'labels.npz'),MP10_loaded=False,graph_components=len(g.connected_components()),graph_isolates=int((np.array(g.degree())==0).sum()),code_hash=sha(__file__))
 js(folder/'selection.json',choice)
 # Display reconstruction has returned before mutation of frozen results.
 mark(folder,'NICHE_FROZEN',[folder/x for x in ['labels.npz','selection.json','profiles.parquet','reference.npz','reference_allocation.tsv','candidates.tsv','stability.tsv','coverage.parquet','old_niche_mapping.tsv','ring_plan.tsv']]);clean(cid);log(cid+' NICHE_FROZEN K='+str(K))
if __name__=='__main__':
 import fcntl
 folder=OUT/'combinations'/sys.argv[1];folder.mkdir(exist_ok=True);lock=open(folder/'.refit.lock','w');fcntl.flock(lock,fcntl.LOCK_EX)
 try:run(sys.argv[1],rebuild='--rebuild' in sys.argv)
 except ResourcePaused as e:
  js(OUT/'combinations'/sys.argv[1]/'state.json',dict(stage='RESOURCE_PAUSED',error=str(e)));sys.exit(75)
