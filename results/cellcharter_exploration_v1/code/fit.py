from common import *
import numpy as np,pandas as pd,anndata as ad,cellcharter as cc,torch,importlib.metadata as md,pickle,gc
from scipy import sparse
from scipy.spatial import Delaunay
from scipy.sparse.csgraph import connected_components
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'savefig.transparent':True,'figure.facecolor':'none'})
torch.set_num_threads(8);torch.cuda.set_per_process_memory_fraction(.4,0)
o,z,xy=load();out=pd.read_parquet(ROOT/'cellcharter_assignments.parquet');thresholds={};graphrows=[];hoprows=[]
params=dict(random_state=123,batch_size=16384,covariance_regularization=1e-5,trainer_params={'accelerator':'gpu','devices':1,'max_epochs':100,'enable_progress_bar':False,'logger':False})
js(ROOT/'metadata/parameters.json',dict(versions={k:md.version(k) for k in ['cellcharter','squidpy','anndata','torchgmm','torch']},seed=123,model=params,graph='Per-section scipy Delaunay (same triangulation underlying generic Squidpy); PRIMARY cutoff min(99th percentile,40um); identical cutoff for sensitivities',aggregate='CellCharter aggregate_neighbors n_layers=3 mean sample_key=section_id',gpu_fraction=.4))
def build(name):
 ix=np.flatnonzero(out[name+'_included']);obs=out.iloc[ix][['cell_id','section_id']].copy();obs['section_id']=obs.section_id.astype('category');a=ad.AnnData(sparse.csr_matrix((len(ix),1)),obs=obs);a.obsm['X_scVI']=z[ix];a.obsm['spatial']=xy[ix];rr=[];ccs=[]
 for sec in ['Cervilla','D1_Johnson']:
  si=np.flatnonzero(obs.section_id.eq(sec));p=xy[ix[si]];duplicates=len(p)-len(np.unique(p,axis=0));assert duplicates==0,('duplicate_coordinates',sec,duplicates)
  tri=Delaunay(p).simplices;edges=np.unique(np.sort(np.vstack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[0,2]]]),axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1)
  if name=='PRIMARY':thresholds[sec]=min(float(np.quantile(length,.99)),40.)
  keep=length<=thresholds[sec];ed=edges[keep];rr.extend([si[ed[:,0]],si[ed[:,1]]]);ccs.extend([si[ed[:,1]],si[ed[:,0]]]);g=sparse.csr_matrix((np.ones(2*len(ed)),(np.r_[ed[:,0],ed[:,1]],np.r_[ed[:,1],ed[:,0]])),shape=(len(si),len(si)));degree=np.diff(g.indptr);nc,lab=connected_components(g);sizes=np.bincount(lab)
  out.loc[out.index[ix[si]],name+'_context']=np.where(degree==0,'isolated_self_expression_only','spatial_context');graphrows.append(dict(version=name,section_id=sec,n=len(si),threshold_um=thresholds[sec],edges=len(ed),removed_edges=int((~keep).sum()),edge_median_um=float(np.median(length[keep])),edge_p95_um=float(np.quantile(length[keep],.95)),neighbors_median=float(np.median(degree)),isolated=int((degree==0).sum()),components=nc,largest_component=int(sizes.max()),largest_components=json.dumps(sorted(sizes.tolist(),reverse=True)[:10])))
  rng=np.random.default_rng(123)
  for center in rng.choice(len(si),min(300,len(si)),replace=False):
   seen={int(center)};front={int(center)}
   for hop in [1,2,3]:
    nxt=set()
    for n in front:nxt.update(g.indices[g.indptr[n]:g.indptr[n+1]])
    nxt-=seen;seen|=nxt;front=nxt;reach=np.linalg.norm(p[list(nxt)]-p[center],axis=1) if nxt else np.array([]);hoprows.append(dict(version=name,section_id=sec,cell_id=obs.iloc[si[center]].cell_id,hop=hop,neighbor_count=len(nxt),max_distance_um=float(reach.max()) if len(reach) else 0,median_distance_um=float(np.median(reach)) if len(reach) else 0))
  if name=='PRIMARY':
   fig,ax=plt.subplots(figsize=(7,7));ax.scatter(*p.T,s=.1,c='gray');long=ed[np.argsort(length[keep])[-100:]]
   for e in long:ax.plot(p[e,0],p[e,1],c='red',lw=.3)
   ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(sec+' longest100 retained edges; cap %.2fum'%thresholds[sec]);ax.set_xlabel('um');ax.set_ylabel('um');fig.savefig(ROOT/f'figures/{sec}_graph_long_edges.png',dpi=160,bbox_inches='tight');plt.close(fig)
 a.obsp['spatial_connectivities']=sparse.csr_matrix((np.ones(sum(len(v) for v in rr),dtype='float32'),(np.concatenate(rr),np.concatenate(ccs))),shape=(len(ix),len(ix)));ar,ac=a.obsp['spatial_connectivities'].nonzero();assert (obs.section_id.to_numpy()[ar]==obs.section_id.to_numpy()[ac]).all()
 cc.gr.aggregate_neighbors(a,n_layers=3,aggregations='mean',use_rep='X_scVI',sample_key='section_id',connectivity_key='spatial_connectivities');assert a.obsm['X_cellcharter'].shape==(len(ix),80);assert np.isfinite(a.obsm['X_cellcharter']).all();pd.DataFrame(graphrows).to_csv(ROOT/'metadata/graph_summary.tsv',sep='\t',index=False);pd.DataFrame(hoprows).to_parquet(ROOT/'metadata/hop_sample.parquet');js(ROOT/'metadata/physical_thresholds.json',thresholds);print('AGGREGATED',name,a.shape,a.obsm['X_cellcharter'].shape,flush=True);return a,ix

if (ROOT/'models/primary_gmm.pkl').exists():
 K=json.loads((ROOT/'metadata/K_selection.json').read_text())['K'];ix=np.flatnonzero(out.PRIMARY_included);labels=out.iloc[ix].PRIMARY_niche.str.replace('niche','',regex=False).astype(int).to_numpy();thresholds=json.loads((ROOT/'metadata/physical_thresholds.json').read_text());graphrows=pd.read_csv(ROOT/'metadata/graph_summary.tsv',sep='\t').query('version == "PRIMARY"').to_dict('records');hoprows=pd.read_parquet(ROOT/'metadata/hop_sample.parquet').query('version == "PRIMARY"').to_dict('records');print('RESUME_PRIMARY_COMPLETE',K,flush=True)
else:
 a,ix=build('PRIMARY');np.save(ROOT/'metadata/primary_aggregate.npy',a.obsm['X_cellcharter']);auto=cc.tl.ClusterAutoK(n_clusters=(2,12),max_runs=5,model_params=params.copy());auto.fit(a,use_rep='X_cellcharter')
 def curve(auto,stage):
  rows=[]
  for i,k in enumerate(auto.n_clusters):
   labels=np.asarray(auto.predict(a,k=k)).astype(int);counts=np.bincount(labels,minlength=k);rows.append(dict(stage=stage,K=k,stability=float(np.mean(auto.stability[i])),sd=float(np.std(auto.stability[i])),runs=len(auto.labels[k]),min_cluster_n=int(counts.min()),min_cluster_fraction=float(counts.min()/len(labels)),cluster_sizes=json.dumps(counts.tolist())))
  return rows
 rows=curve(auto,'initial');peaks=list(auto.peaks)
 if auto.best_k==12 or not peaks:
  # One modest extension, never a parameter sweep; reuse existing fitted K values through a separate bounded candidate range.
  ext=cc.tl.ClusterAutoK(n_clusters=(12,16),max_runs=5,model_params=params.copy());ext.fit(a,use_rep='X_cellcharter');rows+=curve(ext,'extension');candidates=[(auto,k) for k in auto.peaks]+[(ext,k) for k in ext.peaks]
 else:candidates=[(auto,k) for k in peaks]
 if not candidates:candidates=[(auto,auto.best_k)]+([(ext,ext.best_k)] if 'ext' in locals() else [])
 def score(pair):
  obj,k=pair;v=float(np.mean(obj.stability[obj.n_clusters.index(k)]));lab=np.asarray(obj.predict(a,k=k)).astype(int);small=np.bincount(lab,minlength=k).min()/len(lab);return (small>=.005,v)
 chosen,K=max(candidates,key=score);K=int(K);labels=np.asarray(chosen.predict(a,k=K)).astype(int);out.loc[out.index[ix],'PRIMARY_niche']=['niche'+str(k) for k in labels];out.loc[out.index[ix],'PRIMARY_raw_cluster']=labels
 with open(ROOT/'models/primary_gmm.pkl','wb') as f:pickle.dump(chosen.best_models[K],f)
 pd.DataFrame(rows).to_csv(ROOT/'metadata/autok_stability.tsv',sep='\t',index=False);js(ROOT/'metadata/K_selection.json',dict(K=K,peaks_initial=[int(k) for k in peaks],extended='ext' in locals(),reason='Highest official stability among internal peaks without clusters below0.5%; if none eligible retain best available and flag ambiguity. No biological target used.',run_seeds=[int(s.generate_state(1)[0]) for s in np.random.SeedSequence(123).spawn(5)]));out.to_parquet(ROOT/'cellcharter_assignments.parquet',compression='zstd');print('CHOSEN_K',K,flush=True);del a,auto;gc.collect()
with open(ROOT/'models/primary_gmm.pkl','rb') as f:
 selected_model=pickle.load(f)
params['random_state']=int(selected_model.random_state)
metadata=json.loads((ROOT/'metadata/parameters.json').read_text())
metadata['autok_base_seed']=123
metadata['selected_primary_and_sensitivity_seed']=params['random_state']
metadata['model']=params
js(ROOT/'metadata/parameters.json',metadata)
print('FIXED_SELECTED_MODEL_SEED',params['random_state'],flush=True)
del selected_model
primary=pd.Series(labels,index=out.index[ix]);sens=[]
for name in ['EXPANDED_LOWINFO','NO_CURRENT_CONFLICT']:
 a,ii=build(name);model=cc.tl.GaussianMixture(n_clusters=K,**params);model.fit(a.obsm['X_cellcharter']);raw=model.predict(a.obsm['X_cellcharter']).astype(int);new=pd.Series(raw,index=out.index[ii]);shared=primary.index.intersection(new.index);ct=pd.crosstab(primary.loc[shared],new.loc[shared]).reindex(index=range(K),columns=range(K),fill_value=0);r,c=linear_sum_assignment(-ct.to_numpy());mapping={int(k):int(v) for k,v in zip(c,r)};out.loc[out.index[ii],name+'_niche']=['niche'+str(mapping[int(k)]) for k in raw];out.loc[out.index[ii],name+'_raw_cluster']=raw;js(ROOT/f'metadata/{name}_display_mapping.json',mapping)
 for sec in ['joint','Cervilla','D1_Johnson']:
  ids=shared if sec=='joint' else shared[out.loc[shared,'section_id'].eq(sec)];cross=pd.crosstab(primary.loc[ids],new.loc[ids]);cross.to_csv(ROOT/f'metadata/{name}_{sec}_raw_crosstab.tsv',sep='\t');sens.append(dict(version=name,section_id=sec,common_cells=len(ids),ARI=adjusted_rand_score(primary.loc[ids],new.loc[ids]),mapped_agreement=float((primary.loc[ids].to_numpy()==new.loc[ids].map(mapping).to_numpy()).mean())))
 out.to_parquet(ROOT/'cellcharter_assignments.parquet',compression='zstd');pd.DataFrame(sens).to_csv(ROOT/'metadata/sensitivity.tsv',sep='\t',index=False);del a,model;gc.collect();print('SENSITIVITY_DONE',name,flush=True)
print('FITTING_COMPLETE',flush=True)
