"""Two-section raw-count scVI, section-balanced HVG selection, expression graphs only."""
from common import *
setup_gpu_env()
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,scvi,torch,time,gc
from scipy import sparse
scvi.settings.seed=123;scvi.settings.dl_num_workers=0;np.random.seed(123);torch.set_num_threads(8)
assert torch.cuda.is_available(),'CUDA must be available before training';torch.cuda.set_per_process_memory_fraction(.35,0)
names=['Cervilla','D1_Johnson'];data=[ad.read_h5ad(ROOT/n/'cells.h5ad') for n in names];shared=data[0].var_names.intersection(data[1].var_names,sort=False);assert len(shared)>15000
ranks=[];core_sizes={};hvgrows=[]
for name,a in zip(names,data):
 core=a.obs.model_eligible.to_numpy();core_sizes[name]=int(core.sum());v=a[core,shared].copy();sc.pp.normalize_total(v,target_sum=1e4);sc.pp.log1p(v);sc.pp.highly_variable_genes(v,n_top_genes=2000,flavor='seurat');score=np.nan_to_num(v.var.dispersions_norm.to_numpy(),nan=-np.inf);rank=pd.Series(-score).rank(method='average').to_numpy();ranks.append(rank);hvgrows.append(v.var.highly_variable.to_numpy());del v;gc.collect()
# Equal per-section rank contribution; first prioritize HVGs found in both sections.
nvotes=np.sum(hvgrows,axis=0);meanrank=np.mean(ranks,axis=0);order=np.lexsort((meanrank,-nvotes));hv=shared[order[:2000]]
pd.DataFrame({'feature_id':shared,'section_votes':nvotes,'mean_rank':meanrank,'selected':shared.isin(hv)}).to_csv(ROOT/'metadata/scvi_hvg.tsv',sep='\t',index=False)
for name,a in zip(names,data):
 hvg_counts=np.asarray(a[:,hv].X.sum(1)).ravel();a.obs['scvi_hvg_counts']=hvg_counts;a.obs['scvi_latent_mapped']=hvg_counts>0;a.obs['model_eligible']=a.obs.model_eligible&(hvg_counts>0);a.obs['model_role']=a.obs.model_role.astype(object);a.obs.loc[hvg_counts==0,'model_role']='prediction_only_no_hvg_counts';core_sizes[name]=int(a.obs.model_eligible.sum())
views=[a[a.obs.model_eligible,hv].copy() for a in data];train=ad.concat(views,join='inner',merge='same');assert train.obs_names.is_unique and train.obs.section_id.nunique()==2
print('SCVI_INPUT',train.shape,train.obs.section_id.value_counts().to_dict(),flush=True)
base=train.copy();sc.pp.normalize_total(base,target_sum=1e4);sc.pp.log1p(base);sc.tl.pca(base,n_comps=30,random_state=123);baseline=base.obsm['X_pca'].astype('float32');del base
scvi.model.SCVI.setup_anndata(train,batch_key='section_id');model=scvi.model.SCVI(train,n_latent=20,n_hidden=128,n_layers=2,gene_likelihood='nb',dispersion='gene');start=time.monotonic()
if (ROOT/'models/scvi_two_sections/model.pt').exists():
 model=scvi.model.SCVI.load(ROOT/'models/scvi_two_sections',adata=train,accelerator='gpu',device=0);prior=json.loads((ROOT/'metadata/scvi_training_summary.json').read_text());training_seconds=prior['training_seconds'];trained_epochs=prior['trained_epochs'];print('REUSE_FINISHED_SCVI_MODEL',flush=True)
else:
 model.train(max_epochs=50,accelerator='gpu',devices=1,batch_size=1024,train_size=.9,early_stopping=True,early_stopping_patience=10,check_val_every_n_epoch=1,enable_progress_bar=False,enable_checkpointing=False,logger=None,default_root_dir=str(ROOT/'logs/scvi_training'),plan_kwargs={'lr':1e-3,'n_epochs_kl_warmup':10})
 training_seconds=time.monotonic()-start;trained_epochs=len(model.history.get('elbo_train',[]));print('SCVI_TRAINING_DONE',trained_epochs,training_seconds,flush=True)
 model.save(ROOT/'models/scvi_two_sections',overwrite=False,save_anndata=False)
 for key,h in model.history.items():h.to_csv(ROOT/f'logs/scvi_{key}.csv')
 jwrite(ROOT/'metadata/scvi_training_summary.json',dict(training_seconds=training_seconds,trained_epochs=trained_epochs,training_cells=train.n_obs,sections=core_sizes))
latents=[]
for a in data:
 mapped=a.obs.scvi_latent_mapped.to_numpy();z=np.full((a.n_obs,20),np.nan,dtype='float32');z[mapped]=model.get_latent_representation(adata=a[mapped,hv].copy(),batch_size=2048).astype('float32');latents.append(z)
corelat=np.vstack([z[a.obs.model_eligible] for a,z in zip(data,latents)]);s=ad.AnnData(np.zeros((train.n_obs,1),dtype='float32'),obs=train.obs.copy());s.obsm['X_scVI']=corelat
print('EXPRESSION_GRAPH',s.n_obs,flush=True);sc.pp.neighbors(s,use_rep='X_scVI',n_neighbors=15,random_state=123);sc.tl.leiden(s,resolution=.8,random_state=123,flavor='igraph',n_iterations=2,key_added='expression_cluster');sc.tl.umap(s,random_state=123)
s.obs['lineage_cluster']='not_clustered'
for major in sorted(s.obs.predicted_major.unique()):
 ix=np.flatnonzero(s.obs.predicted_major.eq(major).to_numpy())
 if len(ix)<80:continue
 l=s[ix].copy();sc.pp.neighbors(l,use_rep='X_scVI',n_neighbors=min(15,len(ix)-1),random_state=123);sc.tl.leiden(l,resolution=.7,random_state=123,flavor='igraph',n_iterations=2,key_added='within_lineage');s.obs.loc[s.obs_names[ix],'lineage_cluster']=[major+'::'+x for x in l.obs.within_lineage.astype(str)];print('LINEAGE_CLUSTER',major,len(ix),l.obs.within_lineage.nunique(),flush=True)
starti=0
for name,a,z in zip(names,data,latents):
 core=a.obs.model_eligible.to_numpy();n=int(core.sum());ids=a.obs_names[core];assert np.array_equal(ids,s.obs_names[starti:starti+n]);u=np.full((a.n_obs,2),np.nan,dtype='float32');u[core]=s.obsm['X_umap'][starti:starti+n];pca=np.full((a.n_obs,30),np.nan,dtype='float32');pca[core]=baseline[starti:starti+n]
 a.obsm['X_scVI']=z;a.obsm['X_umap']=u;a.obsm['X_pca_unintegrated']=pca;a.obs['expression_cluster']='not_in_structure_core';a.obs['lineage_cluster']='not_in_structure_core';a.obs.loc[ids,'expression_cluster']=s.obs.loc[ids,'expression_cluster'].astype(str).to_numpy();a.obs.loc[ids,'lineage_cluster']=s.obs.loc[ids,'lineage_cluster'].astype(str).to_numpy();a.obs['scvi_training_core']=core;a.obs['umap_mapped']=core;a.var['scvi_hvg']=a.var_names.isin(hv);a.uns['scvi_note']='2 sections; section_id batch; expression graph only; physical coordinates unchanged; noncore latent inferred only with positive measured HVG counts; zero-HVG latent NA; all noncore UMAP NA';a.write_h5ad(ROOT/name/'cells.h5ad',compression='gzip');starti+=n
config=dict(seed=123,n_sections=2,sections=names,structure_core_by_section=core_sizes,total_training_cells=train.n_obs,training_genes=len(hv),hvg_selection='section-wise seurat dispersion ranks, prioritize shared HVGs then equal-weight mean rank; master all features retained',raw_counts=True,batch_key='section_id',likelihood='nb',latent=20,max_epochs=50,kl_warmup_epochs=10,early_stopping_patience=10,batch_size=1024,device=torch.cuda.get_device_name(0),cuda=torch.version.cuda,torch=torch.__version__,scvi=scvi.__version__,training_seconds=training_seconds,trained_epochs=trained_epochs,runtime_seconds=time.monotonic()-start,noncore='latent inference with positive HVG counts, zero-HVG latent NA; all noncore UMAP NA; all retained in master/spatial plots',n_global_clusters=s.obs.expression_cluster.nunique(),lineage_clusters=s.obs.groupby('predicted_major',observed=True).lineage_cluster.nunique().to_dict(),limitation='One section per source: section/cohort/tissue effects confounded; no acceptance criterion based on complete mixing')
jwrite(ROOT/'metadata/scvi_config.json',config);print('SCVI_DONE',config,flush=True)
