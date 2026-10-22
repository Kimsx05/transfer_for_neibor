# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Single available HD section: scVI expression structure, not demonstrated cross-section integration."""
import os
from pathlib import Path
root=Path(__file__).resolve().parents[1]
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
if hasattr(os,'sched_setaffinity'):
 cpus=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,cpus[8:16] if len(cpus)>=16 else cpus[:8])
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,scvi,torch,json
scvi.settings.seed=123;scvi.settings.dl_num_workers=0;torch.set_num_threads(8);np.random.seed(123)
a=ad.read_h5ad(root/'cell_annotation.h5ad');core=a.obs.model_eligible.to_numpy();v=a[core].copy();sc.pp.normalize_total(v,target_sum=1e4);sc.pp.log1p(v);sc.pp.highly_variable_genes(v,n_top_genes=2000,flavor='seurat');hv=v.var.highly_variable.to_numpy();sc.tl.pca(v,n_comps=30,use_highly_variable=True,random_state=123);baseline=v.obsm['X_pca'].copy();del v
modeldata=a[core,hv].copy();print('STRUCTURE CORE',modeldata.shape,flush=True)
scvi.model.SCVI.setup_anndata(modeldata) # one section, no invented batch or donor covariate
model=scvi.model.SCVI(modeldata,n_latent=20,n_hidden=128,n_layers=2,gene_likelihood='nb',dispersion='gene')
model.train(max_epochs=50,accelerator='cpu',devices=1,batch_size=512,train_size=.9,early_stopping=True,early_stopping_patience=10,check_val_every_n_epoch=1,enable_progress_bar=False,enable_checkpointing=False,logger=False,plan_kwargs={'lr':1e-3})
model.save(root/'models/scvi_hd',overwrite=False,save_anndata=False)
allview=a[:,hv].copy();latent=model.get_latent_representation(adata=allview,batch_size=1024);assert np.isfinite(latent).all()
struct=ad.AnnData(np.zeros((core.sum(),1),dtype='float32'),obs=a.obs.loc[core].copy());struct.obsm['X_scVI']=latent[core];sc.pp.neighbors(struct,use_rep='X_scVI',n_neighbors=15,random_state=123);sc.tl.leiden(struct,resolution=.6,random_state=123,flavor='igraph',n_iterations=2);sc.tl.umap(struct,random_state=123)
umap=np.full((a.n_obs,2),np.nan,dtype='float32');umap[core]=struct.obsm['X_umap'];clusters=np.repeat('not_in_structure_core',a.n_obs).astype(object);clusters[core]=struct.obs.leiden.astype(str).to_numpy();pd.DataFrame({'cell_id':a.obs_names,'expression_cluster':clusters}).to_parquet(root/'metadata/expression_clusters.parquet',index=False,compression='zstd')
np.savez_compressed(root/'metadata/expression_embeddings.npz',cell_ids=a.obs_names.to_numpy(dtype=str),X_scVI=latent.astype('float32'),X_umap=umap,baseline_PCA_core=baseline.astype('float32'),core_mask=core)
a.var.loc[hv,['feature_id','gene_symbol']].to_csv(root/'metadata/scvi_hvg.tsv',sep='\t',index=False)
for name,hist in model.history.items():hist.to_csv(root/f'logs/scvi_{name}.csv')
config=dict(seed=123,core_n=int(core.sum()),retained_n=a.n_obs,hvg_n=int(hv.sum()),n_latent=20,n_hidden=128,n_layers=2,likelihood='nb',max_epochs=50,early_stopping_patience=10,batch_size=512,train_size=.9,accelerator='cpu',batch_key=None,n_available_sections=1,batch_interpretation='Only Cervilla reconstructable; no cross-section batch integration claim; section/donor effects cannot be separated',low_information_role='inference_only for latent; excluded from graph, Leiden and UMAP definition; classification retained',leiden_resolution=.6,n_expression_clusters=struct.obs.leiden.nunique(),counts='Original integer aggregate X, modeling view only selects HVGs',version=scvi.__version__)
(root/'metadata/scvi_config.json').write_text(json.dumps(config,indent=2));print('SCVI DONE',config,flush=True)
