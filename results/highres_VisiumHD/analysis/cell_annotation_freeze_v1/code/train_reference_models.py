# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Custom sparse CellTypist models using audited v2 labels and conservative source-group holdout."""
import os
from pathlib import Path
root=Path(__file__).resolve().parents[1];project=root.parents[2]
os.environ.update(CELLTYPIST_FOLDER=str(root/'models/celltypist_cache'),OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:8])
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,celltypist,json,gc
np.random.seed(123);rng=np.random.default_rng(123)
manifest=json.loads((root/'metadata/reference_manifest.json').read_text());a=ad.read_h5ad(manifest['object_path'],backed='r');obs=a.obs.copy();obs['major_label']=obs.major_celltype.astype(str);obs.loc[obs.reference_celltype_v2.eq('Epithelial'),'major_label']='Epithelial'
# This old project contributes ONLY documented source-group metadata, never labels or scores.
group_source=project/'mp10_k32_reference_identifiability_v5/scope_frozen_resume_v1/01_scope/patient_mapping_status.tsv'
g=pd.read_csv(group_source,sep='\t');gm=g.set_index('sample')['group'];obs['source_group']=(obs.dataset.astype(str)+'|'+obs['orig.ident'].astype(str)).map(gm);assert obs.source_group.notna().all()
g.to_csv(root/'metadata/reference_source_group_provenance.tsv',sep='\t',index=False)
holdgroups=[]
for ds,frame in obs.groupby('dataset',observed=True):
 groups=np.array(sorted(frame.source_group.unique()));holdgroups.extend(rng.choice(groups,max(1,round(.2*len(groups))),replace=False).tolist())
obs['split']=np.where(obs.source_group.isin(holdgroups),'heldout_source_group','train')
selected=[]
for (split,group,label),inds in obs.groupby(['split','source_group','reference_celltype_v2'],observed=True).indices.items():
 cap=40 if split=='heldout_source_group' else 100
 selected.extend(rng.choice(inds,min(cap,len(inds)),replace=False))
selected=np.array(sorted(selected));o=obs.iloc[selected].copy();assert set(o[o.split=='train'].source_group).isdisjoint(set(o[o.split!='train'].source_group))
assert set(o[o.split=='train'].reference_celltype_v2)==set(obs.reference_celltype_v2)
gm=pd.read_csv(root/'metadata/gene_mapping_Cervilla.tsv',sep='\t').dropna(subset=['reference_symbol']);genes=pd.Index(sorted(gm.reference_symbol.unique()));vi=a.var_names.get_indexer(genes);assert (vi>=0).all()
print('READ stratified reference',len(selected),'genes',len(genes),flush=True)
b=a[selected,vi].to_memory();b.obs=o;b.X=b.X.tocsr().astype('float32');b.X.indices=b.X.indices.astype('int32');b.X.indptr=b.X.indptr.astype('int32');assert b.X.data.min()>=0 and np.equal(b.X.data,np.floor(b.X.data)).all()
train=b[b.obs.split.eq('train')].copy();held=b[~b.obs.split.eq('train')].copy();held.write_h5ad(root/'inputs/reference_holdout_temporary.h5ad',compression='gzip');del b;gc.collect()
train.obs[['reference_celltype_v2','major_label','dataset','orig.ident','source_group','split']].to_parquet(root/'metadata/reference_training_selection.parquet',compression='zstd');held.obs.to_parquet(root/'metadata/reference_holdout_selection.parquet',compression='zstd')
sc.pp.normalize_total(train,target_sum=1e4);sc.pp.log1p(train)
(root/'models').mkdir(exist_ok=True)
models={}
for name,label,mask in [('major','major_label',np.ones(train.n_obs,bool)),('subtype','reference_celltype_v2',np.ones(train.n_obs,bool)),('b_plasma','reference_celltype_v2',train.obs.major_label.isin(['B','Plasma']).to_numpy())]:
 dest=root/f'models/celltypist_{name}.pkl'
 if dest.exists():print('REUSE current-run model',name,flush=True);continue
 view=train[mask].copy();print('TRAIN',name,view.shape,flush=True)
 model=celltypist.train(view,labels=label,with_mean=False,use_SGD=True,mini_batch=False,max_iter=500,n_jobs=8,random_state=123,class_weight='balanced',alpha=.0001,version='HD_freeze_v1_scRNA32_seed123',details='Audited frozen v2 labels; source-group holdout; sparse log1p(1e4) common measured genes; no spatial information',source=manifest['object_path']);model.write(str(dest));del model,view;gc.collect()
config=dict(reference_version=manifest['reference_version'],source_object=manifest['object_path'],label_field='reference_celltype_v2',n_train=train.n_obs,n_heldout=held.n_obs,n_common_genes=len(genes),heldout_groups=holdgroups,n_train_source_groups=train.obs.source_group.nunique(),n_heldout_source_groups=held.obs.source_group.nunique(),verified_independent_donor_holdout=False,validation_limit='Conservative source bundles from recorded naming rules; not verified independent donors, and not cross-platform accuracy',sampling='max100 per train source_group x subtype; max40 per heldout stratum; rare strata kept completely',normalization='Restrict both reference and HD modeling views to same mapped measured genes, raw-count sum for duplicate identities, normalize_total 1e4 then log1p; master HD retains every input feature',training='CellTypist1.7.1 SGD, with_mean=False, alpha1e-4, balanced classes, max_iter500, random_state123',gene_list=genes.tolist(),models=['major','subtype','b_plasma'])
(root/'metadata/classifier_config.json').write_text(json.dumps(config,indent=2));print('TRAINING DONE',train.n_obs,held.n_obs,flush=True)
