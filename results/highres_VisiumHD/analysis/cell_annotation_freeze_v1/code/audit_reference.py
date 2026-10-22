# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Read-only audit of the explicitly frozen v2 scRNA taxonomy and HD gene coverage."""
from pathlib import Path
import anndata as ad,pandas as pd,numpy as np,h5py,json,hashlib
out=Path(__file__).resolve().parents[1];hd=out.parents[1];project=hd.parent
src=project/'subtype/v2/reference_training_v2/11_final_spatial_reference/reference_v2_final_30type_regression_input.h5ad'
a=ad.read_h5ad(src,backed='r');o=a.obs.copy();v=a.var.copy()
assert o.index.is_unique and v.index.is_unique and len(o)==300033
label='reference_celltype_v2'
o['major_label']=o.major_celltype.astype(str)
o.loc[o[label].eq('Epithelial'),'major_label']='Epithelial' # explicit frozen Cycling->Epithelial reference policy
counts=o.groupby([label,'major_label'],observed=True).size().rename('n_cells').reset_index()
counts['n_source_specimens']=counts[label].map(o.groupby(label,observed=True)['orig.ident'].nunique())
counts['n_datasets']=counts[label].map(o.groupby(label,observed=True)['dataset'].nunique())
counts.to_csv(out/'metadata/reference_label_dictionary.tsv',sep='\t',index=False)
o.groupby(['orig.ident','dataset',label],observed=True).size().rename('n_cells').reset_index().to_csv(out/'metadata/reference_specimen_subtype_counts.tsv',sep='\t',index=False)
mapping=pd.read_csv(project/'subtype/v2/final_annotation_pre_reference/frozen_annotation_barcode_mapping_v2.csv',index_col=0)
inc=mapping[mapping.spatial_reference_include.eq(True)]
assert set(o.index)==set(inc.index)
assert (o[label].astype(str)==inc.loc[o.index,label].astype(str)).all()
# Match by stable Ensembl identity first; unambiguous exact current symbol only as fallback.
ens={}
for name,row in v.iterrows():
 for e in str(row.ensembl_gene_id).split('|'):
  if e.startswith('ENSG'):ens.setdefault(e.split('.')[0],set()).add(name)
res=[]
paths={
'Cervilla':hd/'data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um/filtered_feature_bc_matrix.h5',
'D1_Johnson':hd/'data/GSE300692/D1_Johnson/outs/binned_outputs/square_002um/filtered_feature_bc_matrix.h5',
'UC_SCC':hd/'data/GSE299573/GSM9040999/square_008um/filtered_feature_bc_matrix.h5'}
for dataset,p in paths.items():
 with h5py.File(p) as f:
  ids=f['matrix/features/id'].asstr()[:];names=f['matrix/features/name'].asstr()[:]
 rows=[]
 for i,(e,n) in enumerate(zip(ids,names)):
  cand=ens.get(e.split('.')[0],set());method='Ensembl_ID'
  if not cand and n in v.index:cand={n};method='exact_symbol_fallback'
  if len(cand)!=1:target='';method='ambiguous_or_not_in_reference'
  else:target=next(iter(cand))
  rows.append(dict(dataset_id=dataset,hd_feature_index=i,feature_id=e,source_symbol=n,reference_symbol=target,mapping_method=method))
 df=pd.DataFrame(rows);df.to_csv(out/f'metadata/gene_mapping_{dataset}.tsv',sep='\t',index=False)
 res.append(dict(dataset_id=dataset,hd_features=len(ids),mapped_hd_features=int(df.reference_symbol.ne('').sum()),unique_reference_genes=df.loc[df.reference_symbol.ne(''),'reference_symbol'].nunique(),duplicate_hd_symbols=int(pd.Index(names).duplicated().sum()),matrix=str(p)))
pd.DataFrame(res).to_csv(out/'metadata/reference_platform_overlap.tsv',sep='\t',index=False)
excluded=mapping[~mapping.spatial_reference_include.eq(True)].groupby(['major_celltype','subtype_v2','spatial_reference_exclude_reason'],dropna=False).size().rename('n_cells').reset_index();excluded.to_csv(out/'metadata/reference_excluded_populations.tsv',sep='\t',index=False)
manifest=dict(object_path=str(src),object_use='read_only; raw integer count X; labels selected from preserved 32-type scRNA field',reference_version='cell2location_v2_frozen_20260916__scRNA32',label_field=label,major_field='major_celltype with frozen Cycling-to-Epithelial policy when reference_celltype_v2=Epithelial',n_cells=a.n_obs,n_genes=a.n_vars,n_subtypes=32,n_source_specimens=o['orig.ident'].nunique(),datasets=o.dataset.value_counts().to_dict(),donor_id='not present; orig.ident is specimen, not a verified independent donor',membership_exact_match_to_frozen_mapping=True,label_exact_match_to_frozen_mapping=True,gene_harmonization=a.uns.get('minimal_harmonization',{}),spatial_30_merge_status='PENDING_REVIEW in final 30-type report; do not silently replace preserved scRNA labels',known_missing_or_unsupported=['IFN-stimulated B excluded by frozen decision','Independent Cycling B not an allowed identity; cycling GC B remains GC B','FOLR2 TAM vs C1QC TAM and Naive B vs Memory B have known spatial identifiability risks','vCAF provisionally included; must assess against Pericyte','Epithelial only broad identity; no validated malignancy field used'],counts_provenance='Derived minimally harmonized raw counts then gene-filtered; source scripts/count-conservation audits retained in subtype/v2; no learned signature or latent used as counts')
(out/'metadata/reference_manifest.json').write_text(json.dumps(manifest,indent=2,default=lambda x:x.tolist() if hasattr(x,'tolist') else str(x)))
print(json.dumps({k:manifest[k] for k in ['n_cells','n_genes','n_subtypes','n_source_specimens']},indent=2));print(pd.DataFrame(res).drop(columns='matrix').to_string(index=False))
