from common import *
import numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
rng=np.random.default_rng(123);checks=[]
for name in ['Cervilla','D1_Johnson']:
 meta,base,out=section(name);a=ad.read_h5ad(out/'cells.h5ad');o=a.obs;s=json.loads((out/'metadata/reconstruction_qc_summary.json').read_text());t=json.loads((out/'metadata/transform.json').read_text());mp=pd.read_parquet(out/'inputs/bin_to_cell.parquet');labels=pd.read_parquet(out/'cell_labels.parquet')
 def check(key,yes,detail=''):
  assert bool(yes),(name,key,detail)
  checks.append(dict(section=name,check=key,passed=True,detail=str(detail)))
 check('unique_prefixed_cell_ids',a.obs_names.is_unique and all(x.startswith(name+'__v2__cell_') for x in a.obs_names));check('labels_exact_alignment',a.obs_names.equals(labels.index));check('raw_integer_nonnegative',sparse.issparse(a.X) and np.issubdtype(a.X.dtype,np.integer) and a.X.data.min()>=0)
 check('raw_row_sums_and_gene_counts',np.array_equal(np.asarray(a.X.sum(1)).ravel(),o.n_counts) and np.array_equal(a.X.getnnz(1),o.n_genes));check('retained_total_counts',int(a.X.sum())==s['retained_counts']);check('qc_floor',((o.n_counts>=10)&(o.n_genes>=10)).all());check('low_info_preserved',int(o.low_information.sum())==s['low_information_retained']);check('retained_denominator',len(o)==s['retained']);check('all_gene_features',a.n_vars==18085)
 check('all_count_partition',s['assigned_counts']+s['unassigned_counts']==s['source_counts'] and s['retained_counts']+s['excluded_counts']==s['assigned_counts'])
 check('noncore_UMAP_explicit',np.isnan(a.obsm['X_umap'][~o.umap_mapped]).all() and np.isfinite(a.obsm['X_umap'][o.umap_mapped]).all());check('latent_all_retained',a.obsm['X_scVI'].shape==(len(o),20) and np.isfinite(a.obsm['X_scVI'][o.scvi_latent_mapped]).all() and np.isnan(a.obsm['X_scVI'][~o.scvi_latent_mapped]).all());check('spatial_physical_unchanged',np.allclose(a.obsm['spatial_um'],a.obsm['spatial']*t['native_mpp']))
 check('crop_applied_once',np.allclose(a.obsm['spatial_image'],a.obsm['spatial']-np.array(t['fullres_crop_xyxy'][:2])));check('model_transform',np.allclose(a.obsm['spatial_model'],(a.obsm['spatial_image']+.5)*t['local_to_model_scale_xy']-.5))
 with h5py.File(out/'cells.h5ad','r') as disk:
  check('master_no_duplicate_raw_layer',a.raw is None and len(disk['layers'])==0,'AnnData runtime layers[None] aliases X; on-disk layers group is empty')
 chosen=np.unique(np.r_[rng.choice(a.n_obs,min(25,a.n_obs),replace=False),int(np.argmax(o.n_counts.to_numpy()))]);source_ids=o.source_cell_id.astype(int).to_numpy();group={i:mp.loc[mp.cell_label.eq(source_ids[i]),'source_bin_index'].to_numpy() for i in chosen}
 with h5py.File(base/'filtered_feature_bc_matrix.h5') as f:
  g=f['matrix'];check('original_feature_order',np.array_equal(g['features/id'].asstr()[:],a.var_names));ptr=g['indptr'][:]
  for i,bins in group.items():
   row=np.zeros(a.n_vars,dtype=np.int64)
   for b in bins:
    lo,hi=ptr[b:b+2];np.add.at(row,g['indices'][lo:hi],g['data'][lo:hi])
   check('sample_original_integer_bin_sum_'+str(source_ids[i]),np.array_equal(row,a.X[i].toarray().ravel()),f'{len(bins)} bins; {int(row.sum())} original UMIs')
 for rel in ['inputs/nuclei.npz','inputs/cell_mask_2um_array.npz','metadata/annotation_summary.json','metadata/cluster_marker_evidence.tsv','figures/QC.png','figures/expression_UMAP.png','figures/spatial_labels.png','figures/marker_dotplot.png']:
  check('file_'+rel,(out/rel).is_file())
 print('VALIDATED',name,a.shape,flush=True)
cfg=json.loads((ROOT/'metadata/scvi_config.json').read_text());assert cfg['n_sections']==2 and set(cfg['sections'])=={'Cervilla','D1_Johnson'};pd.DataFrame(checks).to_csv(ROOT/'metadata/delivery_validation.tsv',sep='\t',index=False);jwrite(ROOT/'metadata/delivery_validation.json',dict(checks_passed=len(checks),sections=cfg['sections'],seed=123,raw_check='All row sums/gene counts and stored count conservation plus 25 random cells and highest-count cell per section independently reaggregated from original 2um H5; not a full duplicated recomputation'))
