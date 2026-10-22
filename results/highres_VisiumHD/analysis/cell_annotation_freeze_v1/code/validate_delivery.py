# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Independent checks against the original integer bins, every mapped cell and final labels."""
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
from pathlib import Path
import json,hashlib,numpy as np,pandas as pd,anndata as ad,h5py
from scipy import sparse
root=Path(__file__).resolve().parents[1];assert not (root/'metadata/FROZEN.json').exists();rows=[]
def check(name,truth,detail=''):
 rows.append(dict(check=name,passed=bool(truth),detail=str(detail)));print(name,bool(truth),detail,flush=True)
 pd.DataFrame(rows).to_csv(root/'metadata/delivery_validation.tsv',sep='\t',index=False)
 assert truth,name

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
a=ad.read_h5ad(root/'cell_annotation.h5ad');o=pd.read_parquet(root/'frozen_cell_labels.parquet');ex=pd.read_parquet(root/'excluded_cell_metadata.parquet');q=json.loads((root/'metadata/reconstruction_qc_summary.json').read_text());manifest=json.loads((root/'metadata/input_manifest.json').read_text());t=json.loads((root/'metadata/Cervilla_transform.json').read_text())
required='dataset_id section_id donor_id cell_id source_cell_id unit_type x_um y_um n_counts n_genes n_bins area_um2 segmentation_source retained_cell qc_flags reconstruction_flag low_information model_eligible model_role major_label subtype_label individual_prediction prediction_score prediction_margin annotation_status label_evidence unresolved_reason reference_version annotation_version'.split()
check('required_fields',set(required)<=set(o.columns))
check('ids_unique_aligned',a.obs_names.is_unique and o.index.is_unique and a.obs_names.equals(o.index) and np.array_equal(o.cell_id,o.index))
pd.testing.assert_frame_equal(a.obs.astype({c:'object' for c in a.obs.select_dtypes('category').columns}),o.astype({c:'object' for c in o.select_dtypes('category').columns}),check_dtype=False)
check('full_metadata_matches_h5ad',True)
check('main_sparse_raw_integer',sparse.issparse(a.X) and np.issubdtype(a.X.dtype,np.integer) and (a.X.data>=0).all(),str(a.shape)+' '+str(a.X.dtype))
check('raw_count_and_gene_QC',np.array_equal(a.X.sum(1).A1,o.n_counts) and np.array_equal(a.X.getnnz(1),o.n_genes))
check('raw_retention_floor',((o.n_counts>=10)&(o.n_genes>=10)&o.retained_cell).all())
check('low_information_exact',np.array_equal(o.low_information,o.n_genes.between(10,30)),int(o.low_information.sum()))
check('single_actual_section',o.section_id.nunique()==1 and o.dataset_id.eq('Cervilla').all())
check('uncertainty_retained',len(o)==q['retained'] and o.annotation_status.eq('Unknown').sum()==50767 and len(o)+len(ex)==q['reconstructed'])
check('excluded_not_in_main',len(o.index.intersection(ex.index))==0 and (~ex.retained_cell).all())
check('all_input_genes_preserved',a.n_vars==18085 and a.var_names.is_unique)
check('latent_finite_and_roles',np.isfinite(a.obsm['X_scVI']).all() and o.loc[o.low_information,'model_role'].eq('inference_only').all() and not o.loc[o.low_information,'model_eligible'].any())
check('UMAP_only_core',np.isfinite(a.obsm['X_umap'][o.scvi_training_core]).all() and np.isnan(a.obsm['X_umap'][~o.scvi_training_core]).all(),int(o.scvi_training_core.sum()))
for inp in manifest:check('unchanged_input_sha256_'+inp['dataset_id'],sha(inp['matrix'])==inp['matrix_sha256'],inp['matrix'])
bm=pd.read_parquet(root/'inputs/Cervilla_bin_to_cell.parquet');p=Path(manifest[0]['matrix']);print('Read original counts for full independent aggregation',flush=True)
with h5py.File(p) as f:
 g=f['matrix'];bar=g['barcodes'].asstr()[:];ids=g['features/id'].asstr()[:];shape=g['shape'][:];raw=sparse.csr_matrix((g['data'][:].astype('int32'),g['indices'][:],g['indptr'][:]),shape=(shape[1],shape[0]))
check('mapping_original_barcode_order',np.array_equal(bar,bm.barcode) and np.array_equal(bm.source_bin_index,np.arange(len(bar))))
check('original_feature_identity',np.array_equal(ids,a.var_names))
check('bin_raw_count_metadata',np.array_equal(raw.sum(1).A1,bm.raw_counts))
check('uncovered_bins_not_assigned',(bm.loc[~bm.image_covered,'cell_label']==0).all())
assigned=bm.cell_label.to_numpy()>0;labels=bm.cell_label.to_numpy();objects=np.unique(labels[assigned]);joined=pd.concat([o,ex]).set_index('source_cell_id').loc[objects.astype(str)];ix=np.searchsorted(objects,labels[assigned]);assignment=sparse.csr_matrix((np.ones(assigned.sum(),dtype='int32'),(ix,np.flatnonzero(assigned))),shape=(len(objects),len(bar)));cells=(assignment@raw).tocsr();cells.eliminate_zeros()
check('every_cell_mapping_present',len(objects)==len(joined)==q['reconstructed'])
check('every_cell_raw_qc_reproduced',np.array_equal(cells.sum(1).A1,joined.n_counts) and np.array_equal(cells.getnnz(1),joined.n_genes))
sel=np.searchsorted(objects,o.source_cell_id.astype(int));delta=cells[sel]-a.X;delta.eliminate_zeros();check('full_raw_cell_matrix_reproduced',delta.nnz==0,'Every retained cell/gene exactly equals independent original-bin SUM')
check('per_gene_conservation',np.array_equal(raw[assigned].sum(0),cells.sum(0)))
check('assigned_unassigned_count_conservation',int(raw.sum())==int(cells.sum())+int(raw[~assigned].sum()),int(raw.sum()))
check('cell_bin_counts_and_area',np.array_equal(assignment.sum(1).A1,joined.n_bins) and np.array_equal(joined.n_bins*4,joined.area_um2))
pos=pd.read_parquet(p.parent/'spatial/tissue_positions.parquet',filters=[('in_tissue','==',1)]).set_index('barcode').loc[bar];xy=pos[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();sf=json.loads((p.parent/'spatial/scalefactors_json.json').read_text());cent=(assignment@xy)/joined.n_bins.to_numpy()[:,None];check('all_cell_coordinates_reproduced',np.allclose(cent*sf['microns_per_pixel'],joined[['x_um','y_um']].to_numpy(),rtol=0,atol=1e-8))
check('main_image_coordinate_obsm',np.allclose(a.obsm['spatial'],cent[sel]) and np.allclose(a.obsm['spatial_um'],o[['x_um','y_um']]))
check('annotation_summary_exact',o.annotation_status.value_counts().to_dict()==json.loads((root/'metadata/annotation_summary.json').read_text())['annotation_status'])
# Reference model shared-gene mapping, source-group split and panel provenance.
conf=json.loads((root/'metadata/classifier_config.json').read_text());check('reference_model_genes_unique',len(conf['gene_list'])==len(set(conf['gene_list']))==14327)
train=pd.read_parquet(root/'metadata/reference_training_selection.parquet');hold=pd.read_parquet(root/'metadata/reference_holdout_selection.parquet');check('training_holdout_cell_disjoint',len(train.index.intersection(hold.index))==0)
print('SELECTION_COLUMNS',train.columns.tolist(),flush=True)
for c in ['source_group','validation_group','patient_group']:
 if c in train and c in hold:check('training_holdout_group_disjoint',not set(train[c])&set(hold[c]),c)
print('VALIDATION COMPLETE',len(rows),flush=True)
