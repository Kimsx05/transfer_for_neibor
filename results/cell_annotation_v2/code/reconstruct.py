"""Bin2cell expansion, original-count sparse aggregation, and permissive QC.
Raw integer X is restored before aggregation; adjusted counts are never rounded.
"""
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:12])
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,bin2cell as b2c,h5py,json,gc
from scipy import sparse
from common import ROOT,HD,section
import sys
name=sys.argv[1];meta,p,root=section(name);hd=HD;np.random.seed(123)
t=json.loads((root/'metadata/transform.json').read_text());sf=json.loads((p/'spatial/scalefactors_json.json').read_text())
assert not (root/'metadata/FROZEN.json').exists()
print('READ original integer bins',flush=True)
with h5py.File(p/'filtered_feature_bc_matrix.h5') as f:
 g=f['matrix'];bar=g['barcodes'].asstr()[:];ids=g['features/id'].asstr()[:];symbols=g['features/name'].asstr()[:];shape=tuple(g['shape'][:]);raw=sparse.csr_matrix((g['data'][:].astype('int32'),g['indices'][:],g['indptr'][:]),shape=(shape[1],shape[0]))
assert raw.has_canonical_format and raw.data.min()>=0
pos=pd.read_parquet(p/'spatial/tissue_positions.parquet',filters=[('in_tissue','==',1)]).set_index('barcode')
assert pos.index.is_unique
ix=pos.index.get_indexer(bar);assert (ix>=0).all();pos=pos.iloc[ix].copy();assert np.array_equal(pos.index,bar)
xy=pos[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();hires=xy-np.array(t['fullres_crop_xyxy'][:2]);mxy=np.rint((hires+.5)*t['local_to_model_scale_xy']-.5).astype('int32')
h,w=t['model_image_shape_yx'];inframe=(hires[:,0]>=0)&(hires[:,0]<t['source_image_shape_yx'][1])&(hires[:,1]>=0)&(hires[:,1]<t['source_image_shape_yx'][0])&(mxy[:,0]>=0)&(mxy[:,0]<w)&(mxy[:,1]>=0)&(mxy[:,1]<h)
mask=sparse.load_npz(root/'inputs/nuclei.npz');nseg=int(mask.max());nuclear_area=np.bincount(mask.data,minlength=nseg+1)*np.prod(t['actual_model_mpp_xy'])
labels=np.zeros(len(bar),dtype='int32');labels[inframe]=np.asarray(mask[mxy[inframe,1],mxy[inframe,0]]).ravel()
# Border labels have incomplete image context; retain as a review flag unless a severe error is verified.
border=np.unique(np.concatenate([mask[:2,:].data,mask[-2:,:].data,mask[:,:2].data,mask[:,-2:].data]));del mask
pos['labels_he']=labels;pos['n_counts']=np.asarray(raw.sum(1)).ravel();a=ad.AnnData(raw,obs=pos,var=pd.DataFrame({'gene_symbol':symbols},index=ids));a.obsm['spatial']=xy
print('DESTRIPE diagnostics',flush=True)
b2c.destripe(a,adjust_counts=False)
factor=np.divide(a.obs.n_counts_adjusted.to_numpy(),a.obs.n_counts.to_numpy(),out=np.ones(a.n_obs),where=a.obs.n_counts.to_numpy()>0);factor[~np.isfinite(factor)]=1
# Only the temporary tie-breaking expression view uses destriping.
a.X=sparse.diags(factor.astype('float32')).dot(raw)
original_pca=sc.pp.pca
def seeded_pca(*args,**kwargs):
 kwargs.setdefault('random_state',123);kwargs.setdefault('n_comps',30);return original_pca(*args,**kwargs)
sc.pp.pca=seeded_pca
print('EXPAND nuclei by <=2 array bins (4um)',flush=True)
b2c.expand_labels(a,labels_key='labels_he',expanded_labels_key='labels_cell',max_bin_distance=2,k=4,subset_pca=True)
sc.pp.pca=original_pca;a.X=raw;labels=a.obs.labels_cell.to_numpy().astype('int32');labels[~inframe]=0
assigned=labels>0;objects=np.unique(labels[assigned]);mapped=np.searchsorted(objects,labels[assigned]);assign=sparse.csr_matrix((np.ones(assigned.sum(),dtype='int32'),(mapped,np.flatnonzero(assigned))),shape=(len(objects),len(bar)))
print('AGGREGATE from original integer counts',len(objects),flush=True)
cellraw=(assign@raw).tocsr();cellraw.sum_duplicates();cellraw.eliminate_zeros();assert cellraw.data.min()>=0 and np.issubdtype(cellraw.dtype,np.integer)
expected=np.asarray(raw[assigned].sum(0)).ravel();actual=np.asarray(cellraw.sum(0)).ravel();assert np.array_equal(expected,actual)
nb=np.asarray(assign.sum(1)).ravel();centroid=(assign@xy)/nb[:,None]
obs=pd.DataFrame(index=[f'{name}__v2__cell_{i:07d}' for i in objects]);obs['cell_id']=obs.index;obs['source_cell_id']=objects.astype(str);obs['dataset_id']='Cervilla' if name=='Cervilla' else 'GSE300692';obs['section_id']=name;obs['donor_id']='unknown';obs['unit_type']='reconstructed_cell';obs['x_um']=centroid[:,0]*sf['microns_per_pixel'];obs['y_um']=centroid[:,1]*sf['microns_per_pixel'];obs['n_counts']=np.asarray(cellraw.sum(1)).ravel();obs['n_genes']=cellraw.getnnz(1);obs['n_bins']=nb;obs['area_um2']=nb*4;obs['nuclear_area_um2']=nuclear_area[objects];obs['segmentation_source']='H&E_StarDist_2D_versatile_he';obs['image_native_mpp']=t['native_mpp'];obs['coverage_status']='provided_counts_intersect_recovered_image'
obs['reconstruction_flag']='none_confirmed';obs['image_edge_nucleus']=np.isin(objects,border) # touching image border is a review flag, not proof of severe reconstruction error
flags=[[] for _ in objects]
for col in ['n_counts','n_genes','area_um2','nuclear_area_um2']:
 v=np.log1p(obs[col].to_numpy());med=np.median(v);mad=np.median(np.abs(v-med));hit=np.abs(v-med)>5*1.4826*mad if mad>0 else np.zeros(len(v),bool)
 for j in np.flatnonzero(hit):flags[j].append(col+'_5MAD')
obs['qc_flags']=[';'.join(x) if x else 'none' for x in flags];obs['low_information']=obs.n_genes.between(10,30);obs['retained_cell']=(obs.n_counts>=10)&(obs.n_genes>=10)&obs.reconstruction_flag.eq('none_confirmed');obs['model_eligible']=obs.retained_cell&(obs.n_genes>30)&obs.qc_flags.eq('none')&~obs.image_edge_nucleus;obs['model_role']=np.where(obs.model_eligible,'structure_core',np.where(obs.retained_cell,'prediction_only','excluded'))
# No morphology-outlier flag alone removes a cell. MT ratio is observed probe-panel fraction only.
mt=np.char.startswith(symbols.astype(str),'MT-');obs['mt_probe_count_fraction']=np.asarray(cellraw[:,mt].sum(1)).ravel()/np.maximum(obs.n_counts.to_numpy(),1)
obs['exclusion_reason']=np.where(~obs.reconstruction_flag.eq('none_confirmed'),'confirmed_image_edge_truncation',np.where((obs.n_counts<10)|(obs.n_genes<10),'below_raw_count_gene_floor',''))
obs.loc[~obs.retained_cell].to_parquet(root/'excluded_cell_metadata.parquet',index=True,compression='zstd')
# Include every source bin, with explicit unassigned/missing-image disposition, compactly.
map_df=pd.DataFrame({'source_bin_index':np.arange(len(bar),dtype='int32'),'barcode':bar,'nucleus_label':a.obs.labels_he.to_numpy().astype('int32'),'cell_label':labels,'image_covered':inframe,'raw_counts':a.obs.n_counts.to_numpy().astype('int32'),'array_row':pos.array_row.to_numpy().astype('int16'),'array_col':pos.array_col.to_numpy().astype('int16')})
map_df.to_parquet(root/'inputs/bin_to_cell.parquet',index=False,compression='zstd')
keep=obs.retained_cell.to_numpy();c=ad.AnnData(cellraw[keep],obs=obs.loc[keep].copy(),var=pd.DataFrame({'feature_id':ids,'gene_symbol':symbols,'assay':'Visium HD v1','probe_set':'Visium Human Transcriptome Probe Set v2.0'},index=ids));c.obsm['spatial']=centroid[keep];c.obsm['spatial_image']=centroid[keep]-np.array(t['fullres_crop_xyxy'][:2]);c.obsm['spatial_model']=((c.obsm['spatial_image']+.5)*t['local_to_model_scale_xy']-.5);c.obsm['spatial_um']=centroid[keep]*sf['microns_per_pixel'];c.uns.update(seed=123,analysis_status='reconstructed_not_annotated',counts_provenance='SUM of original Space Ranger integer counts by stored final bin-to-cell map; destriping used only for expansion tie breaking',input_matrix=str(p/'filtered_feature_bc_matrix.h5'),reference_version='cell2location_v2_frozen_20260916__scRNA32',coordinate_transform=json.dumps(t),mt_note='Measured MT probe counts / all counts; assay coverage incomplete and not used as a hard filter',coverage_status='provided_counts_intersect_recovered_image; image-only regions are not zero-expression tissue')
c.write_h5ad(root/'cells.h5ad',compression='gzip')
summary=dict(segmentation_nuclei=nseg,nuclei_with_assigned_bins=len(objects),nuclei_without_assigned_bins=nseg-len(objects),reconstructed=int((obs.n_counts>0).sum()),zero_expression_nuclei_with_bins=int((obs.n_counts==0).sum()),qc_input_objects=len(obs),retained=int(keep.sum()),excluded=int((~keep).sum()),low_information_retained=int(obs.loc[keep,'low_information'].sum()),structure_core=int(obs.model_eligible.sum()),source_bins=len(bar),image_uncovered_bins=int((~inframe).sum()),assigned_bins=int(assigned.sum()),unassigned_bins=int((~assigned).sum()),source_counts=int(raw.sum()),assigned_counts=int(cellraw.sum()),unassigned_counts=int(raw[~assigned].sum()),retained_counts=int(c.X.sum()),excluded_counts=int(cellraw[~keep].sum()),per_gene_count_conservation=True,mt_features=symbols[mt].tolist(),reconstruction='Bin2cell 0.3.5 max_bin_distance=2; seed123 PCA30; no GEX salvage in dense tissue; original counts aggregated')
assert summary['assigned_counts']+summary['unassigned_counts']==summary['source_counts']
(root/'metadata/reconstruction_qc_summary.json').write_text(json.dumps(summary,indent=2));obs.groupby(['retained_cell','exclusion_reason'],dropna=False).size().rename('n').reset_index().to_csv(root/'metadata/qc_disposition.tsv',sep='\t',index=False)
# No-bin nuclei are retained as compact geometry metadata, not called cells.
geom=np.load(root/'inputs/nucleus_geometry.npz');assert len(geom['points'])==nseg
absent=np.setdiff1d(np.arange(1,nseg+1),objects);points=geom['points'][absent-1];pd.DataFrame({'source_nucleus_label':absent,'model_y':points[:,0],'model_x':points[:,1],'reason':'no_assigned_in_tissue_2um_bins'}).to_parquet(root/'inputs/unassigned_nuclei.parquet',index=False,compression='zstd')
# Compact final cell assignment on the original 2um array grid, not an invented pixel boundary.
grid=sparse.coo_matrix((labels,(pos.array_row.to_numpy(),pos.array_col.to_numpy())),shape=(int(pos.array_row.max())+1,int(pos.array_col.max())+1)).tocsr();grid.eliminate_zeros();sparse.save_npz(root/'inputs/cell_mask_2um_array.npz',grid)
# Label fields are not yet attached; raw matrix contains every original feature exactly once.
print(json.dumps(summary,indent=2),flush=True)
