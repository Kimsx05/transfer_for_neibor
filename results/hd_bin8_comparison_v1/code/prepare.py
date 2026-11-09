from pathlib import Path
import numpy as np,pandas as pd,h5py,json,sys,os,gc
from scipy import sparse
root=Path(__file__).resolve().parents[1];hd=root.parents[1];s=sys.argv[1];o=root/s;o.mkdir(exist_ok=True)
if (o/'input_summary.json').exists():print('INPUT_CACHED',s);sys.exit()
v2=hd/'analysis/segmentation_assignment_v2'/s;old=hd/'analysis/hd_major_annotation_v1'/s;diag=hd/'analysis/hd_annotation_diagnosis_singler_v1'/s
b2=hd/('data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um' if s=='Cervilla' else 'data/GSE300692/D1_Johnson/outs/binned_outputs/square_002um');b8=o/'native8' if s=='Cervilla' else hd/'data/GSE300692/D1_Johnson/outs/binned_outputs/square_008um'
def read10(p):
 with h5py.File(p) as f:
  g=f['matrix'];shape=g['shape'][:];x=sparse.csr_matrix((g['data'][:].astype(np.int32),g['indices'][:],g['indptr'][:]),shape=(shape[1],shape[0]));bc=g['barcodes'].asstr()[:];ids=g['features/id'].asstr()[:];sy=g['features/name'].asstr()[:]
 x.eliminate_zeros();return x,bc,ids,sy
raw,bc,ids,sy=read10(b2/'filtered_feature_bc_matrix.h5');mp=pd.read_parquet(v2/'bin_to_cell.parquet');assert np.array_equal(bc,mp.barcode);assert np.array_equal(np.asarray(raw.sum(1)).ravel(),mp.raw_counts)
roi=pd.read_csv(old/'spatial_assets/rois.tsv',sep='\t');roi.to_csv(o/'fixed_rois.tsv',sep='\t',index=False);tr=json.loads((v2/'metadata/transform.json').read_text());(o/'transform.json').write_text(json.dumps(tr,indent=2));cr=tr['fullres_crop_xyxy'];sc=tr['local_to_model_scale_xy'];mpp=tr['native_mpp'];mx=(mp.pxl_col_in_fullres-cr[0]+.5)*sc[0]-.5;my=(mp.pxl_row_in_fullres-cr[1]+.5)*sc[1]-.5
rn=np.full(len(mp),'',object)
for r in roi.itertuples():rn[(mx>=r.x)&(mx<r.X)&(my>=r.y)&(my<r.Y)]=r.roi
pan=pd.read_csv(old/'tables/marker_panel_genes.tsv',sep='\t');print('PANELS',pan.columns.tolist(),flush=True)
# Exact original panels, no new lineage thresholds.
panels={k:g.symbol.tolist() for k,g in pan.groupby('major')}
keys=(mp.array_row.values//4)*1000+mp.array_col.values//4;uk,rev=np.unique(keys,return_inverse=True);A=sparse.csr_matrix((np.ones(len(mp),np.int32),(rev,np.arange(len(mp)))),shape=(len(uk),len(mp)));summed=(A@raw).tocsr();assert np.array_equal(np.asarray(summed.sum(0)).ravel(),np.asarray(raw.sum(0)).ravel())
native=(b8/'filtered_feature_bc_matrix.h5').exists()
if native:
 x,bar,gi,gs=read10(b8/'filtered_feature_bc_matrix.h5');assert np.array_equal(ids,gi);pos=pd.read_parquet(b8/'spatial/tissue_positions.parquet').set_index('barcode').loc[bar].reset_index();ky=pos.array_row.values*1000+pos.array_col.values
else:
 assert s=='Cervilla';x=summed;ky=uk;bar=np.array([f's_008um_{k//1000:05d}_{k%1000:05d}-1' for k in ky]);sample=np.linspace(0,len(mp)-1,2000,dtype=int);coef=np.linalg.lstsq(np.c_[mp.array_row.values[sample],mp.array_col.values[sample],np.ones(len(sample))],mp[['pxl_col_in_fullres','pxl_row_in_fullres']].values[sample],rcond=None)[0];xy=np.c_[(ky//1000)*4+1.5,(ky%1000)*4+1.5,np.ones(len(ky))]@coef;pos=pd.DataFrame(dict(barcode=bar,array_row=ky//1000,array_col=ky%1000,pxl_col_in_fullres=xy[:,0],pxl_row_in_fullres=xy[:,1],in_tissue=np.asarray(A@mp.in_tissue.values).ravel()>0));b8.mkdir(exist_ok=True)
 with h5py.File(b8/'filtered_feature_bc_matrix.h5','w') as f:
  g=f.create_group('matrix');g['shape']=[len(ids),len(bar)]
  for k,v in [('data',x.data),('indices',x.indices),('indptr',x.indptr),('barcodes',bar.astype('S'))]:g.create_dataset(k,data=v,compression='gzip')
  gg=g.create_group('features');gg['id']=ids.astype('S');gg['name']=sy.astype('S');gg['feature_type']=np.array(['Gene Expression']*len(ids),dtype='S')
 idx=pd.Index(uk).get_indexer(ky)
idx=pd.Index(uk).get_indexer(ky);inside=idx>=0;meta=pos.copy();meta['section_id']=s;meta['bin_id']=bar;meta['nCount_RNA']=np.asarray(x.sum(1)).ravel();meta['nFeature_RNA']=x.getnnz(1);meta['qc_keep']=(meta.in_tissue==1)&(meta.nCount_RNA>5);meta['qc_reason']=np.where(meta.in_tissue!=1,'outside_tissue',np.where(meta.nCount_RNA<=5,'nCount_le5','pass'))
meta['x_um']=meta.pxl_col_in_fullres*mpp;meta['y_um']=meta.pxl_row_in_fullres*mpp
for field,values in [('supplied_2um_bins',np.ones(len(mp),np.int32)),('common_2um_UMI',mp.raw_counts.values),('cell_assigned_UMI',mp.raw_counts.values*(mp.cell_label.values>0)),('cell_unassigned_UMI',mp.raw_counts.values*((mp.cell_label.values==0)&mp.status.str.startswith('unassigned').values)),('cell_excluded_UMI',mp.raw_counts.values*mp.status.str.startswith('excluded').values)]:
 vals=np.asarray(A@values).ravel();v=np.zeros(len(ky),np.int64);v[inside]=vals[idx[inside]];meta[field]=v
qc=pd.read_parquet(v2/'object_qc.parquet');labels=qc.object_id.values;keepmap=np.zeros(labels.max()+1,bool);keepmap[labels]=qc.qc_keep_analysis.values;u=mp.raw_counts.values*keepmap[mp.cell_label.values];val=np.asarray(A@u).ravel();meta['cell_analysis_UMI']=np.where(inside,val[np.maximum(idx,0)],0);meta['coverage_fraction_2um']=meta.supplied_2um_bins/16;meta['partial_2um_coverage']=(meta.supplied_2um_bins>0)&(meta.supplied_2um_bins<16);meta['native_minus_supplied_2um_UMI']=meta.nCount_RNA-meta.common_2um_UMI
# Grid coordinate relationship, verified from actual shared centers, never image dimensions.
ix=np.linspace(0,len(mp)-1,2000,dtype=int);coef=np.linalg.lstsq(np.c_[mp.array_row.values[ix],mp.array_col.values[ix],np.ones(len(ix))],mp[['pxl_col_in_fullres','pxl_row_in_fullres']].values[ix],rcond=None)[0];fit=np.c_[pos.array_row.values*4+1.5,pos.array_col.values*4+1.5,np.ones(len(pos))]@coef;err=np.max(np.abs(fit-pos[['pxl_col_in_fullres','pxl_row_in_fullres']].values));assert err*mpp<.1,err
meta['roi']='';bx=(meta.pxl_col_in_fullres-cr[0]+.5)*sc[0]-.5;by=(meta.pxl_row_in_fullres-cr[1]+.5)*sc[1]-.5
for r in roi.itertuples():meta.loc[(bx>=r.x)&(bx<r.X)&(by>=r.y)&(by<r.Y),'roi']=r.roi
meta.to_csv(o/'bin_metadata.tsv.gz',sep='\t',index=False);pd.DataFrame(dict(feature_id=ids,symbol=sy)).to_csv(o/'gene_mapping.tsv',sep='\t',index=False)
hit=rn!='';r2=mp.loc[hit,['barcode','array_row','array_col','raw_counts','cell_label','v1_cell_label','status']].copy();r2['roi']=rn[hit];r2['x_um']=mp.pxl_col_in_fullres.values[hit]*mpp;r2['y_um']=mp.pxl_row_in_fullres.values[hit]*mpp;r2['bin8_id']=[f's_008um_{k//1000:05d}_{k%1000:05d}-1' for k in keys[hit]];rawroi=raw[hit];coverage=[]
for k,genes in panels.items():
 ii=np.flatnonzero(np.isin(sy,genes));r2['raw_'+k]=np.asarray(rawroi[:,ii].sum(1)).ravel();r2['detected_'+k]=rawroi[:,ii].getnnz(1);coverage.append(dict(major=k,requested=len(genes),present=','.join(g for g in genes if g in sy),missing=','.join(g for g in genes if g not in sy)))
r2.to_parquet(o/'roi_2um_signals.parquet',compression='zstd',index=False);pd.DataFrame(coverage).to_csv(o/'marker_coverage.tsv',sep='\t',index=False)
# All supplied 2um->bin8 relations derive from row//4,col//4; no nearest-centroid pairing.
commons=np.isin(uk,ky);common2genes=np.asarray(summed[commons].sum(0)).ravel();nativecommon=np.asarray(x[inside].sum(0)).ravel();pd.DataFrame(dict(feature_id=ids,all_supplied_2um=np.asarray(raw.sum(0)).ravel(),shared_grid_supplied_2um=common2genes,native8_shared_grid=nativecommon,mask_difference=nativecommon-common2genes)).to_csv(o/'gene_count_accounting.tsv',sep='\t',index=False)
info=dict(section=s,input_kind='native_square_008um' if native else 'supplied_2um_fixed_grid_4x4',matrix_path=str((b8/'filtered_feature_bc_matrix.h5').relative_to(hd)),genes=len(ids),before_bins=len(meta),before_UMI=int(meta.nCount_RNA.sum()),QC_bins=int(meta.qc_keep.sum()),QC_UMI=int(meta.loc[meta.qc_keep,'nCount_RNA'].sum()),supplied2_UMI=int(raw.sum()),shared_grid_2um_UMI=int(common2genes.sum()),supplied2_outside_8grid_UMI=int(summed[~commons].sum()),native8_shared_grid_UMI=int(nativecommon.sum()),native8_only_grid_UMI=int(x[~inside].sum()),native_minus_supplied_mask_UMI=int(meta.native_minus_supplied_2um_UMI.sum()),QC_common2_UMI=int(meta.loc[meta.qc_keep,'common_2um_UMI'].sum()),QC_cell_unassigned_UMI=int(meta.loc[meta.qc_keep,'cell_unassigned_UMI'].sum()),QC_cell_analysis_UMI=int(meta.loc[meta.qc_keep,'cell_analysis_UMI'].sum()),partial_coverage_bins=int(meta.partial_2um_coverage.sum()),grid_center_max_error_pixels=float(err),per_gene_2um_aggregation_conserved=True,grid_origin=[0,0],physical_mpp=mpp,normalization='not yet; raw counts only',fixed_roi_names=roi.roi.tolist())
(o/'input_summary.json').write_text(json.dumps(info,indent=2));print(info,flush=True)
