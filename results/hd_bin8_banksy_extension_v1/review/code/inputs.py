"""Read original supplied matrices; fixed-origin aggregation. Night-only entry in rerun.sh."""
from pathlib import Path
import sys,json,re,h5py,numpy as np,pandas as pd
from scipy import sparse,ndimage
r=Path(__file__).resolve().parents[1];hd=r.parents[1];s=sys.argv[1];o=r/s;o.mkdir(exist_ok=True)
paths={'UC_SCC':'data/GSE299573/GSM9040999/square_008um/filtered_feature_bc_matrix.h5','MPBC1':'data/GSE318916/GSM9505917_mpbc4423filtered_feature_bc_matrix.h5','MPBC2':'data/GSE318916/GSM9505918_mpbc4538filtered_feature_bc_matrix.h5'}
f=hd/paths[s];sig=dict(path=paths[s],bytes=f.stat().st_size,mtime_ns=f.stat().st_mtime_ns)
if (o/'input_summary.json').exists():
 assert json.loads((o/'input_summary.json').read_text())['source_signature']==sig
 print('INPUT_CACHE',s);sys.exit()
with h5py.File(f) as h:
 g=h['matrix'];assert np.issubdtype(g['data'].dtype,np.integer);sh=tuple(int(v) for v in g['shape'][:]);bc=g['barcodes'][:].astype('U');ids=g['features/id'][:].astype('U');sym=g['features/name'][:].astype('U');ft=g['features/feature_type'][:].astype('U');assert np.all(ft=='Gene Expression') and len(set(ids))==len(ids)
 # Fixed ASCII fields inspected first; full vector format and unique grid keys then checked.
 a=pd.Series(bc).str.extract(r'^s_(\d{3})um_(\d{5})_(\d{5})-(\d+)$');assert a.notna().all().all();a=a.astype(np.int32);res,row,col,lib=(a[i].to_numpy() for i in range(4));assert len(np.unique(lib))==1 and np.all(res==(8 if s=='UC_SCC' else 2));assert row.min()>=0 and col.min()>=0
 key=(lib.astype(np.int64)*100000+row)*100000+col;assert np.unique(key).size==len(key)
 mat=sparse.csc_matrix((g['data'][:].astype(np.int64),g['indices'][:],g['indptr'][:]),shape=sh);assert np.all(mat.data>=0)
 info=dict(software_version=str(h.attrs.get('software_version','')),chemistry=str(h.attrs.get('chemistry_description','')),library_ids=[x.decode() for x in h.attrs.get('library_ids',[])])
mat.eliminate_zeros()
original_gene=np.asarray(mat.sum(axis=1)).ravel();input_umi=int(original_gene.sum());zero=np.diff(mat.indptr)==0
if s!='UC_SCC':
 row8=row//4;col8=col//4;parent=(lib.astype(np.int64)*100000+row8)*100000+col8;u,inv=np.unique(parent,return_inverse=True);rr=(u//100000)%100000;cc=u%100000;ll=u//100000**2
 H=sparse.csr_matrix((np.ones(len(inv),np.int8),(np.arange(len(inv)),inv)),shape=(len(inv),len(u)))
 b=(mat@H).tocsc();b.sum_duplicates();b.eliminate_zeros();del H,mat
 assert np.array_equal(original_gene,np.asarray(b.sum(axis=1)).ravel());bc8=np.array([f's_008um_{x:05d}_{y:05d}-{z}' for x,y,z in zip(rr,cc,ll)])
 n2=np.bincount(inv,minlength=len(u));z2=np.bincount(inv,weights=zero.astype(int),minlength=len(u)).astype(int)
 d=pd.DataFrame(dict(bin_id=bc8,array_row=rr,array_col=cc,array_id=ll,x_um=(cc+.5)*8,y_um=(rr+.5)*8,supplied_bin=True,in_tissue_original=np.nan,supplied_2um_n=n2,supplied_zero_2um_n=z2,unprovided_nominal_2um_positions=16-n2,nominal_4x4_supplied_fraction=n2/16,expected_array_clipped_positions=np.nan,coverage_state=np.where(n2==16,'full_nominal_4x4','partial_nominal_4x4_array_edge_unknown')))
 dest=o/'counts8.h5'
 with h5py.File(dest,'w') as h:
  q=h.create_group('matrix')
  for name,v in [('data',b.data),('indices',b.indices),('indptr',b.indptr),('shape',np.array(b.shape)),('barcodes',bc8.astype('S'))]:q.create_dataset(name,data=v,compression='gzip')
  fg=q.create_group('features')
  for name,v in [('id',ids),('name',sym),('feature_type',ft),('genome',np.repeat('GRCh38',len(ids)))]:fg.create_dataset(name,data=v.astype('S'),compression='gzip')
 counts_path='counts8.h5';transform=dict(coordinate_type='relative array centers; x=8*(col8+.5), y=8*(row8+.5); no image registration',origin='unchanged 0-based array origin',units='um',y_direction='increasing array row; anatomical orientation unknown',array_geometry='not independently supplied; 16 positions is nominal 4x4, not a validated clipped array boundary',in_tissue='unknown')
else:
 b=mat;d=pd.read_parquet(f.parent/'spatial/tissue_positions.parquet').set_index('barcode').loc[bc].reset_index().rename(columns={'barcode':'bin_id','in_tissue':'in_tissue_original'});assert np.array_equal(d.array_row,row) and np.array_equal(d.array_col,col)
 scales=json.loads((f.parent/'spatial/scalefactors_json.json').read_text());mpp=scales['microns_per_pixel'];assert scales['bin_size_um']==8
 d['x_um']=d.pxl_col_in_fullres*mpp;d['y_um']=d.pxl_row_in_fullres*mpp;d['supplied_bin']=True;d['array_id']=lib;d['coverage_state']='supplied_native8';d['supplied_2um_n']=np.nan;d['supplied_zero_2um_n']=np.nan
 counts_path=str(Path('../../../')/paths[s]);transform=dict(coordinate_type='native fullres pixel centers times declared microns_per_pixel',units='um',y_direction='down image',microns_per_pixel=mpp,scalefactors=scales,hires_path=str(Path('../../../')/f.parent.relative_to(hd)/'spatial/tissue_hires_image.png'),in_tissue='original barcode-aligned spatial table')
 # Affine relation checks row/col orientation and measured 8um lattice step without inventing metadata.
 V=np.column_stack([d.array_col,d.array_row,np.ones(len(d))]);coef=np.linalg.lstsq(V,d[['x_um','y_um']],rcond=None)[0];err=np.linalg.norm(V@coef-d[['x_um','y_um']],axis=1);transform['array_to_physical_affine']=coef.tolist();transform['affine_max_residual_um']=float(err.max());transform['step_col_row_um']=np.linalg.norm(coef[:2],axis=1).tolist();assert err.max()<.1
nc=np.asarray(b.sum(axis=0)).ravel();nf=np.diff(b.indptr);d['nCount_RNA']=nc;d['nFeature_RNA']=nf;d['qc_keep_expression']=nc>5;d['qc_keep']=d.qc_keep_expression if s!='UC_SCC' else d.qc_keep_expression&(d.in_tissue_original==1);d['qc_reason']=np.where(d.qc_keep,'kept_expression_only_tissue_unknown' if s!='UC_SCC' else 'kept',np.where(nc<=5,'nCount_le5','outside_original_tissue'));d['section_id']=s
assert int(nc.sum())==input_umi
pd.DataFrame(dict(feature_id=ids,symbol=sym,feature_type=ft)).to_csv(o/'gene_mapping.tsv',sep='\t',index=False)
pd.DataFrame(dict(feature_id=ids,input_counts=original_gene,aggregated_counts=np.asarray(b.sum(axis=1)).ravel())).to_csv(o/'gene_count_conservation.tsv',sep='\t',index=False)
d.to_csv(o/'bin_metadata.tsv.gz',sep='\t',index=False);d[d.qc_keep].to_csv(o/'analysis_metadata.tsv.gz',sep='\t',index=False)
(o/'transform.json').write_text(json.dumps(transform,indent=2));summary=dict(section=s,source_signature=sig,counts_relative=counts_path,input_resolution_um=int(res[0]),original_bins=sh[1],bins8=len(d),genes=sh[0],input_UMI=input_umi,QC_bins=int(d.qc_keep.sum()),QC_UMI=int(nc[d.qc_keep].sum()),excluded_UMI=int(nc[~d.qc_keep].sum()),gene_conservation=True,barcode_unique=True,array_id=np.unique(lib).tolist(),row2_or8_range=[int(row.min()),int(row.max())],col2_or8_range=[int(col.min()),int(col.max())],in_tissue_known=s=='UC_SCC',original_supplied_zero_bins=int(zero.sum()),missing_positions_are_zero=False,**info)
(o/'input_summary.json').write_text(json.dumps(summary,indent=2));print('INPUT_COMPLETE',summary,flush=True)
