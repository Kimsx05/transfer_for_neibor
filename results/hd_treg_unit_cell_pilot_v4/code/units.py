from base import *
from scipy.spatial import cKDTree
import anndata as ad

def read_rows(g, rows):
 """Read contiguous HDF5 slabs only for blocks containing requested barcodes."""
 ptr=g['indptr'][:]; genes=len(g['features/id']);parts=[];barcode=[]
 for group in np.split(rows,np.flatnonzero(np.diff(rows//4096))+1):
  a,b=int(group[0]),int(group[-1])+1;lo,hi=ptr[a],ptr[b]
  mat=sparse.csr_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr[a:b+1]-lo),shape=(b-a,genes))
  parts.append(mat[group-a]);barcode.extend(g['barcodes'][a:b][group-a].astype(str))
 return sparse.vstack(parts,format='csr').astype('int64'),barcode

def geometry(s):
 root=R/s;roi=rd(R/'roi_manifest.tsv').query('section==@s');audit=js(N/s/'metadata/input_audit.json');mpp=audit['native_mpp'];aff=np.array(audit['array_col_row_affine_fullres'])
 cells=pd.read_parquet(S/s/'object_qc.parquet',columns=['cell_id','original_nucleus_id','nucleus_x_um','nucleus_y_um','image_confirmed_invalid','geometry_review']).reset_index().rename(columns={'original_nucleus_id':'nid','nucleus_x_um':'x_um','nucleus_y_um':'y_um'});cells['roi']=roi_for(cells.x_um,cells.y_um,roi);cells=cells[cells.roi.ne('')].reset_index(drop=True)
 source=S/s/'bin_to_cell.parquet';filters=[]
 for q in roi.itertuples():
  filters.append([('pxl_col_in_fullres','>=',(q.x0-20)/mpp),('pxl_col_in_fullres','<=',(q.x1+20)/mpp),('pxl_row_in_fullres','>=',(q.y0-20)/mpp),('pxl_row_in_fullres','<=',(q.y1+20)/mpp)])
 d=pd.read_parquet(source,filters=filters)
 # All cell members irrespective of location. Integer labels avoid assumptions about nucleus radii.
 extra=pd.read_parquet(source,filters=[('cell_label','in',cells.nid.tolist())]);d=pd.concat([d,extra]).drop_duplicates('source_bin_index').sort_values('source_bin_index').reset_index(drop=True)
 d['x_um']=d.pxl_col_in_fullres*mpp;d['y_um']=d.pxl_row_in_fullres*mpp
 with h5py.File(audit['input_matrix'])as f:
  g=f['matrix'];raw,bc=read_rows(g,d.source_bin_index.to_numpy());assert np.array_equal(bc,d.barcode)
  genes=g['features/id'][:].astype(str);symbols=g['features/name'][:].astype(str)
 assert np.array_equal(np.asarray(raw.sum(1)).ravel(),d.raw_counts);put(pd.DataFrame({'gene_id':genes,'symbol':symbols}),root/'genes.tsv')
 sel=(d.in_tissue.eq(1)&~d.status.str.startswith('excluded')).to_numpy();native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['bin_id','array_row','array_col','x_um','y_um','component','type_blca','type_hpca','eligible_pool']);native=native.set_index(['array_row','array_col']);parent=native.reindex(pd.MultiIndex.from_arrays([d.array_row//4,d.array_col//4])).reset_index(drop=True)
 d['parent8']=parent.bin_id;d['component']=parent.component;d['input_usable']=sel
 def grid_xy(factor):
  rc=d.array_row//factor;cc=d.array_col//factor
  xy=np.c_[cc*factor+(factor-1)/2,rc*factor+(factor-1)/2,np.ones(len(d))]@aff*mpp
  return rc,cc,xy
 _,_,xy=grid_xy(4);err=np.hypot(xy[:,0]-parent.x_um,xy[:,1]-parent.y_um);assert np.nanquantile(err,.99)<1, np.nanquantile(err,.99)
 d['parent_position_error_um']=err
 for g in T+core+epi+fib+mye:d[g]=np.asarray(raw[:,symbols==g].sum(1)).ravel()
 d.to_parquet(root/'sources.parquet',index=False);summary=[]
 for kind in ['T8_from2','T16_from2','cell_nucleus','cell_conservative','existing_cell8']:
  if kind.startswith('T'):
   f=4 if kind=='T8_from2' else 8;rr,cc,xy=grid_xy(f);ids=rr.astype(str)+'_'+cc.astype(str);rois=roi_for(xy[:,0],xy[:,1],roi)
   obj=pd.DataFrame(dict(id=ids,x_um=xy[:,0],y_um=xy[:,1],roi=rois)).loc[rois!=''].drop_duplicates('id').reset_index(drop=True)
  else:
   key={'cell_nucleus':'nucleus_label','cell_conservative':'v1_cell_label','existing_cell8':'cell_label'}[kind];ids=d[key].astype(str);obj=cells.copy();obj['id']=obj.nid.astype(str)
  idx=pd.Index(obj.id).get_indexer(ids);ok=(idx>=0)&sel;aggregate=sparse.csr_matrix((np.ones(ok.sum(),int),(idx[ok],np.flatnonzero(ok))),shape=(len(obj),len(d)))@raw
  residual=np.asarray(raw.sum(0)-aggregate.sum(0)-raw[~ok].sum(0));assert not residual.any()
  obj['branch']=kind;obj['section']=s;obj['n_supplied_usable']=np.bincount(idx[ok],minlength=len(obj));obj['nCount']=np.asarray(aggregate.sum(1)).ravel();obj['nFeature']=aggregate.getnnz(1)
  if kind in ['cell_conservative','existing_cell8']:
   old=ad.read_h5ad((N if kind=='cell_conservative' else S)/s/'cells.h5ad',backed='r');ix=old.obs_names.get_indexer(obj.cell_id);assert min(ix)>=0 and np.array_equal(old.var_names,genes);hist=sparse.csr_matrix(old.X[ix,:]);old.file.close();assert (hist!=aggregate).nnz==0;aggregate=hist # reuse verified history
  for g in T+core+epi+fib+mye:obj[g]=np.asarray(aggregate[:,symbols==g].sum(1)).ravel()
  for label,panel in [('T',T),('Epi',epi),('Fib',fib),('Myeloid',mye)]:obj[label+'_UMI']=obj[panel].sum(1);obj[label+'_width']=obj[panel].gt(0).sum(1)
  obj['FOXP3_raw_screen']=obj.FOXP3.gt(0);obj['FOXP3_T_screen']=obj.FOXP3.gt(0)&obj.T_UMI.gt(0)
  mem=d.loc[ok,['barcode','source_bin_index','parent8','component']].copy();mem['id']=obj.id.to_numpy()[idx[ok]];mem.to_parquet(root/(kind+'_members.parquet'),index=False)
  obj.to_parquet(root/(kind+'.parquet'),index=False);sparse.save_npz(root/(kind+'.npz'),aggregate)
  summary.append(dict(branch=kind,section=s,read_scope=len(d),all_UMI=int(raw.sum()),object_UMI=int(aggregate.sum()),excluded_UMI=int(raw[~sel].sum()),outside_selected_objects_or_unassigned_UMI=int(raw[sel&(idx<0)].sum()),max_per_gene_residual=0,n=len(obj),median_UMI=obj.nCount.median(),FOXP3_objects=int(obj.FOXP3_raw_screen.sum()),T_and_FOXP3=int(obj.FOXP3_T_screen.sum()),parent8_coordinate_error_p99=float(np.nanquantile(err,.99))))
 put(pd.DataFrame(summary),root/'count_accounting.tsv')
if __name__=='__main__':
 for s in sections:geometry(s)
