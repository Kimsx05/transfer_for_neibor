from common5 import *
def run(s):
 out=R/s;d=pd.read_parquet(out/'all_cells.parquet');old=pd.read_parquet(out/'historical_candidates.parquet');v4=pd.read_parquet(out/'old_model_candidates.parquet');new=d.loc[d.T_hint,['id','x','y']].copy();new['source']='new_cell_hint'
 raw=pd.read_parquet(out/'raw_candidates.parquet');assign=pd.read_parquet(out/'assignments.parquet');raw['id']=assign.new_id.to_numpy()[raw.source_bin_index];raw['source']='raw2';raw['unit']='2:'+raw.source_bin_index.astype(str)
 native=pd.read_parquet(V/s/'bin_evidence.parquet',columns=['bin_id','array_row','array_col','component','x_um','y_um','type_blca','type_hpca','eligible_pool']);nb=native[native[['type_blca','type_hpca']].eq('T/NK').any(axis=1)].rename(columns={'x_um':'x','y_um':'y','bin_id':'historical_unit_id'}).copy();nb['source']='native8_T';old=pd.concat([old,v4,new,nb],ignore_index=True)
 old['unit']='n:'+old.id.astype('Int64').astype(str)
 for prefix,src in [('r:','historical_T'),('8:','native8_T'),('8grid:','v4_T8'),('16grid:','v4_T16')]:
  pick=old.source.str.startswith(src);old.loc[pick,'unit']=prefix+old.loc[pick,'historical_unit_id'].astype(str)
 cand=pd.concat([old[['source','unit','id','x','y']],raw[['source','unit','id','x','y']]],ignore_index=True);fl=cand.groupby('unit',sort=False).source.agg(lambda v:'|'.join(sorted(set(v))));cand=cand.drop_duplicates('unit').copy();cand['sources']=cand.unit.map(fl)
 tree=cKDTree(cand[['x','y']]);xy=d[['x','y']].to_numpy();d['in_D']=tree.query(xy,distance_upper_bound=200,workers=4)[0]<=200
 # Radius union grouped on a regular physical lattice. Exact distance determines D membership.
 origin=np.floor((cand[['x','y']].min().to_numpy()-224)/8)*8;nx,ny=(np.ceil((cand[['x','y']].max().to_numpy()+224-origin)/8).astype(int)+1)
 gx,gy=np.meshgrid(np.arange(nx),np.arange(ny));mx=np.c_[gx.ravel(),gy.ravel()]*8+origin
 domain=(tree.query(mx,workers=4,distance_upper_bound=200.000001)[0]<=200).reshape(ny,nx)
 regions,n=ndimage.label(domain,structure=np.ones((3,3)));np.savez_compressed(out/'D_geometry.npz',roi=regions.astype('uint16'),origin=origin,step_um=8)
 def locate(pos):
  ij=np.rint((pos-origin)/8).astype(int);ij=np.clip(ij,[0,0],[nx-1,ny-1]);return regions[ij[:,1],ij[:,0]]
 cand['roi']=locate(cand[['x','y']].to_numpy());assert cand.roi.min()>0
 d['roi']=locate(xy);zero=d.in_D&d.roi.eq(0)
 if zero.any():d.loc[zero,'roi']=cand.roi.to_numpy()[tree.query(xy[zero])[1]]
 # Component is determined at the nucleus center, from the original supplied native grid.
 au=js(N/s/'metadata/input_audit.json');aff=np.asarray(au['array_col_row_affine_fullres']);rc=((xy/au['native_mpp']-aff[2])@np.linalg.inv(aff[:2]));nrc=np.floor((rc+.5)/4).astype(int)
 ni=pd.MultiIndex.from_arrays([native.array_row,native.array_col]);ii=ni.get_indexer(pd.MultiIndex.from_arrays([nrc[:,1],nrc[:,0]]));d['component']=np.where(ii>=0,native.component.to_numpy()[np.maximum(0,ii)],-1)
 # Multiple target cells per historical region, based on its actual member 8um bins.
 parent=(assign.array_row.to_numpy()//4).astype('int64')*100000+(assign.array_col.to_numpy()//4);par=pd.DataFrame(dict(key=parent,id=assign.new_id));par=par[par.id>0].drop_duplicates()
 nk=native.array_row.astype('int64')*100000+native.array_col;native['key']=nk
 targets=[]
 for size in ['T16','T32']:
  wanted=old[old.source.eq('historical_'+size)].historical_unit_id
  mm=pd.read_parquet(V/s/(size+'_members.parquet'));mm=mm[mm.region_id.isin(wanted)].merge(native[['bin_id','key']],on='bin_id').merge(par,on='key')
  mm['unit']='r:'+mm.region_id;targets.append(mm[['unit','id']].drop_duplicates())
 nn=nb.assign(key=nb.array_row.astype('int64')*100000+nb.array_col).merge(par,on='key')
 nn['unit']='8:'+nn.historical_unit_id;targets.append(nn[['unit','id']])
 for k,f,prefix in [('T8_from2',4,'8grid:'),('T16_from2',8,'16grid:')]:
  m=pd.read_parquet(U/s/(k+'_members.parquet'),columns=['id','source_bin_index']);m['unit']=prefix+m.id;m=m[m.unit.isin(cand.unit)];m['id']=assign.new_id.to_numpy()[m.source_bin_index];targets.append(m[m.id.gt(0)][['unit','id']].drop_duplicates())
 targets.append(cand[cand.id.gt(0)][['unit','id']]);links=pd.concat(targets,ignore_index=True).drop_duplicates();links['id']=links.id.astype(int);links=links.merge(d[['id','status']],on='id',validate='many_to_one')
 links.to_parquet(out/'candidate_cell_links.parquet',index=False);summ=links.groupby('unit').status.agg(lambda v:'|'.join(sorted(set(v))));cand['disposition']=cand.unit.map(summ).fillna('no_assignable_cell');cand.to_parquet(out/'candidate_ledger.parquet',index=False);d.to_parquet(out/'cells_pre_MP10.parquet',index=False)
 counts=cand.groupby('roi').agg(candidates=('unit','size'),unassigned=('disposition',lambda v:v.eq('no_assignable_cell').sum()),xmin=('x','min'),xmax=('x','max'),ymin=('y','min'),ymax=('y','max')).reset_index();counts['section']=s;put(counts,out/'ROI_manifest.tsv');print(s,len(cand),counts.shape,d.status.value_counts().to_dict(),flush=True)
if __name__=='__main__':
 for s in sections:run(s)
