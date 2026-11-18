from common import *
from scipy.ndimage import distance_transform_edt,map_coordinates
rows=[];audit=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet');t=t[t.T_eligible].reset_index(drop=True);e=pd.read_parquet(R/s/'Epi_regions.parquet');tm=pd.read_parquet(V/s/'T32_members.parquet');em=pd.read_parquet(V/s/'Epi64_members.parquet');overlap=tm[tm.region_id.isin(t.region_id)].merge(em,on='bin_id',suffixes=('_T','_E'));bad=set(overlap.region_id_E);eok=e[~e.region_id.isin(bad)].copy();b=pd.read_parquet(V/s/'bin_evidence.parquet');b=b[b.eligible_pool].copy();audit.append(dict(section=s,T_eligible=len(t),Epi_before=len(e),Epi_global_source_disjoint=len(eok),Epi_removed=len(e)-len(eok),original_mask_known=s in SECTIONS[:3],window='supplied in-tissue8um observation union'if s in SECTIONS[:3]else'unknown tissue, supplied expression support only'))
 if s not in SECTIONS[:3]:
  for label in LABELS:
   for h in HIGHS:
    for r in range(0,201,10):rows.append(dict(section=s,component='all',label_definition=label,high_definition=h,distance_um=r,K=np.nan,L=np.nan,L_minus_r=np.nan,status='not_run_original_tissue_window_unknown_no_bounding_rectangle'))
  continue
 # Verify8um array is isometric to original physical plane before using a pixel window.
 v=b.iloc[::max(1,len(b)//10000)];X=np.column_stack([v.array_row,v.array_col,np.ones(len(v))]);coef=np.linalg.lstsq(X,v[['x_um','y_um']],rcond=None)[0];res=np.max(np.linalg.norm(X@coef-v[['x_um','y_um']].to_numpy(),axis=1));norm=np.linalg.norm(coef[:2],axis=1);ortho=abs(np.dot(coef[0],coef[1]))/np.prod(norm);isometric=res<.5 and np.max(abs(norm-8))<.08 and ortho<.01
 audit[-1].update(affine_max_residual_um=res,array_pitch_um=norm.tolist(),orthogonality_cos=ortho,isometric=isometric)
 bc=b.set_index('bin_id');coords=tm[tm.region_id.isin(t.region_id)].merge(b[['bin_id','array_row','array_col']],on='bin_id').groupby('region_id')[['array_row','array_col']].mean();t['row_mean']=coords.loc[t.region_id,'array_row'].to_numpy();t['col_mean']=coords.loc[t.region_id,'array_col'].to_numpy()
 for component,ts in t.groupby('component'):
  win=b[b.component.eq(component)];es=eok[eok.component.eq(component)];r0,c0=int(win.array_row.min()),int(win.array_col.min());mask=np.zeros((int(win.array_row.max())-r0+3,int(win.array_col.max())-c0+3),bool);mask[win.array_row.to_numpy(int)-r0+1,win.array_col.to_numpy(int)-c0+1]=True;dist=np.maximum(distance_transform_edt(mask)*8-4,0);edge=map_coordinates(dist,np.array([ts.row_mean.to_numpy()-r0+1,ts.col_mean.to_numpy()-c0+1]),order=1,mode='constant');area=int(mask.sum())*64
  for label in LABELS:
   aa=ts[ts[label]];aedge=edge[ts[label].to_numpy()]
   for h in HIGHS:
    bb=es[es[h]];near=cKDTree(bb[['x_um','y_um']]).query_ball_point(aa[['x_um','y_um']],200)if len(bb)else[[]for _ in range(len(aa))];dd=[np.linalg.norm(bb[['x_um','y_um']].to_numpy()[ix]-pos,axis=1)for ix,pos in zip(near,aa[['x_um','y_um']].to_numpy())]
    for r in range(0,201,10):
     inside=aedge>=r;na=int(inside.sum());pairs=sum(int((d<=r+1e-8).sum())for d,k in zip(dd,inside)if k);ok=isometric and len(aa)>=10 and len(bb)>=2 and na>=5;K=area*pairs/(na*len(bb))if ok else np.nan;L=np.sqrt(K/np.pi)if ok else np.nan;why='descriptive_homogeneous_border_estimator_nonuniformity_not_resolved'if ok else'nonisometric_window'if not isometric else'insufficient_type_or_interior_points'
     rows.append(dict(section=s,component=component,label_definition=label,high_definition=h,distance_um=r,n_A=len(aa),n_B=len(bb),interior_A=na,observed_cross_pairs=pairs,observation_area_um2=area,K=K,L=L,L_minus_r=L-r,status=why,border_pixel_resolution_um=8,window='supplied in-tissue observation union; not total histologic area',source_subset='globally remove Epi sharing any eligible T32 source',P=np.nan))
 print('K_BORDER_AUDIT',s,flush=True)
put(pd.DataFrame(rows),'standard_cross_KL.tsv');put(pd.DataFrame(audit),'standard_cross_KL_window_audit.tsv')
