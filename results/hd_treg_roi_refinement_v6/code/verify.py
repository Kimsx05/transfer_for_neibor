from base6 import *
from scipy.signal import fftconvolve
from load_result import load_section
rows=[]
for s in sections:
 o=R/s;branch='conservative' if R==ROOT else 'proseg_observed';res=load_section(s,branch,True);d=res['cells'];X=res['counts'];src=res['sources'];review=res['review'];nb=res['neighborhoods'];m=res['matches']
 assert np.array_equal(d.index,np.arange(1,len(d)+1)) and d[['x','y']].notna().all().all()
 b=rd(o/'count_balance.tsv.gz');assert len(b)==18085 and np.array_equal(b.total,b.assigned+b.excluded+b.ambiguous+b.unassigned) and np.array_equal(np.asarray(X.sum(0)).ravel(),b.assigned) and np.array_equal(np.asarray(X.sum(1)).ravel(),d.nCount)
 assert (d['main']==(d.A|d.B)).all() and (~d['main']|d.extended).all() and (~d.strict|d.A).all() and not (d['main']&d.control).any()
 assert (d.loc[d['main'],'nCount']>=25).all() and (d.loc[d['main'],'direct_width']>=1).all()
 assert not src.unit.duplicated().any() and not review.review_id.duplicated().any() and not nb.duplicated(['id','radius']).any()
 assert d.loc[~d.score_ok,'MP10_cell_score'].isna().all() and nb.loc[nb.scored==0,'mean_MP10'].isna().all()
 if len(m):assert m.distance_um.le(1000).all() and m.groupby(['case','radius','control_pool']).size().le(3).all() and not m.duplicated(['case','control','radius','control_pool']).any()
 links=pd.read_parquet(o/'source_links.parquet');assert not links.duplicated().any() and links.id.isin(d.index).all()
 a=pd.read_parquet(o/'assignment.parquet',columns=['source_bin_index','array_row','array_col','usable','id']);assert np.array_equal(a.source_bin_index,np.arange(len(a))) and a.id.between(0,len(d)).all()
 au=js(N/s/'metadata/input_audit.json');aff=np.array(au['array_col_row_affine_fullres']);ij=np.rint((src[['x','y']].to_numpy()/au['native_mpp']-aff[2])@np.linalg.inv(aff[:2])).astype(int)
 xx=a.array_col.to_numpy();yy=a.array_row.to_numpy();origin=np.minimum([xx.min(),yy.min()],ij.min(0))-1;xx=xx-origin[0];yy=yy-origin[1];ij-=origin;shape=np.maximum([xx.max(),yy.max()],ij.max(0))+2;grid=np.zeros(shape[::-1],dtype=np.float32)
 cc=src[['unit','roi']].copy()
 for rad in [50,200]:
  width=rad//2;u,v=np.ogrid[-width:width+1,-width:width+1];kernel=(u*u+v*v<=width*width).astype(np.float32)
  for tag,mask in [('supplied',1),('usable',a.usable.to_numpy())]:
   grid[yy,xx]=mask;cc[f'{tag}_{rad}']=np.rint(fftconvolve(grid,kernel,mode='same')[ij[:,1],ij[:,0]]).astype(int)
  cc[f'theoretical_{rad}']=int(kernel.sum())
 cc.to_parquet(o/'candidate_grid_coverage.parquet',index=False)
 geom=rd(o/'priority_final.tsv');raw=pd.read_parquet(o/'raw_signal.parquet')
 row=dict(section=s,source_entries=len(src),unlocalized_sources=int((src.n_final_cells==0).sum()),review_units=len(review),reviewed=int(review.visual_note.notna().sum()),pending=int(review.visual_note.isna().sum()),roi=src.roi.nunique(),nuclei=len(d),in_D=d.in_D.sum(),nuclear_bins_zero=d.nuclear_bins.eq(0).sum(),supplied_UMI=int(b.total.sum()),excluded=int(b.excluded.sum()),unassigned=int(b.unassigned.sum()),ambiguous=int(b.ambiguous.sum()),assigned=int(b.assigned.sum()),input_bins=len(a),raw_signal_bins=len(raw),raw_unassigned=(raw.id==0).sum())
 row.update({k:int(d[k].sum()) for k in ['A','B','main','strict','extended','control','epi','score_ok']});rows.append(row)
tab(pd.DataFrame(rows),'completion.tsv')
