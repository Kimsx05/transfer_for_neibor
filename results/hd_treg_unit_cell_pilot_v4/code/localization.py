from base import *
from scipy.spatial import cKDTree
rows=[]
for s in sections:
 d=pd.read_parquet(R/s/'sources.parquet');rois=rd(R/'roi_manifest.tsv').query('section==@s');d['roi']=roi_for(d.x_um,d.y_um,rois);d['T_UMI']=d[T].sum(axis=1)
 target=d[d.input_usable&d.T_UMI.gt(0)].reset_index(drop=True);tree=cKDTree(target[['x_um','y_um']])
 for q in d[d.input_usable&d.FOXP3.gt(0)&d.roi.ne('')].itertuples():
  near=tree.query_ball_point([q.x_um,q.y_um],12);n=target.iloc[near];dist=np.hypot(n.x_um-q.x_um,n.y_um-q.y_um)
  a=dict(section=s,roi=q.roi,barcode=q.barcode,x_um=q.x_um,y_um=q.y_um,FOXP3=q.FOXP3,T_same_capture_UMI=q.T_UMI,local_T_capture_bins=len(near),nearest_T_um=dist.min()if len(n)else np.nan)
  for key in ['nucleus_label','v1_cell_label','cell_label']:
   label=getattr(q,key);same=n[key].eq(label)&(label>0);other=n[key].gt(0)&~same
   a[key]=label;a[key+'_T_same_nucleus_within12']=n.loc[same,'T_UMI'].sum();a[key+'_T_other_nuclei_within12']=n.loc[other,key].nunique()
  rows.append(a)
put(pd.DataFrame(rows),'FOXP3_capture_localization.tsv')
