from base_data import *
from skimage.measure import label
from scipy import ndimage
sections=[]
for s in HD:
 d=paths(s);a=dataset(s);b=pd.read_parquet(d/'members.parquet',columns=['array_row','array_col','final_label','cell_label','image_covered']);r,c=b.array_row.values,b.array_col.values;g=np.zeros((int(r.max())+1,int(c.max())+1),'i4');g[r,c]=b.final_label
 comp=label(g,background=0,connectivity=2);valid=g.ravel()>0;num=int(comp.max());member_owner=np.zeros(num+1,'i4');member_owner[comp.ravel()[valid]]=g.ravel()[valid];u=a.obs.label.to_numpy();n_parts=np.bincount(member_owner[1:],minlength=int(u.max())+1)[u]
 masks=pd.DataFrame({'components_8_connected':n_parts,'observed_bin_fragmentation_flag':n_parts>1,'source':a.obs.source.to_numpy(),'area_um2':a.obs.area_um2,'initial_or_image_suspect':a.obs.geometry_flag},index=a.obs_names)
 masks.to_parquet(d/'geometry_qc.parquet');assert not ((b.final_label>0)&~b.image_covered).any()
 for k,ids in masks.groupby('source').indices.items():
  z=masks.iloc[ids];sections.append((s,k,len(z),int(z.observed_bin_fragmentation_flag.sum()),np.median(z.area_um2),int(z.initial_or_image_suspect.sum()),0))
write(pd.DataFrame(sections,columns=['sample','source','objects','multi_component_observed_2um_footprint','median_area_um2','initial_or_image_suspect','assigned_outside_image']), 'hd_geometry_QC.tsv')

