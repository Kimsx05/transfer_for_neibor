from base_data import *
# Compare registered RGB intensity at assigned 2 µm bin centres to section-specific true H&E nuclear bins.
for s in HD:
 d=paths(s);b=pd.read_parquet(H/s/'inputs/bin_to_cell.parquet',columns=['array_row','array_col','nucleus_label','image_covered'])
 members=pd.read_parquet(d/'members.parquet',columns=['final_label','cell_label']);sel=(members.final_label>0)&(members.cell_label==0)
 t=json.loads((H/s/'metadata/transform.json').read_text());tx=np.array(json.loads((d/'conservation.json').read_text())['grid_to_fullres'])
 im=np.load(H/s/'inputs/model_image.npy',mmap_mode='r')
 def lum(tab):
  xy=np.c_[tab.array_row,tab.array_col,np.ones(len(tab))]@tx
  mp=np.rint((xy-np.array(t['fullres_crop_xyxy'][:2])+.5)*t['local_to_model_scale_xy']-.5).astype(int)
  return np.mean(im[np.clip(mp[:,1],0,im.shape[0]-1),np.clip(mp[:,0],0,im.shape[1]-1)],axis=1)
 ref=lum(b[(b.nucleus_label>0)&b.image_covered].iloc[::23]);threshold=np.quantile(ref,.90)
 table=pd.DataFrame({'label':members.final_label[sel],'nuclear_range':lum(b[sel])<=threshold,'intensity':lum(b[sel])})
 v=table.groupby('label').agg(fraction_nuclear_range=('nuclear_range','mean'),median_intensity=('intensity','median'))
 v['section_nuclear_P90']=threshold;v['image_suspect']=v.fraction_nuclear_range<.25
 v.index=[f'{s}__retry__{i}' for i in v.index];v.to_parquet(d/'secondary_image_support.parquet')
 print(s,len(ref),threshold,len(v),v.image_suspect.sum())

