from base6 import *
from PIL import Image,ImageDraw
Image.MAX_IMAGE_PIXELS=None
rows=[]
for s in sections:
 o=ROOT/s;d=pd.read_parquet(o/'cells_geometry.parquet').set_index('id');ov=pd.read_parquet(o/'nuclear_overlap.parquet');tr=js(N/s/'metadata/transform.json');m=tr['native_mpp'];sc=np.array(tr['local_to_model_scale_xy']);org=np.array(tr['fullres_crop_xyxy'][:2])
 new=sparse.load_npz(o/'final_nuclei.npz').tocsr();old=sparse.load_npz(S/s/'inputs/nuclei.npz').tocsr();he=Image.open(tr['source_image']);he.load()
 near=ov[(ov.overlap_old>.15)&(ov.overlap_new>.45)&(ov.intersection>=20)]
 split=near[near.old_id.isin(near.old_id.value_counts().loc[lambda x:x>=2].index)].sort_values(['old_id','intersection'],ascending=[True,False])
 nm=ov[(ov.overlap_old>.65)&(ov.overlap_new>.20)];merge=nm[nm.new_id.isin(nm.new_id.value_counts().loc[lambda x:x>=2].index)].sort_values('new_id')
 choose=[]
 for k,ids in [('split',split.new_id),('merged',merge.new_id),('new',d.index[d.nuclear_action.eq('new_detection')]),('rescued',d.index[d.nuclear_action.eq('retain_old_visible_unmatched')])]:
  ids=pd.Index(ids).unique();z=d.loc[ids];z=z[(z.nCount>=25)&(z.nuclear_area_um2>=10)&(z.nuclear_area_um2<150)]
  if len(z):choose.append((k,z.sort_values('nuclear_area_um2').iloc[len(z)//2]))
 out=Image.new('RGB',(640,204*len(choose)),'white');dr=ImageDraw.Draw(out)
 for r,(k,z) in enumerate(choose):
  i=z.name;xy=z[['x','y']].to_numpy(float);mid=xy/m-org;sz=192;span=42;dy,dx=np.indices((sz,sz));pix=(np.stack([dx,dy],axis=2)+.5-sz/2)*span/(sz*m)+mid;mp=np.rint((pix+.5)*sc-.5).astype(int);x,y=mp[:,:,0],mp[:,:,1]
  a=np.array(he.crop(tuple(np.rint(np.r_[mid-span/(2*m),mid+span/(2*m)]).astype(int))).convert('RGB').resize((sz,sz)));out.paste(Image.fromarray(a),(0,r*204+12))
  ids=[]
  for col,ma in enumerate([old,new]):
   x0,x1=x.min(),x.max()+1;y0,y1=y.min(),y.max()+1;v=ma[y0:y1,x0:x1].toarray()[y-y0,x-x0];edge=(v!=np.roll(v,1,axis=0))|(v!=np.roll(v,1,axis=1));edge[0]=False;edge[:,0]=False
   b=a.copy();b[edge]=[0,200,220] if col==0 else [255,160,0];out.paste(Image.fromarray(b),((col+1)*sz,r*204+12))
  prev=ov[ov.new_id.eq(i)&ov.substantial].old_id.tolist();dr.text((2,r*204),f'{k}: new {i}, overlapping old {prev}   |   H&E   |   previous   |   final',fill='black')
  rows.append((s,k,i,','.join(map(str,prev)),z.x,z.y))
 out.save(ROOT/'figures'/(s+'_nuclear_changes.jpg'),quality=94)
tab(pd.DataFrame(rows,columns=['section','action','id','overlap_old_ids','x','y']),'nuclear_examples.tsv')
