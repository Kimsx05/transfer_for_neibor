from base6 import *
from PIL import Image,ImageDraw,ImageFont
Image.MAX_IMAGE_PIXELS=None
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
def font(n):
 try:return ImageFont.truetype(FONT,n)
 except OSError:return ImageFont.load_default()
class Cards:
 def __init__(self,s):
  self.s=s;o=R/s;self.o=o;tr=js(N/s/'metadata/transform.json');self.m=tr['native_mpp'];self.org=np.array(tr['fullres_crop_xyxy'][:2]);self.sc=np.array(tr['local_to_model_scale_xy']);self.aff=np.array(js(N/s/'metadata/input_audit.json')['array_col_row_affine_fullres']);self.he=Image.open(tr['source_image']);self.he.load()
  self.n=sparse.load_npz(o/'final_nuclei.npz').tocsr();self.old=sparse.load_npz(S/s/'inputs/nuclei.npz').tocsr();self.panel=sparse.load_npz(o/'observed_panel.npz').tocsr()
  self.native_polys=None
  if R!=ROOT:
   import gzip,json
   from shapely.geometry import shape
   pp=ROOT/s/'native_full';off=pp/'native_offset_um.txt';self.shift=float(off.read_text()) if off.exists() else 0.
   m=pd.read_parquet(pp/'registered/cells.parquet').set_index('cell').original_cell_id
   feats=json.load(gzip.open(pp/'registered/cells.geojson.gz'))['features']
   self.native_polys={int(m.loc[z['properties']['cell']]):shape(z['geometry']) for z in feats}
  self.d=pd.read_parquet(o/'identity_frozen.parquet').set_index('id',drop=False);self.q=pd.read_parquet(o/'review_queue.parquet').set_index('review_id');a=pd.read_parquet(o/'assignment.parquet');self.hard=a.id.to_numpy();self.amb=a.ambiguous.to_numpy();self.lookup=np.full((a.array_row.max()+1,a.array_col.max()+1),-1,'int32');self.lookup[a.array_row,a.array_col]=a.source_bin_index
 def card(self,k):
  q=self.q.loc[k];s=self.s;xy=q[['x','y']].to_numpy(dtype=float);i=int(q.id) if pd.notna(q.id) else 0;z=self.d.loc[i] if i else None
  mid=xy/self.m-self.org;w=178;span=36;scale=w/span
  def cut(span):
   v=span/(2*self.m);return self.he.crop(tuple(np.rint(np.r_[mid-v,mid+v]).astype(int))).convert('RGB').resize((w,w),Image.Resampling.BILINEAR)
  he=cut(span);context=cut(140);dr=ImageDraw.Draw(context);off=(140-span)/280*w;dr.rectangle((off,off,w-off,w-off),outline='orange',width=2)
  ly,lx=np.indices((w,w));px=(np.stack([lx,ly],axis=2)+.5-w/2)/scale/self.m+mid;mp=np.rint((px+.5)*self.sc-.5).astype(int)
  def cropmask(mask):
   xx,yy=mp[:,:,0].clip(0,mask.shape[1]-1),mp[:,:,1].clip(0,mask.shape[0]-1);x0,x1=xx.min(),xx.max()+1;y0,y1=yy.min(),yy.max()+1
   v=mask[y0:y1,x0:x1].toarray()[yy-y0,xx-x0];v[(mp[:,:,0]!=xx)|(mp[:,:,1]!=yy)]=0;return v
  def overlay(mask):
   v=cropmask(mask);edge=(v!=np.roll(v,1,axis=0))|(v!=np.roll(v,1,axis=1));edge[0]=False;edge[:,0]=False;b=np.array(he);b[edge]=[0,183,201];return Image.fromarray(b)
  nimg=overlay(self.n);oldimg=overlay(self.old)
  if i:
   v=cropmask(self.n)==i;b=np.array(nimg);b[v&~ndimage.binary_erosion(v,iterations=2)]=[255,147,0];nimg=Image.fromarray(b)
  rc=np.rint((xy/self.m-self.aff[2])@np.linalg.inv(self.aff[:2])).astype(int)
  cols,rows=np.meshgrid(np.arange(rc[0]-15,rc[0]+16),np.arange(rc[1]-15,rc[1]+16));rows=rows.ravel();cols=cols.ravel();good=(rows>=0)&(cols>=0)&(rows<self.lookup.shape[0])&(cols<self.lookup.shape[1]);rows,cols=rows[good],cols[good];idx=self.lookup[rows,cols];valid=idx>=0;rows,cols,idx=rows[valid],cols[valid],idx[valid]
  points=((np.c_[cols,rows]@self.aff[:2]+self.aff[2])*self.m-xy)*scale+w/2
  v=self.panel[idx].toarray();assign=nimg.copy();d=ImageDraw.Draw(assign)
  for j,(u,y) in enumerate(points):
   if i and self.hard[idx[j]]==i:d.rectangle((u-3,y-3,u+3,y+3),outline='lime',width=2)
   elif self.amb[idx[j]]:d.line((u-3,y-3,u+3,y+3),fill='red',width=1)
  if self.native_polys is not None and i in self.native_polys:
   for poly in self.native_polys[i].geoms:
    line=(np.asarray(poly.exterior.coords)[:,::-1]-self.shift-xy)*scale+w/2
    d.line(list(map(tuple,line)),fill='#efef02',width=2)
  imgs=[context,oldimg,nimg,assign];labs=['H&E 140 um','Old nuclei 36 um','New nuclei / focus orange','Assigned grids / red ambiguity']
  for g,name,col in [(0,'CD3D/E/G or TRAC','#06a94c'),(1,'FOXP3','#f40060'),(2,'IL2RA/CTLA4/IKZF2','#c000f6'),(3,'Epithelial / collagen / myeloid',None)]:
   p=nimg.copy();dr=ImageDraw.Draw(p)
   for j,(u,y) in enumerate(points):
    for gg in ([g] if col else [3,4,5]):
     if v[j,gg]>0:
      c=col or ['#db8210','#005eff','#8d10df'][gg-3];r=2.3+min(3,np.sqrt(v[j,gg]));dr.ellipse((u-r,y-r,u+r,y+r),fill=c,outline='black',width=1)
   imgs.append(p);labs.append(name)
  out=Image.new('RGB',(4*w,2*(w+18)+86),'white');dr=ImageDraw.Draw(out)
  for j,(p,n) in enumerate(zip(imgs,labs)):
   x=(j%4)*w;y=(j//4)*(w+18)+20;out.paste(p,(x,y));dr.text((x+2,y+w+1),n,font=font(10),fill='black')
  dr.text((2,1),f'{s}   {k}   ROI {int(q.roi)}   x,y {xy[0]:.1f}, {xy[1]:.1f} um   {q.outcome}',fill='black',font=font(13))
  if i:
   detail=f'old {z.v5_old_id:g} {z.v5_pred} ({z.v5_pTreg:.2f}) | new {z.major}; {z.T_subtype} | pT {z.p_T:.2f} pReg|T {z.p_Treg_given_T:.2f} | {z.mixed_v6}'
   d2=f'UMI {z.nCount}  T {int(z[direct].sum())}  FOXP3 {z.FOXP3}  IL2RA {z.IL2RA}  CTLA4 {z.CTLA4}  IKZF2 {z.IKZF2}  | {z.nuclear_action}'
  else:detail='Source has no final cell assignment';d2='Signal must remain unlocalized.'
  import textwrap
  dr.multiline_text((3,20+2*(w+18)),'\n'.join(textwrap.wrap(detail,100)),fill='black',font=font(12));dr.text((3,62+2*(w+18)),d2,fill='black',font=font(12));return out
 def save(self,k):
  p=self.o/'cards'/(k.replace(':','_')+'.jpg');p.parent.mkdir(exist_ok=True);self.card(k).save(p,quality=90);return p
if __name__=='__main__':
 for s in sections if sys.argv[1]=='both' else [sys.argv[1]]:
  c=Cards(s);d=c.d
  if sys.argv[2]=='index':
   ids=d.index[d.extended|d.control].tolist();tr=rd(c.o/'priority_trace.tsv');ids+=tr.new_id.dropna().astype(int).tolist()
   keys=['cell:'+str(i) for i in sorted(set(ids))];keys+=list(c.q.groupby('roi').head(1).index);keys+=list(c.q[c.q.id.isna()].groupby('roi').head(1).index)
  else:keys=list(c.q.index) if sys.argv[2]=='all' else sys.argv[2:]
  for k in dict.fromkeys(keys):
   if k in c.q.index:c.save(k)
