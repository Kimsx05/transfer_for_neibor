from common5 import *
from PIL import Image,ImageDraw,ImageFont
from scipy.spatial import cKDTree
from matplotlib import colormaps
import subprocess
FONT=subprocess.check_output(['fc-match','-f','%{file}','Arial'],text=True)
def label(draw,xy,t,size=15):
 draw.text(xy,str(t),fill='black',font=ImageFont.truetype(FONT,size))
class View:
 def __init__(self,s):
  self.s=s;self.out=R/s;self.tr=js(N/s/'metadata/transform.json');self.au=js(N/s/'metadata/input_audit.json');self.mpp=self.tr['native_mpp'];self.crop=np.array(self.tr['fullres_crop_xyxy'][:2]);self.scale=np.array(self.tr['local_to_model_scale_xy'])
  self.img=np.load(N/s/'inputs/model_image.npy',mmap_mode='r');self.nu=sparse.load_npz(S/s/'inputs/nuclei.npz').toarray()
  self.grid=sparse.load_npz(self.out/'inferred_cell_grid.npz').toarray();self.off=np.load(self.out/'grid_origin_col_row.npy')
  self.aff=np.array(self.au['array_col_row_affine_fullres']);self.ia=np.linalg.inv(self.aff[:2]);self.D=np.load(self.out/'D_geometry.npz')
  self.d=pd.read_parquet(self.out/'cells_pre_MP10.parquet');self.c=pd.read_parquet(self.out/'candidate_ledger.parquet')
  self.raw=pd.read_parquet(self.out/'raw_candidates.parquet')
  self.t=np.zeros(self.grid.shape,np.uint8);self.g=self.t.copy()
  ri=self.raw.array_row.to_numpy()-self.off[1];ci=self.raw.array_col.to_numpy()-self.off[0];self.t[ri,ci]=self.raw[T].gt(0).any(axis=1);self.g[ri,ci]=self.raw[reg].gt(0).any(axis=1).astype(int)+self.raw.FOXP3.gt(0).astype(int)
 def sample(self,x,y,extent,n,mode):
  yy,xx=np.mgrid[:n,:n];xy=np.stack((x+(xx-(n-1)/2)*extent/n,y+(yy-(n-1)/2)*extent/n),axis=-1)
  orig=xy/self.mpp;pix=np.rint((orig-self.crop+.5)*self.scale-.5).astype(int);px,py=pix[...,0],pix[...,1];ok=(px>=0)&(py>=0)&(px<self.img.shape[1])&(py<self.img.shape[0]);rgb=np.full((n,n,3),242,np.uint8);rgb[ok]=self.img[py[ok],px[ok],:3];nu=np.zeros((n,n),np.int32);nu[ok]=self.nu[py[ok],px[ok]]
  if mode=='he':return rgb
  rc=np.rint((orig-self.aff[2])@self.ia).astype(int)-self.off;col,row=rc[...,0],rc[...,1];ok=(col>=0)&(row>=0)&(col<self.grid.shape[1])&(row<self.grid.shape[0]);col=np.clip(col,0,self.grid.shape[1]-1);row=np.clip(row,0,self.grid.shape[0]-1);cell=self.grid[row,col];cell[~ok]=0
  def edge(a):
   e=np.zeros(a.shape,bool);e[1:]|=a[1:]!=a[:-1];e[:,1:]|=a[:,1:]!=a[:,:-1];return e
  rgb[edge(cell)]=[230,160,40];rgb[edge(nu)]=[0,185,176];t=(self.t[row,col]>0)&ok;g=self.g[row,col]*ok
  rgb[t]=[30,100,245];rgb[g==1]=[173,51,176];rgb[g==2]=[251,44,35];rgb[t&(g>0)]=[248,245,18]
  return rgb
 def overview(self,r,second=False):
  ca=self.c[self.c.roi.eq(r)];dd=self.d[self.d.roi.eq(r)&self.d.in_D];lo=ca[['x','y']].min().to_numpy()-205;hi=ca[['x','y']].max().to_numpy()+205;cen=(lo+hi)/2;e=max(hi-lo)
  counts_rank=(dd.FOXP3.gt(0).astype(int)*4+dd.core_width+dd.T_width)
  ordidx=np.lexsort((dd.id,-dd.pTreg.fillna(-1),-counts_rank)) if len(dd)else []
  z=dd.iloc[ordidx[0]] if len(dd)else None;f=np.array([z.x,z.y])if z is not None else cen;num=int(z.id)if z is not None else -1
  canvas=Image.new('RGB',(1050,830),'white');draw=ImageDraw.Draw(canvas);title=f'{self.s}   ROI {r}   {len(ca):,} screened locations / objects'
  label(draw,(8,5),title,19)
  for j,(ct,extent,n,which,txt) in enumerate([(cen,e,340,'he','Registered H&E'),(cen,e,340,'he','Candidate core'),(f,90,340,'edge','Focused nuclear / RNA sources'),(f,420,340,'he','50 / 200 µm context'),(f,96,340,'edge','Cell class / source overlay'),(f,48,340,'edge','Detail and cell ID')]):
   img=Image.fromarray(self.sample(*ct,extent,n,which));dr=ImageDraw.Draw(img)
   def pos(x,y):return (int((x-ct[0])*n/extent+n/2),int((y-ct[1])*n/extent+n/2))
   if j==1:
    a=np.asarray(img).copy();p=np.rint((ca[['x','y']].to_numpy()-ct)*n/extent+n/2).astype(int).clip(0,n-1);a[p[:,1],p[:,0]]=[20,144,110];img=Image.fromarray(a);dr=ImageDraw.Draw(img)
   if j in [3,4]:
    near=dd[(dd.x-ct[0]).abs().lt(extent/2)&(dd.y-ct[1]).abs().lt(extent/2)]
    for c in near.itertuples():
     px,py=pos(c.x,c.y)
     color=None
     if second and getattr(c,'score_ok',False):color=tuple((colormaps['viridis'](min(1,c.MP10_cell_score/self.vmax))[:3]*np.array([255]*3)).astype(int))
     elif not second and c.epithelial:color=(244,175,86)
     elif not second and c.status=='supported_other_T':color=(50,179,149)
     elif not second and (c.pred=='Treg' or c.FOXP3>0):color=(167,61,170)
     if color:dr.ellipse([px-2,py-2,px+2,py+2],fill=color)
    if j==3:
     for rr in [50,200]:
      w=rr*n/extent;dr.ellipse([n/2-w,n/2-w,n/2+w,n/2+w],outline='white',width=2)
   if j==5 and z is not None:dr.text((4,4),str(num),fill='black',stroke_width=2,stroke_fill='white')
   loc=((j%3)*350,43+(j//3)*377);canvas.paste(img,loc);label(draw,(loc[0]+2,loc[1]+342),txt if not(second and j in [3,4])else 'Cell MP10 (scoreable epithelium)',14)
  label(draw,(5,801),'Teal nucleus | gold inferred extent | blue T RNA | purple program | red FOXP3 | yellow co-detection',14)
  if second:
   for i in range(125):
    color=tuple(int(c*255) for c in colormaps['viridis'](i/124)[:3]);draw.line((656+i,807,656+i,819),fill=color)
   label(draw,(785,801),f'0 to {self.vmax:.3f} MP10',13)
  path=R/'figures'/f"{self.s}_ROI{r}_{'MP10'if second else 'geometry'}.jpg";canvas.save(path,quality=91);return dict(section=self.s,roi=r,focus_id=num,x=f[0],y=f[1],path=str(path),selection='max FOXP3/core/T width then reference Treg support; independent of MP10')
 def detail(self):
  old=pd.read_parquet(self.out/'old_model_candidates.parquet');oldids=set(old[old.pred.eq('Treg')].id)
  d=self.d[self.d.in_D & (self.d.pred.eq('Treg')|self.d.reg_UMI.gt(0)|self.d.mixed|self.d.id.isin(oldids))][['id','x','y','roi','T_UMI','reg_UMI','FOXP3','Epi_UMI','Fib_UMI','status']].copy()
  raw=self.raw.copy()
  o=pd.read_parquet(self.out/'assignments.parquet',columns=['new_id']).new_id.to_numpy();raw=raw[(o[raw.source_bin_index]==0)&raw[reg].gt(0).any(axis=1)].copy();raw['unit']='2:'+raw.source_bin_index.astype(str);raw=raw.merge(self.c[['unit','roi']],on='unit');raw['id']=-(raw.source_bin_index+1);raw['T_UMI']=raw[T].sum(axis=1);raw['reg_UMI']=raw[reg].sum(axis=1);raw['status']='unassigned'
  d=pd.concat([d,raw[d.columns.intersection(raw.columns)]],ignore_index=True);d=d.sort_values(['roi','y','x']).reset_index(drop=True);entries=[]
  for r,sub in d.groupby('roi'):
   for pg,start in enumerate(range(0,len(sub),400)):
    part=sub.iloc[start:start+400];nx=20;h=67;w=54;canvas=Image.new('RGB',(nx*w,((len(part)+nx-1)//nx)*h+24),'white');draw=ImageDraw.Draw(canvas);label(draw,(3,2),f'{self.s} | ROI {r} | detail {pg+1} | 26 µm per field',14)
    for j,q in enumerate(part.itertuples()):
     xx=j%nx*w;yy=j//nx*h+24;arr=self.sample(q.x,q.y,26,52,'edge');canvas.paste(Image.fromarray(arr),(xx,yy));draw.text((xx,yy+52),str(q.id),fill='black')
    path=R/'figures'/'details'/f'{self.s}_{r:02d}_{pg+1:03d}.jpg';canvas.save(path,quality=88)
    part=part.copy();part['file']=str(path);part['slot']=np.arange(len(part));entries.append(part)
  pd.concat(entries).to_parquet(self.out/'zoom_index.parquet',index=False)
if __name__=='__main__':
 second='second'in sys.argv; (R/'figures/details').mkdir(exist_ok=True);rows=[]
 for s in sections:
  v=View(s)
  if second:v.d=pd.read_parquet(v.out/'cells_scored.parquet');v.vmax=v.d.MP10_cell_score.max()
  for r in sorted(v.c.roi.unique()):rows.append(v.overview(r,second))
  if 'details'in sys.argv:v.detail()
  del v
 put(pd.DataFrame(rows),R/'tables'/('round2_images.tsv'if second else 'round1_images.tsv'))
