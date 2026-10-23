from common import *
import numpy as np,tifffile,zarr
from PIL import Image
from scipy.ndimage import map_coordinates
Image.MAX_IMAGE_PIXELS=None
for name in ['Cervilla','D1_Johnson']:
 m,base,out=section(name);w,h=m['exported_dimensions_xy'];mpp=m['microns_per_pixel'];nw,nh=round(w*mpp/.5),round(h*mpp/.5);sx,sy=nw/w,nh/h;path=out/'inputs/model_image.npy'
 t=dict(section=name,source_image=m['image_path'],source_image_shape_yx=[h,w],model_image_shape_yx=[nh,nw],native_mpp=mpp,model_nominal_mpp=.5,actual_model_mpp_xy=[mpp/sx,mpp/sy],local_to_model_scale_xy=[sx,sy],fullres_crop_xyxy=m['source_fullres_crop_xyxy'],model_xy_formula='((fullres_xy-crop_origin_xy)+0.5)*local_to_model_scale_xy-0.5',physical_xy_formula='original_fullres_xy*native_mpp; no batch correction',original_image_recovered=True,rotation_added=0,flip_added=False,seed=123)
 if not path.exists():
  print('BUILD_VIEW',name,[nh,nw],flush=True);dest=np.lib.format.open_memmap(path,mode='w+',dtype=np.uint8,shape=(nh,nw,3))
  if name=='Cervilla':
   im=Image.open(m['image_path']);im.load();dest[:]=np.asarray(im.resize((nw,nh),Image.Resampling.BILINEAR));del im
  else:
   store=tifffile.imread(m['image_path'],aszarr=True);src=zarr.open(store,mode='r');xs=(np.arange(nw,dtype=np.float32)+.5)/sx-.5
   for y0 in range(0,nh,256):
    y1=min(y0+256,nh);ys=(np.arange(y0,y1,dtype=np.float32)+.5)/sy-.5;lo=max(0,int(np.floor(ys.min()))-1);hi=min(h,int(np.ceil(ys.max()))+2);block=np.asarray(src[lo:hi,:,:]);yy,xx=np.meshgrid(ys-lo,xs,indexing='ij');coords=np.stack([yy,xx])
    for c in range(3):dest[y0:y1,:,c]=map_coordinates(block[:,:,c],coords,order=1,mode='nearest',prefilter=False)
    if y0%2048==0:print(name,'ROWS',y1,'/',nh,flush=True)
   store.close()
  dest.flush();del dest
 im=np.load(path,mmap_mode='r');sample=np.asarray(im[::8,::8]);lo=np.percentile(sample,1,axis=(0,1));hi=np.percentile(sample,99.8,axis=(0,1));t['normalization']={'percentiles':[1,99.8],'stride_sample':8,'low_rgb':lo.tolist(),'high_rgb':hi.tolist(),'scope':'one fixed per-section channel normalization, shared across all tiles'}
 preview=Image.fromarray(np.asarray(im[::8,::8]));preview.thumbnail((1500,1500));preview.save(out/'inputs/analysis_preview.jpg',quality=90)
 # Original v1 Cervilla ROI footprints; D1 chosen from whole image morphology before labels.
 if name=='Cervilla':
  old=np.array([[1500,2100,1800,2400],[2450,2700,2750,3000],[3600,3380,3900,3680]])/m['fullres_to_hires_scale'];boxes=np.rint(old*np.array([sx,sy,sx,sy])).astype(int).tolist()
 else:
  boxes=[[round(x*nw),round(y*nh),round(x*nw)+768,round(y*nh)+768] for x,y in [(.62,.26),(.43,.48),(.12,.53)]]
 t['representative_model_bboxes']=dict(zip(['dense','stroma','edge'],boxes));jwrite(out/'metadata/transform.json',t);print('VIEW_READY',name,flush=True)
