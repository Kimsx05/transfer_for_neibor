import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','TF_NUM_INTRAOP_THREADS','TF_NUM_INTEROP_THREADS']:os.environ[k]='2' if k=='TF_NUM_INTRAOP_THREADS' else '1'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:6])
import sys,json
from pathlib import Path
import pandas as pd,numpy as np
from scipy import sparse,ndimage
from PIL import Image,ImageDraw
import tensorflow as tf
from stardist.models import StarDist2D
R=Path(__file__).resolve().parents[1];A=R.parent;N=A/'segmentation_restart_v1'
tf.config.set_logical_device_configuration(tf.config.list_physical_devices('GPU')[0],[tf.config.LogicalDeviceConfiguration(memory_limit=6000)])
w=A/'cell_annotation_freeze_v1/models/keras_cache/models/StarDist2D/2D_versatile_he/2D_versatile_he_extracted'
model=StarDist2D(None,name=w.name,basedir=w.parent)
def edges(mask):
 x=np.zeros(mask.shape,bool);x[:,1:]|=mask[:,1:]!=mask[:,:-1];x[1:]|=mask[1:]!=mask[:-1];return x
for s in ['Cervilla','D1_Johnson']:
 o=R/s;t=json.loads((N/s/'metadata/transform.json').read_text());im=np.load(N/s/'inputs/model_image.npy',mmap_mode='r')
 low=np.array(t['normalization']['low_rgb'],np.float32);high=np.array(t['normalization']['high_rgb'],np.float32)
 def norm(img):return (img.astype(np.float32)-low)/(high-low)
 pars=dict(prob_thresh=.25,nms_thresh=.40)
 if sys.argv[1]=='cal':
  old=sparse.load_npz(A/'segmentation_assignment_v2'/s/'inputs/nuclei.npz').tocsr()
  sheet=Image.new('RGB',(4*288,3*172),'white');draw=ImageDraw.Draw(sheet)
  for i,z in enumerate(pd.read_csv(R/'tables/calibration_fields.tsv',sep='\t').query('section==@s').itertuples()):
   y=int(np.clip(z.my,80,im.shape[0]-80));x=int(np.clip(z.mx,80,im.shape[1]-80));sl=np.s_[y-80:y+80,x-80:x+80];img=im[sl];lab=model.predict_instances(norm(img),**pars)[0]
   for j,mask in enumerate([old[sl].toarray(),lab]):
    v=img.copy();v[edges(mask)]=[10,200,195] if j==0 else [245,170,15];sheet.paste(Image.fromarray(v).resize((144,144)),((i%4)*288+j*144,(i//4)*172))
   draw.text(((i%4)*288,(i//4)*172+145),z.field+'   old | new',fill='black')
  sheet.save(R/'figures'/(s+'_calibration.jpg'),quality=94)
 else:
  path=o/'nuclei_v6.npz'
  if path.exists():continue
  lab=model.predict_instances_big(norm(im),axes='YXC',block_size=1536,min_overlap=160,context=128,**pars)[0]
  sparse.save_npz(path,sparse.csr_matrix(lab))
  print('FINAL_MASK',s,int(lab.max()),flush=True)
