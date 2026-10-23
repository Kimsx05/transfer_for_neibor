from common import *
setup_gpu_env()
import numpy as np,pandas as pd,sys,random,time
import tensorflow as tf
from stardist.models import StarDist2D
from scipy import sparse
from PIL import Image
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
Image.MAX_IMAGE_PIXELS=None
np.random.seed(SEED);random.seed(SEED);tf.random.set_seed(SEED)
gpus=tf.config.list_physical_devices('GPU');assert gpus,'GPU runtime required before full segmentation'
tf.config.set_logical_device_configuration(gpus[0],[tf.config.LogicalDeviceConfiguration(memory_limit=6144)])
weights=V1/'models/keras_cache/models/StarDist2D/2D_versatile_he/2D_versatile_he_extracted'
model=StarDist2D(None,name=weights.name,basedir=weights.parent)
print('GPU',gpus,'StarDist loaded cached weights',flush=True)
name=sys.argv[1];phase=sys.argv[2];meta,base,out=section(name);t=json.loads((out/'metadata/transform.json').read_text());im=np.load(out/'inputs/model_image.npy',mmap_mode='r');norminfo=t['normalization'];lo=np.array(norminfo['low_rgb'],dtype='float32');hi=np.array(norminfo['high_rgb'],dtype='float32')
def norm(x):return (x.astype('float32')-lo)/(hi-lo+1e-20)
if phase=='pilot':
 fig,axs=plt.subplots(3,3,figsize=(13,12));rows=[]
 pos=pd.read_parquet(base/'spatial/tissue_positions.parquet',columns=['in_tissue','pxl_col_in_fullres','pxl_row_in_fullres']);pos=pos[pos.in_tissue==1];xy=((pos[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()-t['fullres_crop_xyxy'][:2]+.5)*t['local_to_model_scale_xy']-.5)
 for i,(label,box) in enumerate(t['representative_model_bboxes'].items()):
  x0,y0,x1,y1=box;raw=np.array(im[y0:y1,x0:x1]);labels,details=model.predict_instances(norm(raw),prob_thresh=.05,nms_thresh=.3,n_tiles=(2,2,1));overlay=raw.copy();overlay[find_boundaries(labels,mode='inner')]=[20,230,60]
  axs[i,0].imshow(raw);axs[i,0].set_title(f'{label}: source {t["native_mpp"]:.3f} um/px, model ~0.5')
  axs[i,1].imshow(overlay);axs[i,1].set_title(f'{labels.max()} nuclei; prob=.05,nms=.3')
  hh,ww=raw.shape[:2];axs[i,2].imshow(overlay);hit=(xy[:,0]>=x0)&(xy[:,0]<x1)&(xy[:,1]>=y0)&(xy[:,1]<y1);local=xy[hit]-[x0,y0];axs[i,2].scatter(local[:,0],local[:,1],s=1,c='#00bfff',alpha=.35);axs[i,2].set_xlim(ww*.35,ww*.65);axs[i,2].set_ylim(hh*.65,hh*.35);axs[i,2].set_title('Center detail + observed in-tissue bin centers')
  for ax in axs[i]:ax.axis('off')
  areas=np.bincount(labels.ravel())[1:]*np.prod(t['actual_model_mpp_xy']);rows.append(dict(roi=label,model_bbox=box,nuclei=int(labels.max()),bins=int(hit.sum()),nuclear_area_um2_quantiles=np.quantile(areas,[0,.1,.5,.9,.99,1]).tolist()))
  sparse.save_npz(out/f'inputs/pilot_{label}_nuclei.npz',sparse.csr_matrix(labels));print('PILOT',name,label,rows[-1],flush=True)
 fig.tight_layout();fig.savefig(out/'figures/segmentation_pilot.png',dpi=160,transparent=True);plt.close(fig);jwrite(out/'metadata/segmentation_pilot.json',rows)
else:
 assert (out/'metadata/pilot_review.json').exists(),'Review pilot before full inference'
 start=time.monotonic();print('FULL_START',name,im.shape,flush=True);normalized=norm(im)
 labels,details=model.predict_instances_big(normalized,axes='YXC',block_size=2048,min_overlap=128,context=128,prob_thresh=.05,nms_thresh=.3)
 del normalized;sparse.save_npz(out/'inputs/nuclei.npz',sparse.csr_matrix(labels),compressed=True);np.savez_compressed(out/'inputs/nucleus_geometry.npz',points=details['points'],prob=details['prob'])
 summary=dict(section=name,n_nuclei=int(labels.max()),shape=list(labels.shape),model='2D_versatile_he',prob_thresh=.05,nms_thresh=.3,block_size=2048,min_overlap=128,context=128,device='GPU0 RTX6000; TensorFlow memory cap 6144MiB',native_mpp=t['native_mpp'],actual_model_mpp_xy=t['actual_model_mpp_xy'],runtime_seconds=time.monotonic()-start,seed=123)
 jwrite(out/'metadata/segmentation_summary.json',summary);print('SEGMENTATION_DONE',summary,flush=True)
