# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Seeded H&E-only feasibility pilot; masks are NOT accepted cells or frozen labels."""
import os
from pathlib import Path
root=Path(__file__).resolve().parents[1];hd=root.parents[1]
os.environ['KERAS_HOME']=str(root/'models/keras_cache')
os.environ['CUDA_VISIBLE_DEVICES']='-1';os.environ['TF_NUM_INTRAOP_THREADS']='8';os.environ['TF_NUM_INTEROP_THREADS']='2';os.environ['OMP_NUM_THREADS']='8'
import numpy as np, pandas as pd,json,random
import tensorflow as tf
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from scipy import sparse
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops_table
random.seed(123);np.random.seed(123);tf.random.set_seed(123)
(root/'models').mkdir(exist_ok=True);(root/'inputs/pilot_masks').mkdir(exist_ok=True)
model=StarDist2D.from_pretrained('2D_versatile_he')
# 3 representative Cervilla tissue regions, selected before segmentation or labels.
rois=[('Cervilla','dense',1500,2100,1800,2400),('Cervilla','stroma',2450,2700,2750,3000),('Cervilla','edge',3600,3380,3900,3680)]
imgpath=hd/'data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um/spatial/tissue_hires_image.png'
sf=json.loads((imgpath.parent/'scalefactors_json.json').read_text());native_mpp=sf['microns_per_pixel']/sf['tissue_hires_scalef'];target_mpp=.5
im=Image.open(imgpath).convert('RGB');rows=[]
fig,axes=plt.subplots(len(rois),3,figsize=(13,13))
for i,(ds,name,x0,y0,x1,y1) in enumerate(rois):
 raw=np.asarray(im.crop((x0,y0,x1,y1)))
 size=(round((x1-x0)*native_mpp/target_mpp),round((y1-y0)*native_mpp/target_mpp))
 scaled=np.asarray(Image.fromarray(raw).resize(size,Image.Resampling.BILINEAR))
 labels,details=model.predict_instances(normalize(scaled,1,99.8,axis=(0,1)),prob_thresh=.05,nms_thresh=.3,n_tiles=(2,2,1))
 sparse.save_npz(root/f'inputs/pilot_masks/Cervilla_{name}.npz',sparse.csr_matrix(labels))
 props=pd.DataFrame(regionprops_table(labels,properties=['label','area','centroid','eccentricity']));props['area_um2']=props.area*target_mpp**2;props.to_csv(root/f'metadata/pilot_nuclei_{name}.tsv',sep='\t',index=False)
 overlay=scaled.copy();overlay[find_boundaries(labels,mode='inner')]=[255,40,40]
 axes[i,0].imshow(raw,interpolation='nearest');axes[i,0].set_title(f'{name}: native {native_mpp:.3f} um/px')
 axes[i,1].imshow(overlay);axes[i,1].set_title(f'StarDist pilot: {labels.max()} nuclei')
 # Fixed central close-up, not selected by segmentation quality.
 hh,ww=scaled.shape[:2];r=round(55*native_mpp/target_mpp);axes[i,2].imshow(overlay[hh//2-r:hh//2+r,ww//2-r:ww//2+r]);axes[i,2].set_title('Fixed central detail')
 for ax in axes[i]:ax.axis('off')
 rows.append(dict(dataset=ds,roi=name,hires_bbox=[x0,y0,x1,y1],native_mpp=native_mpp,model_input_mpp=target_mpp,resize_scale_x=size[0]/(x1-x0),resize_scale_y=size[1]/(y1-y0),prob_thresh=.05,nms_thresh=.3,n_nuclei=int(labels.max()),area_um2_quantiles=np.quantile(props.area_um2,[0,.1,.5,.9,1]).tolist(),status='pilot_pending_morphology_review',image=str(imgpath),resampling='bilinear; upscaling does not recover lost detail'))
 print(name,rows[-1],flush=True)
fig.tight_layout();fig.savefig(root/'figures/Cervilla_segmentation_pilot.png',dpi=180);plt.close(fig)
(root/'metadata/segmentation_pilot.json').write_text(json.dumps(rows,indent=2))
