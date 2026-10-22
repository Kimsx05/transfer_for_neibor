# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Segment the available registered H&E canvas after the fixed pilot review.
The missing bottom image region remains uncovered; this is not full-section recovery.
"""
import os
from pathlib import Path
root=Path(__file__).resolve().parents[1];hd=root.parents[1]
os.environ.update(KERAS_HOME=str(root/'models/keras_cache'),CUDA_VISIBLE_DEVICES='-1',TF_NUM_INTRAOP_THREADS='8',TF_NUM_INTEROP_THREADS='2',OMP_NUM_THREADS='8')
if hasattr(os,'sched_setaffinity'):os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:12])
import numpy as np,json,random,tensorflow as tf
from PIL import Image
from scipy import sparse
from stardist.models import StarDist2D
from csbdeep.utils import normalize
random.seed(123);np.random.seed(123);tf.random.set_seed(123)
imgpath=hd/'data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um/spatial/tissue_hires_image.png'
sf=json.loads((imgpath.parent/'scalefactors_json.json').read_text());mpp=sf['microns_per_pixel']/sf['tissue_hires_scalef'];target=.5
out=root/'inputs/Cervilla_nuclei.npz'
if out.exists():print('Existing segmentation retained');raise SystemExit
im=Image.open(imgpath).convert('RGB');w,h=im.size;size=(round(w*mpp/target),round(h*mpp/target));sx,sy=size[0]/w,size[1]/h
transform=dict(source_image=str(imgpath),source_image_shape_yx=[h,w],model_image_shape_yx=[size[1],size[0]],native_mpp=mpp,model_nominal_mpp=target,actual_model_mpp_xy=[mpp/sx,mpp/sy],fullres_to_hires_scale=sf['tissue_hires_scalef'],hires_to_model_scale_xy=[sx,sy],hires_to_model_offset_xy=[.5*sx-.5,.5*sy-.5],rotation_degrees_added=0,flip_added=False,crop_offset_xy=[0,0],formula='model_xy=(fullres_xy*tissue_hires_scalef+0.5)*resize_scale_xy-0.5',original_image_recovered=False,coverage_status='available_image_canvas_only; 1274 filtered in-tissue bins beyond bottom cannot be reconstructed',normalization='csbdeep percentile 1/99.8 per RGB channel',seed=123)
(root/'metadata/Cervilla_transform.json').write_text(json.dumps(transform,indent=2))
scaled=np.asarray(im.resize(size,Image.Resampling.BILINEAR));del im
norm=normalize(scaled,1,99.8,axis=(0,1));del scaled
model=StarDist2D.from_pretrained('2D_versatile_he')
print('START canvas',norm.shape,flush=True)
labels,details=model.predict_instances_big(norm,axes='YXC',block_size=2048,min_overlap=128,context=128,prob_thresh=.05,nms_thresh=.3)
del norm
sparse.save_npz(out,sparse.csr_matrix(labels),compressed=True)
np.savez_compressed(root/'inputs/Cervilla_nucleus_geometry.npz',points=details['points'],prob=details['prob'])
(root/'metadata/Cervilla_segmentation_summary.json').write_text(json.dumps(dict(n_nuclei=int(labels.max()),mask_shape=list(labels.shape),pilot_status='algorithmic_morphology_review_pass_with_resolution_limit',model='2D_versatile_he',prob_thresh=.05,nms_thresh=.3,block_size=2048,min_overlap=128,context=128,cell_assignment_status='not_yet_aggregated'),indent=2))
print('DONE nuclei',labels.max(),flush=True)
