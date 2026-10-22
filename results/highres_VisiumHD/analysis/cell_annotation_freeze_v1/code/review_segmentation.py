# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from PIL import Image
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
root=Path(__file__).resolve().parents[1];t=json.loads((root/'metadata/Cervilla_transform.json').read_text());im=Image.open(t['source_image']).convert('RGB');mask=sparse.load_npz(root/'inputs/Cervilla_nuclei.npz');a=ad.read_h5ad(root/'cell_annotation.h5ad',backed='r');sc=np.array(t['hires_to_model_scale_xy']);geom=np.load(root/'inputs/Cervilla_nucleus_geometry.npz')

def overlay(box,only=None):
 x0,y0,x1,y1=box;crop=np.asarray(im.crop(box));p0=np.rint((np.array([x0,y0])+.5)*sc-.5).astype(int);p1=np.rint((np.array([x1,y1])+.5)*sc-.5).astype(int);p0=np.maximum(p0,0);p1=np.minimum(p1,np.array(t['model_image_shape_yx'][::-1]));labels=mask[p0[1]:p1[1],p0[0]:p1[0]].toarray();arr=np.asarray(Image.fromarray(crop).resize((labels.shape[1],labels.shape[0]),Image.Resampling.BILINEAR)).copy();bound=find_boundaries(labels if only is None else np.where(labels==only,labels,0),mode='inner');arr[bound]=[255,20,20];return crop,arr
fig,axes=plt.subplots(3,2,figsize=(12,16))
for i,box in enumerate([[1500,2100,1800,2400],[2450,2700,2750,3000],[3600,3380,3900,3680]]):
 raw,over=overlay(box);axes[i,0].imshow(raw,interpolation='nearest');axes[i,1].imshow(over);axes[i,0].set_title('Native registered image');axes[i,1].set_title('Final whole-canvas StarDist nuclei')
 for ax in axes[i]:ax.axis('off')
fig.tight_layout();fig.savefig(root/'figures/Cervilla_final_nuclei_review.png',dpi=150);plt.close(fig)
selected=a.obs.nlargest(6,'nuclear_area_um2');fig,axes=plt.subplots(3,2,figsize=(12,16));rows=[]
for ax,(cell,o) in zip(axes.ravel(),selected.iterrows()):
 label=int(o.source_cell_id);my,mx=geom['points'][label-1];x,y=(np.array([mx,my])+.5)/sc-.5;r=55;box=[max(0,int(x-r)),max(0,int(y-r)),min(im.width,int(x+r)),min(im.height,int(y+r))];raw,over=overlay(box,only=label);ax.imshow(over);ax.set_title(f'label {label}; nucleus {o.nuclear_area_um2:.1f} um2\ncounts {o.n_counts}, genes {o.n_genes}');ax.axis('off');rows.append(dict(cell_id=cell,source_cell_id=label,nuclear_area_um2=o.nuclear_area_um2,hires_bbox=box,review_status='pending_visual_review'))
fig.tight_layout();fig.savefig(root/'figures/Cervilla_largest_nuclei_review.png',dpi=150);plt.close(fig);(root/'metadata/morphology_review_objects.json').write_text(json.dumps(rows,indent=2))
