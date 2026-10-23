from common import *
import numpy as np,pandas as pd,sys
from scipy import sparse
from skimage.segmentation import find_boundaries
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
name=sys.argv[1];_,base,out=section(name);t=json.loads((out/'metadata/transform.json').read_text());im=np.load(out/'inputs/model_image.npy',mmap_mode='r');mask=sparse.load_npz(out/'inputs/nuclei.npz');b=pd.read_parquet(out/'inputs/bin_to_cell.parquet');p=pd.read_parquet(base/'spatial/tissue_positions.parquet',columns=['barcode','pxl_col_in_fullres','pxl_row_in_fullres']).set_index('barcode');xy=p.loc[b.barcode,['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();xy=(xy-t['fullres_crop_xyxy'][:2]+.5)*t['local_to_model_scale_xy']-.5
fig,axs=plt.subplots(3,2,figsize=(10,14));stats=[]
for (label,box),aa in zip(t['representative_model_bboxes'].items(),axs):
 cx=(box[0]+box[2])//2;cy=(box[1]+box[3])//2;x0,y0=cx-96,cy-96;x1,y1=cx+96,cy+96;raw=np.array(im[y0:y1,x0:x1]);m=mask[y0:y1,x0:x1].toarray();overlay=raw.copy();overlay[find_boundaries(m,mode='inner')]=[20,230,60];aa[0].imshow(overlay);aa[1].imshow(overlay)
 hit=(xy[:,0]>=x0)&(xy[:,0]<x1)&(xy[:,1]>=y0)&(xy[:,1]<y1);local=xy[hit]-[x0,y0];lab=b.cell_label.to_numpy()[hit];c=plt.get_cmap('hsv')((lab*.61803398875)%1);c[lab==0]=[.2,.2,.2,.35];aa[1].scatter(*local.T,s=12,marker='s',c=c,alpha=.45,linewidths=0)
 aa[0].set_title(label+': final nuclei, ~96um window');aa[1].set_title('Measured 2um centers colored by cell; gray=unassigned')
 for ax in aa:ax.set_xlim(0,192);ax.set_ylim(192,0);ax.axis('off')
 stats.append(dict(roi=label,model_bbox=[x0,y0,x1,y1],bins=int(hit.sum()),assigned=int((lab>0).sum())))
fig.tight_layout();fig.savefig(out/'figures/expansion_detail.png',dpi=160,transparent=True);plt.close(fig);jwrite(out/'metadata/expansion_detail_windows.json',stats);print('EXPANSION_DETAIL',name,flush=True)
