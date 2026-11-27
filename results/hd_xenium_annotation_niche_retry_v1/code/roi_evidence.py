from base_data import *
from scipy.ndimage import gaussian_filter
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt,gzip
from skimage.segmentation import find_boundaries
def save(fig,file):
 for ext in ['png','pdf']:fig.savefig(O/'figures'/(file+'.'+ext),dpi=165,bbox_inches='tight')
 plt.close(fig)
for s in XE:
 z=pd.read_parquet(paths(s)/'analysis_index.parquet');idx=pd.read_csv(O/'tables/xenium_review_ROI.tsv',sep='\t').query('sample==@s')
 p=next((V/'xenium/geo').glob('*_'+s+'_cell_boundaries*'))
 with gzip.open(p,'rb') as f:poly=pd.read_parquet(f)
 geom=pd.read_parquet(V/'xenium'/s/'boundary_geometry.parquet');assert geom.index.is_unique
 for r in idx.itertuples():
  hit=z[(z.x_um.between(r.x0_um,r.x0_um+140))&(z.y_um.between(r.y0_um,r.y0_um+140))]
  part=poly[poly.cell_id.isin(hit.index)]
  y=[];colors=[]
  for cid,obj in part.groupby('cell_id',sort=False):
   y.append(obj[['vertex_x','vertex_y']].values);colors.append({'fine_supported':'#228833','major_only':'#4477aa','mixed_or_suspect':'#ee7733','major_unresolved':'#aaaaaa'}[hit.at[cid,'review_level']])
  fig,ax=plt.subplots(figsize=(4.5,4.5));ax.add_collection(PolyCollection(y,facecolors=colors,edgecolors='#444444',linewidths=.4,alpha=.8));ax.autoscale_view();ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (µm)');ax.set_ylabel('y (µm)');ax.set_title(s+' · '+r.rationale+' · original segmentation');save(fig,s+'_ROI_'+r.rationale+'_boundaries')
for s in HD:
 d=paths(s);a=dataset(s);z=a.obs;b=pd.read_parquet(d/'members.parquet');u=z[z.source.eq('GEX_secondary')].sort_values(['n_counts','label']);t=json.loads((H/s/'metadata/transform.json').read_text());im=np.load(H/s/'inputs/model_image.npy',mmap_mode='r');tx=np.array(json.loads((d/'conservation.json').read_text())['grid_to_fullres']);rows=[]
 for rank,name in [(.25,'lower'),(.5,'median'),(.9,'upper')]:
  cell=u.iloc[round(rank*(len(u)-1))];target=b[b.final_label==cell.label];center=target[['array_row','array_col']].mean().values
  near=b[(np.abs(b.array_row-center[0])<=16)&(np.abs(b.array_col-center[1])<=16)]
  xy=np.c_[near.array_row,near.array_col,np.ones(len(near))]@tx
  mp=(xy-np.array(t['fullres_crop_xyxy'][:2])+.5)*t['local_to_model_scale_xy']-.5;x0,y0=np.maximum(mp.min(0).astype(int),0);x1,y1=mp.max(0).astype(int)+1
  fig,ax=plt.subplots(figsize=(4.4,4.4));ax.imshow(im[y0:y1,x0:x1],extent=(x0,x1,y1,y0));hit=near.final_label.eq(cell.label);ax.scatter(*mp[hit].T,s=24,facecolors='none',edgecolors='#f08700',linewidth=.9);ax.scatter(*mp[near.cell_label.gt(0)].T,s=1,c='#39b7b4',alpha=.35);ax.axis('off');ax.set_title(s+' · secondary '+name+' · '+str(int(cell.n_counts))+' UMI');save(fig,s+'_GEX_'+name)
  rows.append((s,name,cell.name,cell.n_counts,cell.nFeature,cell.area_um2))
 write(pd.DataFrame(rows,columns=['sample','selection','original_cell_id','counts','nFeature','area_um2']),s+'_GEX_examples.tsv')

