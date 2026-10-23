from common import *
import numpy as np,pandas as pd,anndata as ad,json,sys
morphology_only="--morphology-only" in sys.argv
requested=[n for n in ["Cervilla","D1_Johnson"] if n in sys.argv] or ["Cervilla","D1_Johnson"]
from scipy import sparse
from PIL import Image
import tifffile,zarr
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from skimage.segmentation import find_boundaries
from marker_rules import MAJOR,SUBTYPE
np.random.seed(123);Image.MAX_IMAGE_PIXELS=None;plt.rcParams.update({'figure.facecolor':'none','savefig.transparent':True,'axes.spines.top':False,'axes.spines.right':False})
palette=json.loads((ROOT/'metadata/review_palette.json').read_text())['colors'];palette.update(lineage_unresolved='#bdbdbd',subtype_unresolved='#4477aa',reference_broad_identity='#aa8866',mixed_or_ambiguous='#ee7733',subtype_supported='#228833',Unknown='#bdbdbd')
objects={n:ad.read_h5ad(ROOT/n/'cells.h5ad') for n in requested}
labels=sorted(set(x for a in objects.values() for x in a.obs.final_label.astype(str)))
for i,label in enumerate(labels):
 if label not in palette:palette[label]=palette.get(label.replace('_unresolved','').replace('T_NK','T/NK'),'#aaaaaa') if label.endswith('_unresolved') else to_hex(plt.get_cmap('tab20')(i%20))
jwrite(ROOT/'metadata/annotation_palette.json',palette)
def cats(ax,coords,labels,title,s=.35):
 valid=np.isfinite(coords).all(1);lab=np.asarray(labels).astype(str)
 for k in sorted(set(lab[valid]),key=lambda x:(x!='Unknown',x)):
  h=valid&(lab==k);ax.scatter(coords[h,0],coords[h,1],s=s,c=palette.get(k,'#999999'),label=f'{k}: {h.sum():,}',linewidths=0,rasterized=True)
 ax.set_title(title,fontsize=10);ax.set_aspect('equal');ax.axis('off');ax.legend(markerscale=6,fontsize=6,loc='upper left',bbox_to_anchor=(1,1),frameon=False)
def save(fig,path):fig.tight_layout();fig.savefig(path,dpi=160,transparent=True,bbox_inches='tight');plt.close(fig);print('FIGURE',path.name,flush=True)
for name,a in objects.items():
 out=ROOT/name;t=json.loads((out/'metadata/transform.json').read_text());o=a.obs;ex=pd.read_parquet(out/'excluded_cell_metadata.parquet');image=np.load(out/'inputs/model_image.npy',mmap_mode='r');mask=sparse.load_npz(out/'inputs/nuclei.npz');bm=pd.read_parquet(out/'inputs/bin_to_cell.parquet');_,base,_=section(name);pos=pd.read_parquet(base/'spatial/tissue_positions.parquet',columns=['barcode','pxl_col_in_fullres','pxl_row_in_fullres']).set_index('barcode');bxy=pos.loc[bm.barcode,['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();bmodel=(bxy-t['fullres_crop_xyxy'][:2]+.5)*t['local_to_model_scale_xy']-.5;del pos
 fig,axs=plt.subplots(3,3,figsize=(13,13));roirows=[]
 if name=='D1_Johnson':store=tifffile.imread(t['source_image'],aszarr=True);native=zarr.open(store,mode='r')
 else:native=Image.open(t['source_image'])
 for i,(label,box) in enumerate(t['representative_model_bboxes'].items()):
  x0,y0,x1,y1=box;raw=np.array(image[y0:y1,x0:x1]);lab=mask[y0:y1,x0:x1].toarray();overlay=raw.copy();overlay[find_boundaries(lab,mode='inner')]=[20,230,60];scales=np.tile(t['local_to_model_scale_xy'],2);nativebox=np.rint((np.array(box)+.5)/scales-.5).astype(int);nx0,ny0,nx1,ny1=nativebox
  original=np.array(native[ny0:ny1,nx0:nx1,:]) if name=='D1_Johnson' else np.asarray(native.crop(tuple(nativebox)))
  axs[i,0].imshow(original);axs[i,0].set_title(f'{label}: native H&E {t["native_mpp"]:.3f} um/px');axs[i,1].imshow(overlay);axs[i,1].set_title('Final nuclear boundaries (~0.5 um/px)');axs[i,2].imshow(raw)
  h=(bmodel[:,0]>=x0)&(bmodel[:,0]<x1)&(bmodel[:,1]>=y0)&(bmodel[:,1]<y1);loc=bmodel[h]-[x0,y0];labs=bm.cell_label.to_numpy()[h];colors=plt.get_cmap('hsv')((labs*.61803398875)%1);colors[labs==0]=[.3,.3,.3,.35];axs[i,2].scatter(*loc.T,s=4,marker='s',c=colors,alpha=.65,linewidths=0);axs[i,2].set_xlim(0,x1-x0);axs[i,2].set_ylim(y1-y0,0);axs[i,2].set_title('2um bin-to-cell; gray=unassigned')
  for ax in axs[i]:ax.axis('off')
  roirows.append(dict(roi=label,model_bbox=box,native_local_bbox=nativebox.tolist(),n_bins=int(h.sum()),assigned_bins=int((labs>0).sum()),n_labels=int(np.unique(labs[labs>0]).size)))
 if name=='D1_Johnson':store.close()
 save(fig,out/'figures/morphology_and_expansion.png');jwrite(out/'metadata/roi_review_stats.json',roirows)
 # Representative morphology flags are retained for review, not called confirmed doublets.
 geom=np.load(out/'inputs/nucleus_geometry.npz');adjudication=out/'metadata/morphology_adjudication.json';selection=[str(z['source_cell_id']) for z in json.loads(adjudication.read_text())['objects']] if adjudication.exists() else [];review_pool=pd.concat([o,ex],sort=False);largest=review_pool[review_pool.source_cell_id.astype(str).isin(selection)].sort_values('nuclear_area_um2',ascending=False) if selection else o.nlargest(6,'nuclear_area_um2');fig,axs=plt.subplots(2,3,figsize=(12,8));review=[]
 for ax,(cid,ob) in zip(axs.ravel(),largest.iterrows()):
  label=int(ob.source_cell_id);yy,xx=geom['points'][label-1];x0=max(0,int(xx)-120);y0=max(0,int(yy)-120);x1=min(image.shape[1],x0+240);y1=min(image.shape[0],y0+240);raw=np.array(image[y0:y1,x0:x1]);lab=mask[y0:y1,x0:x1].toarray();raw[find_boundaries(np.where(lab==label,label,0),mode='inner')]=[20,230,60];ax.imshow(raw);ax.axis('off');ax.set_title(f'{label}: nucleus {ob.nuclear_area_um2:.0f} um2\ncounts {ob.n_counts}; genes {ob.n_genes}; retained={cid in o.index}',fontsize=9);review.append(dict(cell_id=cid,source_cell_id=label,model_bbox=[x0,y0,x1,y1],reason='largest nucleus geometry, not evidence of confirmed error by itself'))
 save(fig,out/'figures/representative_large_nuclei.png');jwrite(out/'metadata/representative_morphology_objects.json',review)
 if morphology_only:continue
 # Full new source image is not embedded in the AnnData. Spatial plots use true microns.
 fig,axs=plt.subplots(2,3,figsize=(16,9))
 for ax,col in zip(axs[0,:2],['n_counts','n_genes']):
  for f,label,c in [(o,'retained','#4477aa'),(ex,'excluded','#cc6677')]:ax.hist(np.log10(f[col]+1),bins=60,histtype='step',label=label,color=c)
  ax.axvline(np.log10(11),ls=':',c='black');ax.set_xlabel('log10('+col+'+1)');ax.legend()
 axs[0,2].scatter(o.n_genes,o.n_counts,s=.2,c=np.where(o.low_information,'#ee7733','#4477aa'));axs[0,2].set(xscale='log',yscale='log',xlabel='genes',ylabel='raw counts',title=f'Low information retained: {o.low_information.sum():,}')
 axs[1,0].hist(np.log10(o.area_um2+1),bins=60,color='#4477aa',alpha=.6,label='cell bin area');axs[1,0].hist(np.log10(o.nuclear_area_um2+1),bins=60,histtype='step',color='#cc6677',label='nucleus');axs[1,0].set_xlabel('log10(area um2+1)');axs[1,0].legend();counts=ex.exclusion_reason.value_counts();axs[1,1].barh(counts.index,counts.values,color='#cc6677');axs[1,1].tick_params(axis='y',labelsize=6);axs[1,1].set_title('Exclusion reasons');summary=json.loads((out/'metadata/reconstruction_qc_summary.json').read_text());axs[1,2].bar(['assigned','unassigned'],[summary['assigned_counts'],summary['unassigned_counts']],color=['#4477aa','#bbbbbb']);axs[1,2].set_title('Observed 2um counts accounting');save(fig,out/'figures/QC.png')
 fig,axs=plt.subplots(1,3,figsize=(21,6))
 for ax,f in zip(axs,['major_label','final_label','annotation_status']):cats(ax,a.obsm['X_umap'],o[f],f+' (structure core only)')
 fig.suptitle(f'{name}: UMAP {o.umap_mapped.sum():,}/{len(o):,}; {len(o)-o.umap_mapped.sum():,} retained cells unmapped',fontsize=12);save(fig,out/'figures/expression_UMAP.png')
 fig,axs=plt.subplots(1,3,figsize=(21,7))
 for ax,f in zip(axs,['major_label','final_label','annotation_status']):cats(ax,a.obsm['spatial_um'],o[f],f+' (all retained cells)',s=.35);ax.invert_yaxis()
 fig.suptitle(f'{name}: all {len(o):,} retained cells; unmeasured image regions are not zero expression');save(fig,out/'figures/spatial_labels.png')
 # Focused programs; identity comes from expression, not these locations.
 focus=['Epithelial','Matrix_CAF','iCAF','Pericyte','SPP1 TAM','B','Plasma','T/NK','Endothelial'];fig,axs=plt.subplots(3,3,figsize=(13,13));xy=a.obsm['spatial_um']
 for ax,label in zip(axs.ravel(),focus):
  ax.scatter(*xy.T,s=.1,c='#dddddd',linewidths=0);hit=(o.major_label.eq(label) if label in MAJOR else o.final_label.eq(label))&~o.mixed_or_ambiguous;ax.scatter(*xy[hit].T,s=.6,c=palette.get(label,'#4477aa'),linewidths=0);ax.set_title(f'{label}: {hit.sum():,}',fontsize=9);ax.set_aspect('equal');ax.invert_yaxis();ax.axis('off')
 save(fig,out/'figures/focused_spatial.png')
 genes=['EPCAM','KRT8','KRT19','COL1A1','DCN','LUM','COL11A1','POSTN','MMP11','CTHRC1','FAP','PI16','DPT','CCL19','RGS5','MYH11','LYZ','LST1','C1QA','FOLR2','SPP1','GPNMB','MS4A1','CD79A','CD19','BCL6','AICDA','RGS13','MZB1','JCHAIN','CD3D','CD3E','TRAC','TRBC1','TRBC2','NKG7','GNLY','PECAM1','VWF','GJA5','EFNB2','HEY1','ACKR1','PROX1']
 canon=a.var.reference_symbol.to_numpy();values=np.column_stack([np.asarray(a.X[:,canon==g].sum(1)).ravel() for g in genes]);groups=[(z,o.final_label.eq(z).to_numpy()) for z in sorted(o.final_label.unique())];fig,ax=plt.subplots(figsize=(21,max(5,len(groups)*.29)))
 dd=[]
 for i,(label,h) in enumerate(groups):
  fraction=(values[h]>0).mean(0);mean=np.log1p(values[h]/np.maximum(o.n_counts.to_numpy()[h,None],1)*1e4).mean(0);pts=ax.scatter(np.arange(len(genes)),np.repeat(i,len(genes)),s=fraction*130+1,c=mean,cmap='viridis',vmin=0,vmax=3)
  for g,f,v in zip(genes,fraction,mean):dd.append(dict(label=label,n=int(h.sum()),gene=g,detection_fraction=f,mean_log1p_10k=v))
 ax.set_xticks(range(len(genes)),genes,rotation=90,fontsize=8);ax.set_yticks(range(len(groups)),[f'{z} (n={h.sum():,})' for z,h in groups],fontsize=8);ax.invert_yaxis();ax.set_title('Measured markers: area=detection fraction; color=mean log1p(10k)');fig.colorbar(pts,ax=ax,shrink=.5);save(fig,out/'figures/marker_dotplot.png');pd.DataFrame(dd).to_csv(out/'metadata/final_label_marker_evidence.tsv',sep='\t',index=False)
 del image,mask,bm,bxy,bmodel
if morphology_only:sys.exit(0)
# One joint graph display plus unintegrated PCA baseline; no tissue coordinates enter embedding.
fig,axs=plt.subplots(1,2,figsize=(12,5))
for name,a in objects.items():
 hit=a.obs.scvi_training_core.to_numpy();c='#4477aa' if name=='Cervilla' else '#cc6677';u=a.obsm['X_umap'][hit];p=a.obsm['X_pca_unintegrated'][hit];axs[0].scatter(*u.T,s=.15,c=c,label=f'{name} ({hit.sum():,})',linewidths=0);axs[1].scatter(p[:,0],p[:,1],s=.15,c=c,label=name,linewidths=0)
axs[0].set_title('Joint scVI UMAP (section_id batch)');axs[1].set_title('Unintegrated log-count PCA baseline')
for ax in axs:ax.legend(markerscale=10,fontsize=8);ax.set_xticks([]);ax.set_yticks([])
save(fig,ROOT/'metadata/two_section_expression_baseline.png')
