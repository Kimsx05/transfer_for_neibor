# FROZEN_VERSION_GUARD: reproduce in a new version directory; never overwrite v1.
from pathlib import Path as _FreezePath
if (_FreezePath(__file__).resolve().parents[1] / 'metadata/FROZEN.json').exists():
    raise RuntimeError('Frozen v1 is immutable; create v2 for changes.')
"""Compact review panels from observed counts and final frozen labels; no physical-neighborhood analysis."""
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
from pathlib import Path
import json,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from PIL import Image
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from skimage.segmentation import find_boundaries
from sklearn.neighbors import NearestNeighbors
root=Path(__file__).resolve().parents[1];hd=root.parents[1];np.random.seed(123)
a=ad.read_h5ad(root/'cell_annotation.h5ad');o=a.obs;ex=pd.read_parquet(root/'excluded_cell_metadata.parquet');t=json.loads((root/'metadata/Cervilla_transform.json').read_text());im=Image.open(t['source_image']).convert('RGB');mpp=t['native_mpp'];xy=o[['x_um','y_um']].to_numpy()/mpp
palette={'Unknown':'#bdbdbd','major_only':'#4477aa','mixed_or_ambiguous':'#ee7733','subtype_supported':'#228833','Epithelial':'#d95f02','Fibroblast':'#7570b3','Plasma':'#e7298a','Myeloid':'#1b9e77','Endothelial':'#66a61e','Mast':'#e6ab02','T/NK':'#1f78b4','B':'#a6761d'}
palette_source=hd.parent/'subtype/v1/final_annotation_hierarchy_v2/03_palette/final_annotation_hierarchy_palette_master_v2.csv'
if palette_source.exists():
 p=pd.read_csv(palette_source)
 for _,row in p.iterrows():
  if row.final_fine_celltype not in palette:palette[row.final_fine_celltype]=row.hex_color
for i,k in enumerate(sorted(o.subtype_label.astype(str).unique())):
 if k not in palette:palette[k]='#cccccc' if k.startswith('unresolved') else to_hex(plt.get_cmap('tab20')(i%20))
(root/'metadata/review_palette.json').write_text(json.dumps(dict(colors=palette,existing_palette_source=str(palette_source),note='Existing fine-label colors reused only for exact matching names; no labels imported'),indent=2))

def cats(ax,coords,labels,title,s=.5,legend=True):
 labels=np.asarray(labels).astype(str);valid=np.isfinite(coords).all(1)
 order=sorted(set(labels[valid]),key=lambda z:(z not in ['Unknown','major_only','mixed_or_ambiguous'] and not z.startswith('unresolved'),z))
 for label in order:
  hit=valid&(labels==label);ax.scatter(coords[hit,0],coords[hit,1],s=s,color=palette.get(label,'#999999'),label=f'{label} ({hit.sum():,})',linewidths=0,rasterized=True)
 ax.set_title(title);ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 if legend:ax.legend(markerscale=5,fontsize=6,loc='upper left',bbox_to_anchor=(1,1),frameon=False)

def save(fig,name):fig.tight_layout();fig.savefig(root/'figures'/name,dpi=155,bbox_inches='tight');plt.close(fig);print(name,flush=True)
fig,axs=plt.subplots(2,3,figsize=(16,9))
for ax,col in zip(axs[0,:2],['n_counts','n_genes']):
 for frame,label,color in [(o,'Retained','#4477aa'),(ex,'Excluded','#cc6677')]:ax.hist(np.log10(frame[col]+1),bins=60,histtype='step',label=label,color=color)
 ax.axvline(np.log10(11),color='black',ls=':');ax.set(xlabel=f'log10({col}+1)',ylabel='Objects');ax.legend()
ax=axs[0,2];ax.scatter(o.n_genes,o.n_counts,s=.2,c=np.where(o.low_information,'#ee7733','#4477aa'));ax.set(xscale='log',yscale='log',xlabel='Observed genes',ylabel='Original counts',title='Low information in orange (retained)')
for ax,hit,title in [(axs[1,0],o.low_information.to_numpy(),'12,798 low-information cells'),(axs[1,1],o.qc_flags.ne('none').to_numpy(),'Flagged cells retained for review')]:
 ax.imshow(im);ax.scatter(xy[hit,0],xy[hit,1],s=.35,c='#e31a1c',linewidths=0);ax.set_title(title);ax.axis('off')
ax=axs[1,2];e=ex.exclusion_reason.value_counts();ax.barh(e.index,e.values,color='#cc6677');ax.set_title('Exclusion reasons');ax.tick_params(axis='y',labelsize=7)
save(fig,'Cervilla_QC_review.png')
fig,axs=plt.subplots(2,2,figsize=(17,12));u=a.obsm['X_umap'];core=o.scvi_training_core.to_numpy()
axs[0,0].scatter(u[core,0],u[core,1],s=.35,c='#4477aa',linewidths=0);axs[0,0].set_title('Sample: Cervilla only (84,985 structure-core cells)');axs[0,0].axis('off')
cats(axs[0,1],u,o.major_label,'Major identity (uncertain candidates remain visible)');cats(axs[1,0],u,o.subtype_label,'Final subtype / unresolved identity');cats(axs[1,1],u,o.annotation_status,'Annotation evidence status')
fig.suptitle('scVI expression UMAP: no cross-section integration claim; inference-only cells have no UMAP',fontsize=12);save(fig,'Cervilla_expression_UMAP_review.png')
fig,axs=plt.subplots(1,3,figsize=(22,7))
for ax,field,title in zip(axs,['major_label','subtype_label','annotation_status'],['Major identity','Subtype / unresolved','Uncertainty status']):
 cats(ax,xy,o[field],title,s=.5);ax.set_xlim(0,im.width);ax.set_ylim(im.height,0)
fig.suptitle('Cervilla: all 97,847 retained cells; native registered-image coordinates');save(fig,'Cervilla_spatial_labels_review.png')
# Direct raw marker availability and focused B/GC evidence.
genes=['KRT8','KRT19','KRT7','COL1A1','COL3A1','DCN','LUM','RGS5','ACTA2','TAGLN','PECAM1','VWF','CCL21','LYZ','C1QA','FOLR2','SPP1','CD3D','CD3E','TRAC','NKG7','MS4A1','CD79A','CD19','BCL6','AICDA','RGS13','JCHAIN','MZB1']
canon=a.var.reference_symbol.to_numpy();mat=np.column_stack([np.asarray(a.X[:,np.flatnonzero(canon==g)].sum(1)).ravel() for g in genes]);diagnostic=pd.DataFrame({'gene':genes,'shared_reference_features':[(canon==g).sum() for g in genes],'raw_detected_cells':(mat>0).sum(0),'raw_count_sum':mat.sum(0)});diagnostic.to_csv(root/'metadata/observed_marker_diagnostics.tsv',sep='\t',index=False)
gc=o.hierarchical_prediction.eq('GC B').to_numpy();groups=[(k,o.major_label.eq(k).to_numpy()&~o.annotation_status.eq('mixed_or_ambiguous').to_numpy()) for k in ['Epithelial','Fibroblast','Endothelial','Myeloid','T/NK','B','Plasma','Unknown']];groups.append(('GC B model candidates',gc))
fig,axs=plt.subplots(1,2,figsize=(19,8));det=np.array([(mat[h]>0).mean(0)*100 if h.any() else np.zeros(len(genes)) for _,h in groups]);v=axs[0].imshow(det,aspect='auto',vmin=0,vmax=100,cmap='Blues');axs[0].set_xticks(range(len(genes)),genes,rotation=90,fontsize=8);axs[0].set_yticks(range(len(groups)),[f'{k} (n={h.sum():,})' for k,h in groups],fontsize=8);axs[0].set_title('Observed raw-count detection %, shared reference genes');fig.colorbar(v,ax=axs[0],shrink=.5)
axs[1].imshow(im);axs[1].scatter(xy[gc,0],xy[gc,1],s=6,facecolors='none',edgecolors='#cc3311',linewidths=.5);b=o.major_label.eq('B').to_numpy();axs[1].scatter(xy[b,0],xy[b,1],s=25,c='#0077bb',marker='x');axs[1].set_title(f'GC B model candidates: {gc.sum()}; accepted GC B: 0\nB major: {b.sum()} (blue crosses); candidates are not labels');axs[1].axis('off');save(fig,'Cervilla_observed_markers_GC_review.png')
pd.DataFrame(det,index=[k for k,h in groups],columns=genes).to_csv(root/'metadata/marker_detection_review_percent.tsv',sep='\t')
# Whole available image coverage, accounting for every original bin.
p=hd/'data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um';bm=pd.read_parquet(root/'inputs/Cervilla_bin_to_cell.parquet');pos=pd.read_parquet(p/'spatial/tissue_positions.parquet',filters=[('in_tissue','==',1)]).set_index('barcode').loc[bm.barcode];bxy=pos[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()*t['fullres_to_hires_scale'];del pos
bins=[np.arange(0,im.width+51,50),np.arange(0,im.height+51,50)];total=np.histogram2d(*bxy.T,bins=bins,weights=bm.raw_counts)[0];missing=np.histogram2d(*bxy.T,bins=bins,weights=bm.raw_counts*(bm.cell_label==0))[0];fraction=np.divide(missing,total,out=np.full(total.shape,np.nan),where=total>0)
fig,axs=plt.subplots(1,3,figsize=(19,7));axs[0].imshow(im);axs[0].set_title('Native registered image / available canvas');axs[0].axis('off');axs[1].imshow(im);v=axs[1].imshow(fraction.T,extent=[0,bins[0][-1],bins[1][-1],0],cmap='magma',vmin=0,vmax=1,alpha=.8);axs[1].set_title('Unassigned raw-count fraction (50 native-pixel tiles)');axs[1].axis('off');fig.colorbar(v,ax=axs[1],shrink=.5)
axs[2].imshow(im);outside=~bm.image_covered.to_numpy();axs[2].scatter(bxy[outside,0],bxy[outside,1],s=1,c='red');axs[2].set_title(f'{outside.sum():,} in-tissue bins beyond image\nFull-section completion not claimed');axs[2].set_ylim(max(bxy[:,1].max()+20,im.height),0);axs[2].axis('off');save(fig,'Cervilla_coverage_unassigned_review.png')
# Expansion overlays in the same fixed representative ROIs; no new segmentation.
mask=sparse.load_npz(root/'inputs/Cervilla_nuclei.npz');sc=np.array(t['hires_to_model_scale_xy']);fig,axs=plt.subplots(3,3,figsize=(14,14));roi_stats=[]
for i,box in enumerate([[1500,2100,1800,2400],[2450,2700,2750,3000],[3600,3380,3900,3680]]):
 x0,y0,x1,y1=box;crop=np.asarray(im.crop(box));p0=np.rint((np.array(box[:2])+.5)*sc-.5).astype(int);p1=np.rint((np.array(box[2:])+.5)*sc-.5).astype(int);labels=mask[p0[1]:p1[1],p0[0]:p1[0]].toarray();over=np.asarray(Image.fromarray(crop).resize((labels.shape[1],labels.shape[0]),Image.Resampling.BILINEAR)).copy();over[find_boundaries(labels,mode='inner')]=[255,20,20];axs[i,0].imshow(crop);axs[i,1].imshow(over)
 hit=(bxy[:,0]>=x0)&(bxy[:,0]<x1)&(bxy[:,1]>=y0)&(bxy[:,1]<y1);loc=bxy[hit]-[x0,y0];lab=bm.loc[hit,'cell_label'].to_numpy();colors=plt.get_cmap('hsv')(((lab*0.618033988749895)%1));colors[lab==0]=[.3,.3,.3,.3];axs[i,2].imshow(crop);axs[i,2].scatter(loc[:,0],loc[:,1],s=3,marker='s',c=colors,alpha=.65,linewidths=0);axs[i,2].set_xlim(0,x1-x0);axs[i,2].set_ylim(y1-y0,0)
 for ax,title in zip(axs[i],['Native H&E','Final nuclear boundaries','Final bin assignment (gray = unassigned)']):ax.set_title(title,fontsize=9);ax.axis('off')
 roi_stats.append(dict(native_bbox=box,bins=int(hit.sum()),assigned_bins=int((lab>0).sum()),cells=int(np.unique(lab[lab>0]).size)))
fig.suptitle('Representative nuclei and cell expansion: <=2 array bins (4 um); image enlargement adds no detail');save(fig,'Cervilla_expansion_review.png');(root/'metadata/expansion_roi_review.json').write_text(json.dumps(roi_stats,indent=2))
# Descriptive expression-structure preservation check, independent of spatial coordinates.
z=np.load(root/'metadata/expression_embeddings.npz');coreidx=np.flatnonzero(core);co=o.iloc[coreidx];valid=~co.annotation_status.isin(['Unknown','mixed_or_ambiguous']);rng=np.random.default_rng(123);references=[];queries=[]
for label in sorted(co.major_label[valid].unique()):
 ix=np.flatnonzero(valid.to_numpy()&co.major_label.eq(label).to_numpy());rng.shuffle(ix)
 if len(ix)<20:continue
 nq=min(200,max(5,len(ix)//5));queries.extend(ix[:nq]);references.extend(ix[nq:nq+3000])
queries=np.array(queries);references=np.array(references);rows=[];ylab=co.major_label.astype(str).to_numpy()
for name,embedding in [('PCA_log_normalized_2000HVG',z['baseline_PCA_core']),('scVI_raw_count_2000HVG',z['X_scVI'][core])]:
 nn=NearestNeighbors(n_neighbors=15,n_jobs=8).fit(embedding[references]);_,neighbors=nn.kneighbors(embedding[queries]);agree=(ylab[references[neighbors]]==ylab[queries,None]).mean(1)
 for lab in sorted(set(ylab[queries])):
  h=ylab[queries]==lab;rows.append(dict(embedding=name,major_label=lab,n_queries=int(h.sum()),mean_same_major_neighbor_fraction=float(agree[h].mean()),interpretation='Descriptive classifier-supported subset, stratified held-out query cells; NOT annotation accuracy or batch-integration validation'))
pd.DataFrame(rows).to_csv(root/'metadata/expression_structure_review.tsv',sep='\t',index=False);print('ALL FIGURES DONE',flush=True)
