from common import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,ndimage
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
out=ROOT/'Cervilla';a1=ad.read_h5ad(V1/'cell_annotation.h5ad',backed='r');a2=ad.read_h5ad(out/'cells.h5ad',backed='r')
m1=pd.read_parquet(V1/'inputs/Cervilla_bin_to_cell.parquet');m2=pd.read_parquet(out/'inputs/bin_to_cell.parquet');assert np.array_equal(m1.barcode,m2.barcode) and np.array_equal(m1.raw_counts,m2.raw_counts)
common=m1.image_covered.to_numpy()&m2.image_covered.to_numpy();rows=[];compositions=[];obslist=[]
for version,a,m,summarypath,excluded in [('v1',a1,m1,V1/'metadata/reconstruction_qc_summary.json',V1/'excluded_cell_metadata.parquet'),('v2',a2,m2,out/'metadata/reconstruction_qc_summary.json',out/'excluded_cell_metadata.parquet')]:
 s=json.loads(summarypath.read_text());o=a.obs.copy();outside=np.unique(m.loc[~common,'cell_label']);o['fully_within_common_bin_coverage']=~o.source_cell_id.astype(int).isin(outside)
 ex=pd.read_parquet(excluded);allobs=pd.concat([o,ex]);whole=~allobs.source_cell_id.astype(int).isin(outside);q=o[o.fully_within_common_bin_coverage].copy();obslist.append(q)
 d=dict(version=version,reconstructed_full=int((allobs.n_counts>0).sum()),assigned_bin_objects_full=len(allobs),zero_expression_objects_full=int(allobs.n_counts.eq(0).sum()),retained_full=len(o),low_information_full=int(o.low_information.sum()),common_bins=int(common.sum()),common_source_counts=int(m.loc[common,'raw_counts'].sum()),common_assigned_counts=int(m.loc[common&(m.cell_label.to_numpy()>0),'raw_counts'].sum()),common_reconstructed_fully_covered=int((whole&(allobs.n_counts>0)).sum()),common_retained_fully_covered=len(q),common_low_information=int(q.low_information.sum()),common_median_counts=float(q.n_counts.median()),common_median_genes=float(q.n_genes.median()),common_unknown=int(q.major_label.eq('Unknown').sum()),common_mixed=int(q.annotation_status.astype(str).str.contains('mixed').sum()))
 d['common_assigned_fraction']=d['common_assigned_counts']/d['common_source_counts'];d['common_low_information_fraction']=d['common_low_information']/len(q);rows.append(d)
 for field in ['major_label','subtype_label','annotation_status']:
  for label,n in q[field].value_counts().items():compositions.append(dict(version=version,field=field,label=label,n=int(n),fraction=float(n/len(q)),denominator=len(q)))
table=pd.DataFrame(rows);table.to_csv(ROOT/'metadata/Cervilla_v1_v2_comparison.tsv',sep='\t',index=False);pd.DataFrame(compositions).to_csv(ROOT/'metadata/Cervilla_v1_v2_composition.tsv',sep='\t',index=False)
jwrite(ROOT/'metadata/Cervilla_common_coverage.json',dict(source_bins=len(m1),common_image_covered_bins=int(common.sum()),excluded_from_common=int((~common).sum()),counts_comparison='Same original source bins, counted only in intersection of image-covered masks',cell_comparison='Retained/reconstructed cells whose entire assigned-bin set is inside common coverage; no cell ID pairing',D1_v1=None,caveat='Both image and annotation rules changed; label changes cannot all be attributed to image recovery. Broad identities in v2 are distinguished from unresolved subtypes.'))
fig,axs=plt.subplots(2,3,figsize=(15,9))
for ax,metric in zip(axs[0],['common_reconstructed_fully_covered','common_retained_fully_covered','common_assigned_fraction']):ax.bar(['v1','v2'],table[metric],color=['#999999','#4477aa']);ax.set_title(metric.replace('common_','').replace('_',' '))
for ax,col in zip(axs[1,:2],['n_counts','n_genes']):
 for q,v in zip(obslist,['v1','v2']):x=np.sort(q[col].to_numpy());ax.plot(x,np.arange(1,len(x)+1)/len(x),label=v)
 ax.set_xscale('log');ax.set_xlabel(col);ax.set_ylabel('retained-cell cumulative fraction');ax.legend()
comp=pd.DataFrame(compositions);z=comp[comp.field.eq('major_label')].pivot(index='label',columns='version',values='fraction').fillna(0);z.plot.bar(ax=axs[1,2]);axs[1,2].set_ylabel('fraction of common-coverage retained cells');fig.tight_layout();fig.savefig(ROOT/'metadata/Cervilla_v1_v2_summary.png',dpi=160,transparent=True);plt.close(fig)
# Side-by-side boundaries at the same physical coordinates; old labels never transferred.
t1=json.loads((V1/'metadata/Cervilla_transform.json').read_text());t2=json.loads((out/'metadata/transform.json').read_text());oldmask=sparse.load_npz(V1/'inputs/Cervilla_nuclei.npz');newmask=sparse.load_npz(out/'inputs/nuclei.npz');im=np.load(out/'inputs/model_image.npy',mmap_mode='r');fig,axs=plt.subplots(3,2,figsize=(10,14));roi=[]
for i,(label,box) in enumerate(t2['representative_model_bboxes'].items()):
 x0,y0,x1,y1=box;yy,xx=np.mgrid[y0:y1,x0:x1];fullx=(xx+.5)/t2['local_to_model_scale_xy'][0]-.5;fully=(yy+.5)/t2['local_to_model_scale_xy'][1]-.5;ox=np.rint((fullx*t1['fullres_to_hires_scale']+.5)*t1['hires_to_model_scale_xy'][0]-.5).astype(int);oy=np.rint((fully*t1['fullres_to_hires_scale']+.5)*t1['hires_to_model_scale_xy'][1]-.5).astype(int);om=np.asarray(oldmask[oy.ravel(),ox.ravel()]).reshape(xx.shape);nm=newmask[y0:y1,x0:x1].toarray();n=[]
 for ax,mask,v in zip(axs[i],[om,nm],['v1','v2']):
  raw=np.array(im[y0:y1,x0:x1]);raw[find_boundaries(mask,mode='inner')]=[20,230,60];ax.imshow(raw);ax.axis('off');n.append(len(np.unique(mask[mask>0])));ax.set_title(f'{label} {v}: {n[-1]} intersecting nuclei')
 roi.append(dict(roi=label,v1_nuclei_intersecting=n[0],v2_nuclei_intersecting=n[1],note='Same physical window on recovered-image background; boundary labels are separate reconstructions'))
fig.tight_layout();fig.savefig(ROOT/'metadata/Cervilla_v1_v2_segmentation.png',dpi=160,transparent=True);plt.close(fig);jwrite(ROOT/'metadata/Cervilla_v1_v2_roi_counts.json',roi);a1.file.close();a2.file.close();print(table.to_string(index=False),flush=True)
