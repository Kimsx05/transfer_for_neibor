from pathlib import Path
import pandas as pd,numpy as np,h5py,anndata as ad,json,sys
from scipy import sparse
root=Path(__file__).resolve().parents[1];hd=root.parents[1];s=sys.argv[1];o=root/s
if (o/'roi_expression.tsv.gz').exists():sys.exit()
tr=json.loads((o/'transform.json').read_text());rois=pd.read_csv(o/'fixed_rois.tsv',sep='\t');cr=tr['fullres_crop_xyxy'];sc=tr['local_to_model_scale_xy'];mpp=tr['native_mpp'];genes=pd.read_csv(o/'gene_mapping.tsv',sep='\t');pan=pd.read_csv(hd/'analysis/hd_major_annotation_v1'/s/'tables/marker_panel_genes.tsv',sep='\t');panels={k:g.symbol.tolist() for k,g in pan.groupby('major')};groups={'Epithelial':panels['Epithelial'],'Fibroblast':panels['Fibroblast'],'Immune':sum([panels[k] for k in ['T/NK','Myeloid','B','Plasma','Mast']],[]),'Endothelial':panels['Endothelial']}
def region(x,y):
 mx=(np.asarray(x)/mpp-cr[0]+.5)*sc[0]-.5;my=(np.asarray(y)/mpp-cr[1]+.5)*sc[1]-.5;r=np.full(len(x),'',object)
 for b in rois.itertuples():r[(mx>=b.x)&(mx<b.X)&(my>=b.y)&(my<b.Y)]=b.roi
 return r
rows=[];r2=pd.read_parquet(o/'roi_2um_signals.parquet');rr=pd.DataFrame(dict(unit='raw_2um',object_id=r2.barcode,roi=r2.roi,x_um=r2.x_um,y_um=r2.y_um,nCount=r2.raw_counts,nFeature=np.nan,analysis=True))
for k in groups:
 cc=['T/NK','Myeloid','B','Plasma','Mast'] if k=='Immune' else [k];rr['raw_'+k]=r2[['raw_'+c for c in cc]].sum(axis=1);rr['detected_'+k]=r2[['detected_'+c for c in cc]].sum(axis=1)
rows.append(rr)
qc=pd.read_parquet(hd/'analysis/segmentation_assignment_v2'/s/'object_qc.parquet');idx=qc.index[qc.qc_keep_analysis];keep=set(idx)
for unit,dr in [('cell4','segmentation_restart_v1'),('cell8','segmentation_assignment_v2')]:
 a=ad.read_h5ad(hd/'analysis'/dr/s/'cells.h5ad',backed='r');obs=a.obs;xy=a.obsm['spatial'];rg=region(xy[:,0]*mpp,xy[:,1]*mpp);hit=(rg!='')&obs.index.isin(keep);ii=np.flatnonzero(hit);c=a.X[ii,:];meta=pd.DataFrame(dict(unit=unit,object_id=obs.index[ii],roi=rg[ii],x_um=xy[ii,0]*mpp,y_um=xy[ii,1]*mpp,nCount=np.asarray(c.sum(1)).ravel(),nFeature=c.getnnz(1),analysis=True));sy=a.var.gene_symbol.values
 for k,g in groups.items():jj=np.flatnonzero(np.isin(sy,g));meta['raw_'+k]=np.asarray(c[:,jj].sum(1)).ravel();meta['detected_'+k]=c[:,jj].getnnz(1)
 rows.append(meta);a.file.close()
b=pd.read_csv(o/'bin_metadata.tsv.gz',sep='\t');info=json.loads((o/'input_summary.json').read_text());hit=b.roi.notna();bi=np.flatnonzero(hit)
with h5py.File(hd/info['matrix_path']) as f:
 g=f['matrix'];shape=g['shape'][:];x=sparse.csr_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=(shape[1],shape[0]));c=x[bi];sy=g['features/name'].asstr()[:]
meta=pd.DataFrame(dict(unit='bin8',object_id=b.bin_id[hit],roi=b.roi[hit],x_um=b.x_um[hit],y_um=b.y_um[hit],nCount=b.nCount_RNA[hit],nFeature=b.nFeature_RNA[hit],analysis=b.qc_keep[hit]));meta=meta.reset_index(drop=True)
for k,g in groups.items():jj=np.flatnonzero(np.isin(sy,g));meta['raw_'+k]=np.asarray(c[:,jj].sum(1)).ravel();meta['detected_'+k]=c[:,jj].getnnz(1)
rows.append(meta);z=pd.concat(rows,ignore_index=True);e=z.raw_Epithelial>0;f=z.raw_Fibroblast>0;z['signal_category']=np.select([e&f,e,f],['Both panels','Epithelial panel','Collagen panel'],default='Neither panel');z['epi_fib_ge2']=(z.detected_Epithelial>=2)&(z.detected_Fibroblast>=2);z.to_csv(o/'roi_expression.tsv.gz',sep='\t',index=False)
# Native grid basis for actual square outlines, using original spatial metadata.
q=b.iloc[::300];coef=np.linalg.lstsq(np.c_[q.array_row,q.array_col,np.ones(len(q))],q[['x_um','y_um']],rcond=None)[0];(o/'grid8_basis_um.json').write_text(json.dumps(coef[:2].tolist()));stats=z.groupby(['roi','unit']).agg(n=('object_id','size'),UMI=('nCount','sum'),median_nFeature=('nFeature','median'),both_panels=('signal_category',lambda q:(q=='Both panels').sum()),epi_fib_ge2=('epi_fib_ge2','sum'));stats.to_csv(o/'roi_signal_summary.tsv',sep='\t');print('ROI_EXPRESSION_COMPLETE',s,flush=True)
