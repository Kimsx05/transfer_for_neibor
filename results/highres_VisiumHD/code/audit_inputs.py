#!/usr/bin/env python3
"""Read 10x CSC matrices in bin blocks; audit exact counts and barcode-joined geometry.
No normalization, cell assignment, segmentation, or biological grouping is performed.
"""
from pathlib import Path
import os,json,hashlib,gc,re,time,sys
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'logs/matplotlib'))
import numpy as np,pandas as pd,h5py,pyarrow.parquet as pq
from scipy import sparse
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];SEED=123
GENES='HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8<<20):h.update(b)
 return h.hexdigest()
def decode(a):return np.asarray(a).astype(str)
def clean(x):
 if isinstance(x,(np.integer,np.floating)):return x.item()
 if isinstance(x,np.ndarray):return [clean(v) for v in x]
 if isinstance(x,bytes):return x.decode()
 return x

def specifications():
 jobs=[]
 for r in [2,8,16]:
  base=ROOT/f'data/GSE300692/D1_Johnson/outs/binned_outputs/square_{r:03d}um'
  jobs.append(dict(dataset='GSE300692',sample='GSM9066380',section='D1_Johnson',resolution_um=r,h5=base/'filtered_feature_bc_matrix.h5',spatial=base/'spatial',image=ROOT/'data/GSE300692/D1_Johnson/outs/spatial/tissue_hires_image.png'))
 base=ROOT/'data/Cervilla/visiumhd_Bladder/binned_outputs/square_002um'
 jobs.append(dict(dataset='Cervilla',sample='visiumhd_Bladder',section='visiumhd_Bladder',resolution_um=2,h5=base/'filtered_feature_bc_matrix.h5',spatial=base/'spatial',image=base/'spatial/tissue_hires_image.png'))
 base=ROOT/'data/GSE299573/GSM9040999/square_008um'
 jobs.append(dict(dataset='GSE299573',sample='GSM9040999',section='UC_SCC',resolution_um=8,h5=base/'filtered_feature_bc_matrix.h5',spatial=base/'spatial',image=base/'spatial/tissue_hires_image.png'))
 for sample,section in [('GSM9505917','MPBC1'),('GSM9505918','MPBC2')]:
  p=list((ROOT/'data/GSE318916').glob(sample+'*.h5'))[0]
  jobs.append(dict(dataset='GSE318916',sample=sample,section=section,resolution_um=2,h5=p,spatial=ROOT/'data/GSE318916'/section/'missing',image=ROOT/'data/GSE318916'/section/'missing.png'))
 return jobs

def plots(job,counts,genes,positions,scale,barcodes):
 """Use fixed seed for overview subsampling; neutral ROI names, no tissue calls."""
 out=ROOT/'figures'/job['dataset']/job['section'];out.mkdir(parents=True,exist_ok=True)
 tag=f"{job['resolution_um']}um";rng=np.random.default_rng(SEED)
 fig,axs=plt.subplots(1,2,figsize=(9,3.5))
 for ax,values,label in zip(axs,[counts,genes],['Gene-expression counts / bin','Detected genes / bin']):
  ax.hist(np.log1p(values),bins=65,color='#3a7891');ax.set_xlabel('log1p('+label+')');ax.set_ylabel('Number of matrix bins')
 fig.suptitle(job['section']+' | '+tag+' | all supplied bins; no QC filtering');fig.tight_layout();fig.savefig(out/f'{tag}_distributions.png',dpi=130);plt.close(fig)
 if positions is None:
  ids=rng.choice(len(barcodes),min(100000,len(barcodes)),replace=False)
  grid=np.array([[int(v) for v in b.split('_')[2:3]+[b.split('_')[3].split('-')[0]]] for b in barcodes[ids]],dtype=float)
  fig,ax=plt.subplots(figsize=(6,5));sc=ax.scatter(grid[:,1],grid[:,0],c=np.log1p(counts[ids]),s=.5,rasterized=True,cmap='viridis');ax.invert_yaxis();ax.set_aspect('equal');ax.set_xlabel('Barcode array column');ax.set_ylabel('Barcode array row');ax.set_title(job['section']+' | relative grid only; H&E registration unavailable');fig.colorbar(sc,ax=ax,label='log1p counts');fig.tight_layout();fig.savefig(out/f'{tag}_relative_grid.png',dpi=140);plt.close(fig);return
 if not job['image'].exists():return
 valid=np.flatnonzero(np.isfinite(positions.pxl_col_in_fullres.to_numpy()) & (positions.in_tissue.to_numpy()==1))
 ids=rng.choice(valid,min(100000,len(valid)),replace=False)
 sf=scale['tissue_hires_scalef'];x=positions.pxl_col_in_fullres.to_numpy()*sf;y=positions.pxl_row_in_fullres.to_numpy()*sf
 image=Image.open(job['image']);w,h=image.size; image.thumbnail((1800,1800));arr=np.asarray(image)
 effective=scale.get('microns_per_pixel',job['resolution_um']/scale['spot_diameter_fullres'])/sf
 fig,axs=plt.subplots(1,2,figsize=(11,5));extent=(0,w,h,0)
 for ax in axs:ax.imshow(arr,extent=extent);ax.set_xlim(0,w);ax.set_ylim(h,0);ax.set_aspect('equal');ax.set_xlabel('hires image x (pixels)');ax.set_ylabel('hires image y (pixels)')
 axs[0].set_title('Tissue preview');sc=axs[1].scatter(x[ids],y[ids],c=np.log1p(counts[ids]),s=.4,alpha=.65,cmap='viridis',vmin=0,vmax=np.quantile(np.log1p(counts[valid]),.99));axs[1].set_title('Counts aligned by barcode');fig.colorbar(sc,ax=axs[1],label='log1p counts',shrink=.7)
 roi=[]
 # Quantile targets are snapped to actual tissue bins. No morphology class is inferred.
 for qx,qy in [(.3,.3),(.6,.55),(.85,.8)]:
  tx=np.quantile(x[ids],qx);ty=np.quantile(y[ids],qy);k=ids[np.argmin((x[ids]-tx)**2+(y[ids]-ty)**2)];radius=max(50,250/effective);roi.append((x[k],y[k],radius))
 for i,(cx,cy,r) in enumerate(roi):axs[0].add_patch(Rectangle((cx-r,cy-r),2*r,2*r,fill=False,edgecolor='crimson',lw=1));axs[0].text(cx-r,cy-r,f'ROI{i+1}',color='crimson')
 for ax in axs:
  length=1000/effective;ax.plot([.05*w,.05*w+length],[.94*h,.94*h],color='black',lw=3);ax.text(.05*w,.92*h,'1 mm',fontsize=9)
 fig.suptitle(job['section']+' | '+tag+' | seed=123 overview subsample');fig.tight_layout();fig.savefig(out/f'{tag}_overview.png',dpi=130);plt.close(fig)
 fig,axes=plt.subplots(2,3,figsize=(12,8))
 for i,(cx,cy,r) in enumerate(roi):
  ax=axes[1,i]
  box=(max(0,int(cx-r)),max(0,int(cy-r)),min(w,int(cx+r)),min(h,int(cy+r)))
  with Image.open(job['image']) as im:crop=np.asarray(im.crop(box))
  raw=axes[0,i];raw.imshow(crop,extent=(box[0],box[2],box[3],box[1]));raw.set_title(f'ROI{i+1} | original hires pixels');raw.set_aspect('equal')
  ax.imshow(crop,extent=(box[0],box[2],box[3],box[1])); sel=valid[(abs(x[valid]-cx)<r)&(abs(y[valid]-cy)<r)]
  if len(sel)>25000:sel=rng.choice(sel,25000,replace=False)
  ax.scatter(x[sel],y[sel],c=np.log1p(counts[sel]),s=1,alpha=.45,cmap='viridis',vmin=sc.norm.vmin,vmax=sc.norm.vmax)
  ax.set_xlim(box[0],box[2]);ax.set_ylim(box[3],box[1]);ax.set_title(f'ROI{i+1} | identity unassigned');ax.set_aspect('equal');length=100/effective;xx=box[0]+.08*(box[2]-box[0]);yy=box[3]-.08*(box[3]-box[1]);ax.plot([xx,xx+length],[yy,yy],color='black',lw=3);ax.text(xx,yy-5,'100 um',fontsize=8)
 fig.tight_layout();fig.savefig(out/f'{tag}_ROIs.png',dpi=140);plt.close(fig)

def audit(job):
 start=time.time();key=f"{job['dataset']}_{job['section']}_{job['resolution_um']}um"
 print('START',key,flush=True)
 base={k:v for k,v in job.items() if k not in ['h5','spatial','image']};base.update(matrix_path=str(job['h5']),patient_id='unknown',platform='Visium HD',preservation='FFPE',measured_bin_um=2,postprocessing_resolution_um=job['resolution_um'])
 if not job['h5'].exists():return dict(base,status='missing_matrix'),[]
 with h5py.File(job['h5'],'r') as f:
  m=f['matrix'];nfeat,nbin=map(int,m['shape'][:]);bar=decode(m['barcodes'][:]);names=decode(m['features/name'][:]);ids=decode(m['features/id'][:]);types=decode(m['features/feature_type'][:]);ge=types=='Gene Expression';ptr=m['indptr'][:]
  base.update(status='audited',features=nfeat,bins=nbin,stored_nnz=len(m['data']),orientation='features_by_bins',sparse_format='CSC',gene_expression_features=int(ge.sum()),feature_types=json.dumps(dict(zip(*np.unique(types,return_counts=True))),default=clean),duplicate_symbols=int(pd.Index(names[ge]).duplicated().sum()),duplicate_feature_ids=int(pd.Index(ids[ge]).duplicated().sum()),duplicate_matrix_barcodes=int(pd.Index(bar).duplicated().sum()),h5_software_version=str(f.attrs.get('software_version','unknown')),chemistry=str(f.attrs.get('chemistry_description','unknown')),library_ids=json.dumps(clean(f.attrs.get('library_ids',[]))),counts_provenance='Space Ranger supplied filtered UMI counts; no inferred expression used',normalization='none applied; source H5 counts per GEO/Space Ranger provenance',sha256=sha(job['h5']))
  assert len(ptr)==nbin+1 and ptr[0]==0 and ptr[-1]==len(m['data']) and np.all(np.diff(ptr)>=0)
  posfile=job['spatial']/'tissue_positions.parquet';positions=None;scale={};tissue=None
  if posfile.exists():
   pos=pq.read_table(posfile).to_pandas();dup=int(pos.barcode.duplicated().sum());base.update(position_rows=len(pos),duplicate_spatial_barcodes=dup,duplicate_array_coordinates=int(pos.duplicated(['array_row','array_col']).sum()),duplicate_pixel_coordinates=int(pos.duplicated(['pxl_row_in_fullres','pxl_col_in_fullres']).sum()),all_spatial_in_tissue_bins=int((pos.in_tissue==1).sum()))
   assert dup==0,'Cannot unambiguously join duplicate position barcodes'
   pos=pos.set_index('barcode');matched=pd.Index(bar).isin(pos.index);positions=pos.reindex(bar);del pos;gc.collect()
   tissue=positions.in_tissue.to_numpy()==1
   base.update(barcode_match_rate=float(matched.mean()),missing_coordinates=int(positions[['pxl_row_in_fullres','pxl_col_in_fullres']].isna().any(axis=1).sum()),matrix_in_tissue_bins=int(tissue.sum()),coordinate_units='array index; full-resolution pixel row=y,column=x; hires=fullres*tissue_hires_scalef',section_match='same supplied output hierarchy; barcode join verified')
   scale=json.loads((job['spatial']/'scalefactors_json.json').read_text());base['scalefactors']=json.dumps(scale)
   # Fit affine grid->fullres using a deterministic sample and independently check residuals.
   rng=np.random.default_rng(SEED);take=rng.choice(nbin,min(nbin,10000),replace=False);p=positions.iloc[take];A=np.c_[p.array_row,p.array_col,np.ones(len(p))];B=p[['pxl_row_in_fullres','pxl_col_in_fullres']].to_numpy();ok=np.isfinite(A).all(1)&np.isfinite(B).all(1);coef=np.linalg.lstsq(A[ok],B[ok],rcond=None)[0];res=B[ok]-A[ok]@coef
   base.update(grid_to_fullres_affine=json.dumps(coef.tolist()),affine_residual_p99_pixel=float(np.quantile(np.linalg.norm(res,axis=1),.99)),affine_determinant=float(np.linalg.det(coef[:2,:])))
   if job['image'].exists():
    with Image.open(job['image']) as im:w,h=im.size;base.update(image_width=w,image_height=h,image_format=im.format)
    sf=scale['tissue_hires_scalef'];x=positions.pxl_col_in_fullres.to_numpy()*sf;y=positions.pxl_row_in_fullres.to_numpy()*sf;inframe=(x>=0)&(x<w)&(y>=0)&(y<h)
    base.update(image_path=str(job['image']),image_class='downsampled tissue_hires; not original H&E',in_tissue_image_bounds_rate=float(inframe[tissue].mean()),original_microns_per_pixel=scale.get('microns_per_pixel','unknown'),hires_microns_per_pixel=scale.get('microns_per_pixel',job['resolution_um']/scale['spot_diameter_fullres'])/sf,registration_status='numerical affine and bounds checked; see image overlays for morphology')
  else:base.update(barcode_match_rate='NA',coordinate_units='barcode-encoded relative array row/column only',registration_status='H&E/positions/scalefactors missing',matrix_in_tissue_bins='NA',section_match='unknown')
  counts=np.zeros(nbin,dtype=np.float64);detected=np.zeros(nbin,dtype=np.int32)
  mpidx=[np.flatnonzero((names==g)&ge) for g in GENES];g_total=np.zeros(20);g_nz=np.zeros(20,dtype=np.int64);g_tissue=np.zeros(20,dtype=np.int64);g_tcounts=np.zeros(20);negative=nonfinite=noninteger=explicitzeros=0;canonical=True
  for lo in range(0,nbin,100000):
   hi=min(nbin,lo+100000);a,b=int(ptr[lo]),int(ptr[hi]);v=m['data'][a:b];ii=m['indices'][a:b]
   assert np.all((ii>=0)&(ii<nfeat))
   negative+=int((v<0).sum());nonfinite+=int((~np.isfinite(v)).sum());noninteger+=int((v!=np.floor(v)).sum());explicitzeros+=int((v==0).sum())
   block=sparse.csc_matrix((v,ii,ptr[lo:hi+1]-a),shape=(nfeat,hi-lo));canonical=canonical and block.has_canonical_format;block.sum_duplicates();block.eliminate_zeros()
   x=block[ge,:];counts[lo:hi]=np.asarray(x.sum(axis=0)).ravel();detected[lo:hi]=x.getnnz(axis=0)
   for k,ix in enumerate(mpidx):
    if not len(ix):continue
    vals=np.asarray(block[ix,:].sum(axis=0)).ravel();g_total[k]+=vals.sum();g_nz[k]+=(vals>0).sum()
    if tissue is not None:mask=tissue[lo:hi];g_tissue[k]+=((vals>0)&mask).sum();g_tcounts[k]+=vals[mask].sum()
  base.update(negative_values=negative,nonfinite_values=nonfinite,noninteger_values=noninteger,explicit_zeros=explicitzeros,canonical_sparse=canonical,total_gene_expression_counts=float(counts.sum()),zero_count_bins=int((counts==0).sum()),counts_quantiles=json.dumps(np.quantile(counts,[0,.25,.5,.75,.9,.99,1]).tolist()),genes_quantiles=json.dumps(np.quantile(detected,[0,.25,.5,.75,.9,.99,1]).tolist()),quantile_probabilities='0,0.25,0.5,0.75,0.9,0.99,1')
  rows=[]
  for k,g in enumerate(GENES):
   ix=mpidx[k];rows.append({**{k:base[k] for k in ['dataset','sample','section','resolution_um']},'gene':g,'feature_present':bool(len(ix)),'feature_matches':len(ix),'ensembl_ids':';'.join(ids[ix]),'mapping_unambiguous':len(ix)==1,'any_nonzero_all_supplied_bins':bool(g_nz[k]),'nonzero_all_supplied_bins':int(g_nz[k]),'total_counts_all_supplied_bins':float(g_total[k]),'all_supplied_bins_denominator':nbin,'in_tissue_nonzero_bins':int(g_tissue[k]) if tissue is not None else 'NA','in_tissue_bin_denominator':int(tissue.sum()) if tissue is not None else 'NA','in_tissue_nonzero_fraction':float(g_tissue[k]/tissue.sum()) if tissue is not None and tissue.sum() else 'NA','in_tissue_total_counts':float(g_tcounts[k]) if tissue is not None else 'NA','tissue_status':'explicit in_tissue from barcode-joined parquet' if tissue is not None else 'unknown; filtered filename does not establish tissue membership'})
  # 2um is the descriptive overview for D1; existing 8/16um receive distributions only.
  if job['dataset']=='GSE300692' and job['resolution_um']!=2:
   j=job.copy();j['image']=ROOT/'nonexistent';plots(j,counts,detected,positions,scale,bar)
  else:plots(job,counts,detected,positions,scale,bar)
 base['elapsed_s']=time.time()-start
 print('DONE',key,nbin,base['total_gene_expression_counts'],flush=True)
 return base,rows

def main():
 p=ROOT.parent/'MP10_spatial_multisection_top20_vs_core19_v1/01_signatures/01_MP10_Top20_vs_Core19_gene_table.csv';f=pd.read_csv(p);frozen=f[f.DirectTop20].sort_values('MP10_rank').gene.tolist();lock={'source':str(p),'sha256':sha(p),'expected':GENES,'frozen':frozen,'exact_match':frozen==GENES};formal=ROOT.parent/'highres_ST/highres_spatial_audit_v3_GSE300692/01_MP10_lock/frozen_EpiMP10.tsv';ff=pd.read_csv(formal,sep='\t');formal_genes=ff[ff.top20].sort_values('MP10_rank').gene_symbol.tolist();lock.update(formal_lock_source=str(formal),formal_lock_sha256=sha(formal),formal_top20=formal_genes,formal_basis_source=str(ROOT.parent/'ST/mapping/07_EpiMP_gene_basis.qs'),formal_match=formal_genes==GENES);lock['exact_match']=lock['exact_match'] and lock['formal_match'];(ROOT/'metadata/mp10_lock.json').write_text(json.dumps(lock,indent=2))
 auditrows=[];genrows=[]
 for job in specifications():
  key=f"{job['dataset']}_{job['section']}_{job['resolution_um']}um";cache=ROOT/'metadata'/f'audit_{key}.json'
  inputs=[job['h5'],job['image'],job['spatial']/'tissue_positions.parquet',job['spatial']/'scalefactors_json.json']
  signature=[(str(p),p.stat().st_size,p.stat().st_mtime_ns) if p.exists() else (str(p),None,None) for p in inputs]
  signature=json.loads(json.dumps(signature))
  # Checkpoint avoids re-reading millions of bins after unrelated download completion.
  r=json.loads(cache.read_text()) if cache.exists() else {}
  if r.get('input_signature')==signature:a,g=r['matrix'],r['genes']
  else:
   a,g=audit(job)
   if a['status']=='audited':cache.write_text(json.dumps({'matrix':a,'genes':g,'input_signature':signature},indent=2,default=clean))
  auditrows.append(a);genrows.extend(g);pd.DataFrame(auditrows).to_csv(ROOT/'matrix_spatial_audit.tsv',sep='\t',index=False);pd.DataFrame(genrows).to_csv(ROOT/'mp10_feature_detection.tsv',sep='\t',index=False);gc.collect()
if __name__=='__main__':main()
