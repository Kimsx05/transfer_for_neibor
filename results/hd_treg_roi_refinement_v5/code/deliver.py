from common5 import *
from scipy.signal import fftconvolve
from matplotlib import pyplot as plt
from PIL import Image
import importlib.metadata as md
from importlib.util import find_spec

def collect(s):
 o=R/s; d=pd.read_parquet(o/'cells_scored.parquet'); c=pd.read_parquet(o/'candidate_ledger.parquet'); raw=pd.read_parquet(o/'raw_candidates.parquet')
 ar,au=inputdata(s); tr=js(N/s/'metadata/transform.json')
 mani=[(k,au[k]) for k in ['input_matrix','scale_source','source_image']]+[(str(i.relative_to(A)),str(i)) for i in (N/s/'metadata/transform.json',N/s/'metadata/input_audit.json',N/s/'inputs/model_image.npy',S/s/'inputs/nuclei.npz',S/s/'bin_to_cell.parquet',S/s/'object_qc.parquet',S/s/'cells.h5ad',U/'model.joblib',V/'code/common.py',V/s/'bin_evidence.parquet',V/s/'Epi64_members.parquet',A/'hd_annotation_diagnosis_singler_v1'/s/'diagnosis.tsv.gz')]+[('tissue_positions',str(p)) for p in (Path(au['input_matrix']).parent/'spatial').glob('tissue_positions*')]
 src=pd.DataFrame(mani,columns=['role','path']);src['readable']=src.path.map(lambda v:Path(v).is_file());src['section']=s;src['um_per_fullres_pixel']=au['native_mpp'];src['fullres_crop_xyxy']=str(tr['fullres_crop_xyxy'])
 assert c.unit.is_unique and d.id.is_unique and np.all(d.id==np.arange(1,len(d)+1)) and c.roi.gt(0).all() and d.loc[~d.score_ok,'MP10_cell_score'].isna().all()
 a=rd(o/'counts_audit.tsv.gz');x=sparse.load_npz(o/'counts.npz');assert np.array_equal(a.input,a.assigned+a.unassigned+a.excluded) and np.array_equal(np.asarray(x.sum(0)).ravel(),a.assigned) and np.array_equal(np.asarray(x.sum(1)).ravel(),d.nCount)
 # True native 2 µm supplied locations for coverage, including explicitly observed zero counts.
 cc=c[['unit','roi','x','y']].copy()
 for rad in (50,200):
  # Count sites at eligible observed grid centers. A tiled convolution reads occupied on-lattice sites;
  # query evaluation uses quantized center; near the .5 boundary do not use to confer eligibility.
  aff=np.array(au['array_col_row_affine_fullres']); ij=np.rint((cc[['x','y']].to_numpy()/au['native_mpp']-aff[2])@np.linalg.inv(aff[:2])).astype(int)
  col=ar.array_col.to_numpy(dtype=np.int64);row=ar.array_row.to_numpy(dtype=np.int64);mn=[min(col.min(),ij[:,0].min())-1,min(row.min(),ij[:,1].min())-1]
  xx=col-mn[0];yy=row-mn[1];ij=ij-mn;grid=np.zeros((max(yy.max(),ij[:,1].max())+2,max(xx.max(),ij[:,0].max())+2),np.float32)
  w=rad//2;rr,zz=np.ogrid[-w:w+1,-w:w+1];kernel=(rr*rr+zz*zz<=w*w).astype(np.float32)
  for tag,mask in [('supplied',np.ones(len(xx),bool)),('usable',ar.usable.values)]:
   grid[yy,xx]=mask; vals=fftconvolve(grid,kernel,mode='same')[ij[:,1],ij[:,0]]
   cc[f'{tag}_{rad}um_sites']=np.rint(vals).astype(int)
  cc[f'disk_{rad}um_sites']=int(kernel.sum())
 cc.to_parquet(o/'candidate_coverage.parquet',index=False)
 un=c.disposition.eq('no_assignable_cell').sum();z=pd.read_parquet(o/'zoom_index.parquet');r=c.sources.str.split('|',regex=False).explode().value_counts().rename_axis('source').reset_index(name='screened_units');r['section']=s
 parts=pd.DataFrame(dict(section=s,kind=['nuclei','cells_in_D','screen_entries','raw2_entries','no_cell_link','ROIs','epithelial','epithelial_scored','Treg','other_T','RNA_candidate','zoom_fields','zoom_pages'],n=[len(d),d.in_D.sum(),len(c),len(raw),un,c.roi.nunique(),d.epi_eligible.sum(),d.score_ok.sum(),(d.group=='Treg').sum(),(d.group=='other_T').sum(),(d.group=='RNA_candidate').sum(),len(z),z.file.nunique()]))
 p=d[d.in_D & d.pred.eq('Treg')];req=(p.T_width>=2)&(p.core_width>=1)
 stages={'first_reference_Treg':len(p),'plus_T_width_2_core_1':sum(req),'plus_reference_or_FOXP3_support':sum(req & ((p.support>=.95)|((p.FOXP3>0)&(p.core_width>=2)&(p.support>=.5)))),'final_after_geometry_mixture':(d.group=='Treg').sum()}
 pd.DataFrame({'stage':stages.keys(),'n':stages.values()}).to_csv(o/'identity_funnel.tsv',sep='\t',index=False)
 ep=d[d.score_ok];fig,ax=plt.subplots(1,2,figsize=(12,6),layout='constrained');u=d[d.in_D]
 ax[0].scatter(u.x,u.y,s=.17,c='#dddddd',rasterized=True);ax[0].scatter(ep.x,ep.y,s=.5,c='#5d8495',label='Scoreable epithelium',rasterized=True);h=ep[ep.H05];ax[0].scatter(h.x,h.y,s=2,c='#b94a34',label='H05',rasterized=True)
 for tag,color in [('RNA_candidate','#914b9b'),('other_T','#07826d')]:
  v=d[d.group.eq(tag)];ax[0].scatter(v.x,v.y,s=8,facecolors='none',edgecolors=color,linewidths=.6,label=tag)
 ax[0].legend(fontsize=8);ax[0].invert_yaxis();ax[0].set(aspect='equal',xlabel='x (µm)',ylabel='y (µm)')
 eff=rd(o/'effects.tsv'); e=eff[eff.metric.eq('H05_fraction')].copy()
 for j,v in enumerate(e.itertuples()):
  num=v.matched_block_mean_difference
  ax[1].text(.04,.85-j*.22,f'{v.group}   {v.radius} µm\n{v.matched_centers} matched centers, {v.matched_blocks} blocks\n'+('Effect unavailable' if pd.isna(num) else f'H05 difference {100*num:+.2f} pp, interval unavailable'),transform=ax[1].transAxes,fontsize=11)
 ax[1].axis('off')
 for ext in ['pdf','png']:fig.savefig(R/'figures'/f'{s}_section.{ext}',dpi=300)
 plt.close(fig)
 return src,r,parts

if __name__=='__main__':
 frames=[collect(s)for s in sections]
 for i,p in enumerate(['input_manifest.tsv','tables/candidate_sources.tsv','tables/section_summary.tsv']):put(pd.concat([v[i]for v in frames]),R/p)
 for name in ['thresholds','effects','ROI_manifest']:
  put(pd.concat([rd(R/s/(name+'.tsv'))for s in sections]),R/'tables'/(name+'.tsv'))
 versions=dict(python=sys.version.split()[0]);versions.update({v:md.version(v)for v in ['numpy','pandas','scipy','anndata','scikit-learn','numba','Pillow','h5py']});put(pd.DataFrame(versions.items(),columns=['software','version']),R/'tables/versions.tsv')
 atlas=[]
 for s in sections:
  for r in sorted(rd(R/s/'ROI_manifest.tsv').roi):
   for p in ['geometry','MP10']:atlas.append(Image.open(R/'figures'/f'{s}_ROI{r}_{p}.jpg').convert('RGB'))
 atlas[0].save(R/'figures/ROI_atlas.pdf',save_all=True,append_images=atlas[1:],resolution=120)
