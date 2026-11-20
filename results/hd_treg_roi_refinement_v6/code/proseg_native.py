from base6 import *
from scipy.io import mmwrite
from threadpoolctl import threadpool_limits
import gzip, subprocess,time
def run(s,full=False):
 o=R/s/('native_full' if full else 'native_interface');o.mkdir(exist_ok=True);out=o/'out';out.mkdir(exist_ok=True)
 base=o/'in/square_002um';matrix=base/'raw_feature_bc_matrix';spa=base/'spatial';matrix.mkdir(exist_ok=True,parents=True);spa.mkdir(exist_ok=True)
 tr=js(N/s/'metadata/transform.json');au=js(N/s/'metadata/input_audit.json');d,_=inputdata(s);scale=np.array(tr['local_to_model_scale_xy']);org=np.array(tr['fullres_crop_xyxy'][:2]);mpp=tr['native_mpp'];mask=sparse.load_npz(R/s/'final_nuclei.npz').tocsr()
 sel=d.usable.to_numpy().copy()
 if full:
  x0=y0=0;y1,x1=mask.shape
 else:
  q=rd(R/'tables/calibration_fields.tsv').query('section==@s').iloc[0];cx,cy=int(q.mx),int(q.my);x0,y0=cx-200,cy-200;x1,y1=cx+200,cy+200
  mp=(d[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy()-org+.5)*scale-.5
  sel&=(mp[:,0]>=x0)&(mp[:,0]<x1)&(mp[:,1]>=y0)&(mp[:,1]<y1)
 idx=np.flatnonzero(sel);assert len(idx)>0
 with h5py.File(au['input_matrix']) as h:
  g=h['matrix'];genes=g['features/id'][:].astype(str);ptr=g['indptr'][:];allbc=g['barcodes'][:].astype(str);bar=allbc[idx];del allbc
  zz=sparse.csr_matrix((g['data'][:],g['indices'][:],ptr),shape=(len(d),len(genes)))[idx].T.tocoo()
 # Export exact supplied count entries and their original capture-position coordinates.
 pd.DataFrame({0:genes,1:genes,2:'Gene Expression'}).to_csv(matrix/'features.tsv.gz',sep='\t',header=False,index=False);pd.Series(bar).to_csv(matrix/'barcodes.tsv.gz',header=False,index=False)
 with gzip.open(matrix/'matrix.mtx.gz','wb',compresslevel=1) as f:
  with threadpool_limits(limits=2):mmwrite(f,zz,field='integer',symmetry='general')
 offset=((np.array([x0,y0])+.5)/scale-.5+org)*mpp
 translation=10. if offset.min()<0 else 0.
 if translation:(o/'native_offset_um.txt').write_text(str(translation))
 import pyarrow.parquet as pq
 pq.write_table(pa.table(dict(barcode=pa.array(bar,type=pa.string()),pxl_row_in_fullres=pa.array(d.pxl_row_in_fullres.to_numpy()[idx]+translation/mpp,type=pa.float64()),pxl_col_in_fullres=pa.array(d.pxl_col_in_fullres.to_numpy()[idx]+translation/mpp,type=pa.float64()))),spa/'tissue_positions.parquet')
 (spa/'scalefactors_json.json').write_text(json.dumps({'microns_per_pixel':mpp}));np.save(o/'mask.npy',mask[y0:y1,x0:x1].astype('uint32').toarray());np.save(o/'source_idx.npy',idx)
 off=offset+translation;step=mpp/scale
 cmd=[str(R/'.tools/prebuilt/bin/proseg'),'--visiumhd','--cellpose-masks',str(o/'mask.npy'),'--cellpose-x-transform','0',str(step[1]),str(off[1]),'--cellpose-y-transform',str(step[0]),'0',str(off[0]),'--nthreads',('12' if full else '8'),'--output-path',str(out),'--output-cell-metadata','cells.parquet','--output-cell-polygons','cells.geojson.gz','--no-diffusion','--exclude-spatialdata-transcripts']
 if full:cmd+=['--burnin-samples','24','--samples','24','--recorded-samples','12','--hillclimb','6','-m','150','--nhidden','32']
 if not full:cmd+=['--burnin-samples','3','--samples','3','--recorded-samples','2','--hillclimb','1','-m','20','--nhidden','16']
 cmd+=[str(base.parent)];(o/'command.json').write_text(json.dumps(cmd,indent=1));(o/'input_summary.json').write_text(json.dumps(dict(bins=len(idx),genes=len(genes),UMI=int(zz.sum()),mask_unique=int(np.unique(mask[y0:y1,x0:x1].data).size),mpp=mpp)))
 with (o/'native.log').open('w')as f:
  re=subprocess.run(['/usr/bin/time','-v',*cmd],stdout=f,stderr=subprocess.STDOUT);print(s,re.returncode,flush=True)
if __name__=='__main__':
 for s in [arg for arg in sys.argv[1:] if arg in sections] or sections:run(s,'full' in sys.argv)
