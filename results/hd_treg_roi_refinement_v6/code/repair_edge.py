from base6 import *
import subprocess,gzip,shapely
def repair(s):
 o=ROOT/s/'native_full';r=o/'edge_refit';r.mkdir(exist_ok=True);out=r/'out';out.mkdir(exist_ok=True)
 cmd=js(o/'command.json');d,au=inputdata(s);tr=js(N/s/'metadata/transform.json');bb=d.loc[d.usable,['pxl_col_in_fullres','pxl_row_in_fullres']].agg(['min','max']).to_numpy();bm=(bb-np.array(tr['fullres_crop_xyxy'][:2])+.5)*np.array(tr['local_to_model_scale_xy'])-.5;lo=np.ceil(bm[0]).astype(int);hi=np.floor(bm[1]).astype(int)
 lab=np.load(o/'mask.npy').copy();lab[:max(0,lo[1]),:]=0;lab[:,:max(0,lo[0])]=0;lab[max(0,hi[1]+1):,:]=0;lab[:,max(0,hi[0]+1):]=0;np.save(r/'mask.npy',lab);del lab,d
 cmd[cmd.index('--cellpose-masks')+1]=str(r/'mask.npy');cmd[cmd.index('--output-path')+1]=str(out);cmd[cmd.index('--nthreads')+1]='7'
 (r/'command.json').write_text(json.dumps(cmd,indent=1))
 with (r/'native.log').open('w') as f:subprocess.run(['/usr/bin/time','-v',*cmd],stdout=f,stderr=subprocess.STDOUT,check=True)
 # Preserve original interior fit. Use corrected fit at both edges of each measured coordinate axis.
 cm=pd.read_parquet(o/'out/cells.parquet');ncm=pd.read_parquet(out/'cells.parquet');shift=float((o/'native_offset_um.txt').read_text()) if (o/'native_offset_um.txt').exists() else 0
 pos=pd.read_parquet(ROOT/s/'cells_geometry.parquet',columns=['id','x','y']).set_index('id');bounds=bb*au['native_mpp']
 def interior(x):
  return ((x>bounds[0]+50)&(x<bounds[1]-50)).all(axis=1)
 switched=set(pos.index[~interior(pos[['x','y']].to_numpy())]);switched.update(cm.original_cell_id[~interior(cm[['centroid_y','centroid_x']].to_numpy()-shift)].astype(int))
 baseline=cm[~cm.original_cell_id.astype(int).isin(switched)].copy();fixed=ncm[ncm.original_cell_id.astype(int).isin(switched)].copy()
 masks=sparse.load_npz(ROOT/s/'final_nuclei.npz').tocsr();final=[];props=[];audit=[]
 for directory,sub,source in [(o/'out',baseline,'interior_original'),(out,fixed,'boundary_refit')]:
  prior=sub.set_index('cell');features=json.load(gzip.open(directory/'cells.geojson.gz'))['features']
  for f in features:
   old=f['properties']['cell']
   if old not in prior.index:continue
   meta=prior.loc[old].to_dict();i=int(meta['original_cell_id']);poly=shapely.geometry.shape(f['geometry'])
   # Full mask-coordinate overlap check, not centroid distance. Include buffered 0.35um for half-pixel touching.
   box=poly.bounds;p=np.array([[box[1],box[0]],[box[3],box[2]]])-shift;pix=(p/au['native_mpp']-tr['fullres_crop_xyxy'][:2]+.5)*tr['local_to_model_scale_xy']-.5
   t0=np.maximum(0,np.floor(pix[0]).astype(int)-1);t1=np.minimum(masks.shape[::-1],np.ceil(pix[1]).astype(int)+2)
   if (t1>t0).all():
    xx0,yy0=t0;xx1,yy1=t1;mm=masks[yy0:yy1,xx0:xx1].tocoo();ii=mm.data==i;native=((np.c_[mm.col[ii]+xx0,mm.row[ii]+yy0]+.5)/tr['local_to_model_scale_xy']-.5+tr['fullres_crop_xyxy'][:2])*au['native_mpp']+shift
    hit=int(shapely.intersects_xy(poly,native[:,1],native[:,0]).sum()) if len(native) else 0
   else:hit=0
   audit.append((i,source,hit))
   if not hit:continue
   j=len(final);f['properties']['cell']=j;final.append(f);meta['cell']=j;props.append(meta)
 finaldir=o/'registered';finaldir.mkdir(exist_ok=True)
 with gzip.open(finaldir/'cells.geojson.gz','wt') as fp:json.dump({'type':'FeatureCollection','features':final},fp,separators=(',',':'))
 pd.DataFrame(props).to_parquet(finaldir/'cells.parquet',index=False);pd.DataFrame(audit,columns=['id','fit','nuclear_pixels_in_polygon']).to_parquet(finaldir/'registration_audit.parquet',index=False)
 print(s,len(switched),len(baseline),len(fixed),len(final),flush=True)
if __name__=='__main__':repair(sys.argv[1])
