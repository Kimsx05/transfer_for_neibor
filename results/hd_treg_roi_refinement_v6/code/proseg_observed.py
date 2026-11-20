from base6 import *
from classifier import predict,to_ref
import gzip,shapely
def run(s):
 root=ROOT;out=R/s;assert R!=root;out.mkdir(exist_ok=True,parents=True)
 for f in ['final_nuclei.npz','nuclei_v6.npz','priority_old.parquet','nuclear_audit_old.parquet','nuclear_overlap.parquet','observed_panel.npz']:
  link=out/f
  if not link.exists():link.symlink_to(root/s/f)
 native=root/s/'native_full/registered';poly=json.load(gzip.open(native/'cells.geojson.gz'))['features'];cm=pd.read_parquet(native/'cells.parquet');cm['id']=cm.original_cell_id.astype(int);assert cm.id.is_unique
 shift=float((root/s/'native_full/native_offset_um.txt').read_text()) if (root/s/'native_full/native_offset_um.txt').exists() else 0;cm.centroid_x-=shift;cm.centroid_y-=shift
 ci=cm.set_index('cell').id;owner=ci.loc[[p['properties']['cell'] for p in poly]].to_numpy();geoms=np.array([shapely.geometry.shape(p['geometry']) for p in poly]);tree=shapely.STRtree(geoms)
 d,audit=inputdata(s);hard=np.zeros(len(d),np.int32);unc=np.zeros(len(d),bool)
 for lo in range(0,len(d),50000):
  hi=min(lo+50000,len(d));block=d.iloc[lo:hi];idx=np.flatnonzero(block.usable);pos=shapely.points(block.y.to_numpy()[idx]+shift,block.x.to_numpy()[idx]+shift);i,j=tree.query(pos,predicate='intersects');multiplicity=np.bincount(i,minlength=len(idx));good=multiplicity[i]==1;hard[lo+idx[i[good]]]=owner[j[good]];unc[lo+idx[multiplicity>1]]=True
 ap=pd.read_parquet(root/s/'assignment.parquet');ap['conservative_id']=ap.id;ap['conservative_ambiguous']=ap.ambiguous;ap['id']=hard;ap['ambiguous']=unc;ap.to_parquet(out/'assignment.parquet',index=False)
 cells=pd.read_parquet(root/s/'cells_geometry.parquet',columns=['id','x','y','nuclear_area_um2','overlap_old','nuclear_action','nuclear_bins']);cells['proseg_modelled']=cells.id.isin(cm.id);assert np.isin(np.unique(hard[hard>0]),cells.id).all()
 cells['proseg_centroid_shift_um']=np.hypot(cells.id.map(cm.set_index('id').centroid_y)-cells.x,cells.id.map(cm.set_index('id').centroid_x)-cells.y)
 n=len(cells);genes,sy=genedata(audit);ng=len(genes);X=sparse.csr_matrix((n,ng),dtype='int64');use=d.usable.to_numpy();acc=np.zeros((4,ng),np.int64)
 for a,b,z in slabs(audit):
  v=hard[a:b];ii=np.flatnonzero(v);X+=sparse.csr_matrix((np.ones(len(ii),int),(v[ii]-1,ii)),shape=(n,len(v)))@z
  for j,f in enumerate([np.ones(len(v),bool),~use[a:b],unc[a:b],(v==0)&use[a:b]&~unc[a:b]]):acc[j]+=np.asarray(z[f].sum(axis=0)).ravel()
 assigned=np.asarray(X.sum(axis=0)).ravel();assert np.array_equal(acc[0],acc[1:].sum(axis=0)+assigned)
 put(pd.DataFrame(dict(gene=genes,total=acc[0],excluded=acc[1],ambiguous=acc[2],unassigned=acc[3],assigned=assigned)),out/'count_balance.tsv.gz');sparse.save_npz(out/'counts.npz',X)
 cells['nCount']=np.asarray(X.sum(axis=1)).ravel();cells['nFeature']=X.getnnz(1)
 for g in T+reg+epi+fib+mye:cells[g]=np.asarray(X[:,sy==g].sum(axis=1)).ravel()
 for key,gs in [('direct',direct),('other_core',core[1:]),('T',T),('Epi',epi),('Fib',fib),('Myeloid',mye)]:cells[key+'_width']=cells[gs].gt(0).sum(axis=1)
 model=joblib.load(R/'model.joblib');q=predict(model,to_ref(X,sy,model[0]));cells=pd.concat([cells,q],axis=1)
 raw=pd.read_parquet(root/s/'raw_signal.parquet');raw['id']=hard[raw.source_bin_index];raw['ambiguous']=unc[raw.source_bin_index];raw.to_parquet(out/'raw_signal.parquet',index=False)
 foreign=raw[(raw.id>0)&(raw.nucleus_id>0)&raw.id.ne(raw.nucleus_id)&raw[direct+core].sum(axis=1).gt(0)];cells['foreign_key_bins']=cells.id.map(foreign.id.value_counts()).fillna(0).astype(int)
 cells.to_parquet(out/'cells_geometry.parquet',index=False);print(s,len(cm),int(assigned.sum()),flush=True)
if __name__=='__main__':
 for s in sys.argv[1:] or sections:run(s)
