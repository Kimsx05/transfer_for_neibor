from base_data import *
import h5py,importlib.metadata as md,subprocess,platform
from anndata.io import read_elem
samples=pd.read_csv(V/'shared/samples.tsv',sep='\t');files=pd.read_csv(V/'shared/xenium_inputs.tsv',sep='\t')
rows=[]
for rec in samples.itertuples():
 s=rec.sample;row=dict(platform=rec.platform,sample_id=s,specimen=rec.specimen,patient_id='unknown',patient_evidence='v1 specimen manifest; no verified patient linkage',coordinates_unit='um')
 if s in XE:
  f=files[files['sample']==s].set_index('role').path
  z=pd.read_parquet(V/'xenium'/s/'labels.parquet');raw=f['expression'];old=str(V/'xenium'/s/'labels.parquet');coords=f['cell centroids + QC']
  with h5py.File(raw) as h:
   m=h['matrix'];ids=m['barcodes'].asstr()[:];genes=m['features/feature_type'].asstr()[:];assert np.array_equal(ids,z.index);n=len(ids);ng=int((genes=='Gene Expression').sum())
  row.update(raw=raw,format='Xenium 10x H5',n_units=n,n_GEX=ng,positions=coords,boundary=f['original cell polygons'],mask=str(V/'xenium'/s/'window.npz'),image_transform='original centroid/polygon coordinates in um; morphology stain image not used',labels=old,state_metadata=old,availability='complete counts/centroids/polygons',unit_version='original segmentation')
 elif s in HD:
  raw=src(s);old=H/s;t=json.loads((old/'metadata/transform.json').read_text())
  with h5py.File(raw) as h:ng,n=h['matrix/shape'][:]
  row.update(raw=str(raw),format='Space Ranger 2um H5',n_units=n,n_GEX=ng,positions=str(raw.parent/'spatial/tissue_positions.parquet'),boundary=str(old/'inputs/bin_to_cell.parquet'),image=t['source_image'],image_transform=str(old/'metadata/transform.json'),mask=str(V/'hd'/s/'window.npz'),labels=str(A/'cell_annotation_v2_1/cell_annotations.parquet'),state_metadata=str(V.parent/'mp10_hd_xenium_colocalization_v2/hd'/s),availability='complete; only observed count footprint',unit_version='old primary + new secondary; historical MP10 is on separate units')
 else:
  row.update(availability='missing original square_002um GEX, positions and verified registered image bundle; 8um cannot substitute',raw='',unit_version='not reconstructed')
 row['pathology_regions']='no verified pathologist-signed same-section layer used'
 rows.append(row)
write(pd.DataFrame(rows),'input_manifest.tsv')
v={p:md.version(p) for p in ['bin2cell','stardist','scvi-tools','torch','celltypist','cellcharter','scanpy','anndata','pandas','numpy','scipy','scikit-learn','tensorflow-cpu','pyarrow']};v.update(tacco='0.5.0 local target, writable NumPy vectors for pandas 3',python=sys.version,python_executable=sys.executable,reference=str(H/'metadata/reference_manifest.json'),crohn_author_commit=(O/'shared/crohn_notebook_settings.txt').read_text().splitlines()[0],heart_author_commit=(O/'shared/heart_notebook_settings.txt').read_text().splitlines()[0]);(O/'shared/versions.json').write_text(json.dumps(v,indent=2))
for cmd,f in [(['nvidia-smi'],'gpu_snapshot.txt'),(['free','-h'],'memory_snapshot.txt'),(['df','-h','.'],'disk_snapshot.txt'),(['lscpu'],'cpu_snapshot.txt'),(['Rscript','--version'],'R_version.txt')]:
 p=subprocess.run(cmd,capture_output=True,text=True);(O/'logs'/f).write_text(p.stdout+p.stderr)

