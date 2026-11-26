from common import *
from importlib.metadata import version
records=[];panels=[];checks=[]
for s in HD+XE:
 v,d=folders(s); z=labels(s);g=pd.read_parquet(v/'graph_nodes.parquet').set_index('id')
 assert g.index.equals(z.index) and z.index.is_unique and np.allclose(g[['x_um','y_um']],z[['x_um','y_um']])
 x,genes,ids=expression(s);assert np.array_equal(z.index,ids)
 assert np.array_equal(np.asarray(x.sum(1)).ravel(),z['nCount' if s in HD else 'transcript_counts'])
 if s in HD:
  m=pd.read_parquet(v/'members.parquet');w=pd.read_parquet(v/'weights.parquet');assert m.bin_id.is_unique and m.groupby('region_id').size().reindex(ids).eq(z.n_bins.values).all() and w.index.is_unique and set(w.index)<=set(ids)
  assert len(genes)==18082 and set(TOP)<=set(genes) and np.allclose(w.sum(1),1)
  p=imagedir(s)/'tissue_hires_image.png';files=[p,v/'members.parquet',v/'weights.parquet',*list((v/'counts').glob('*')),Path(json.loads((v/'conservation.json').read_text())['input'])]
 else:
  files=[Path(a) for a in pd.read_csv(V/'shared/xenium_inputs.tsv',sep='\t').query('sample==@s').path]+[v/'boundary_geometry.parquet',v/'identity.tsv.gz',v/'state_transfer.tsv.gz',v/'ucell.tsv.gz']
  q=pd.read_parquet(v/'boundary_geometry.parquet');assert set(q.index)==set(ids) and set(genes)==set((V/'shared/common_panel.txt').read_text().splitlines())
 for l,row in TAX.iterrows():
  a=row.anchors.split();panels.append((s,l,','.join(k for k in a if k in genes),','.join(k for k in a if k not in genes),'assessable' if set(a)&set(genes) else 'panel-limited'))
 for p in files+[v/'labels.parquet',v/'window.npz',v/'graph_nodes.parquet',v/'graph_50.npz']:
  records.append((s,str(p),p.exists(),p.stat().st_size if p.exists() else None,p.stat().st_mtime_ns if p.exists() else None))
 checks.append((s,len(z),x.shape[1],int(x.sum()),len(w) if s in HD else None,'passed'))
put(pd.DataFrame(records,columns=['sample','path','available','bytes','mtime_ns']),OUT/'audit/input_manifest.tsv')
put(pd.DataFrame(panels,columns=['sample','subtype','measured_anchors','missing_anchors','assessability']),OUT/'audit/panel_anchors.tsv')
put(pd.DataFrame(checks,columns=['sample','units','genes','UMI','fitted','alignment']),OUT/'audit/checks.tsv')
(OUT/'audit/versions.json').write_text(json.dumps({k:version(k) for k in ['squidpy','cellcharter','MENDER','torch','shapely','numpy','pandas','scipy','scikit-learn','pyarrow','numba']},indent=2))
