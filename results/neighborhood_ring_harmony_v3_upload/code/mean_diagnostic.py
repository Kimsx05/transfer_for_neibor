from base import *
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
m,f,ref=inputs();qs=[np.zeros_like(f['q1'],dtype=float),np.zeros_like(f['q2'],dtype=float)];rows=[]
for (ds,s),gg in m.groupby(['dataset_id','section_id']):
 node=pd.read_csv(OLD/'02_graphs'/ds/(s+'__nodes.tsv'),sep='\t');raw=pd.read_csv(ROOT/'cell2location_full_run/v2_30type/02_sections'/s/'abundance_mean.csv.gz',index_col=0);raw.columns=[c.replace('meanscell_abundance_w_sf_','').replace('means_cell_abundance_w_sf_','') for c in raw.columns];a=raw.reindex(index=node.barcode,columns=TYPES).to_numpy(float);valid=np.isfinite(a).all(1)&(a.sum(1)>0);p=np.zeros_like(a);p[valid]=a[valid]/a[valid].sum(1,keepdims=True);ix=pd.Index(node.spot_uid).get_indexer(gg.spot_uid);assert (ix>=0).all()
 for i,ring in enumerate(['R1','R2']):
  g=sparse.load_npz(OLD/'02_graphs'/ds/(s+'__'+ring+'.npz'));den=g@valid.astype(float);q=(g@p)/np.maximum(den[:,None],1);assert (den[ix]>0).all();qs[i][gg.index]=q[ix]
x=np.concatenate([np.sqrt(qs[0]),np.sqrt(qs[1])],1)/np.sqrt(2);z=PCA(n_components=.95,svd_solver='full',whiten=False).fit(x[ref]).transform(x)
for tag,zz in [('q05',np.load(OUT/'01_features/R12_PCA.npy')),('mean',z)]:
 ids=NearestNeighbors(n_neighbors=31,n_jobs=1).fit(zz[ref]).kneighbors(zz[ref],return_distance=False)[:,1:]
 for key in ['dataset_id','sample_uid']:
  a=m.iloc[ref][key].to_numpy();rows.append(dict(environment=tag,key=key,n=len(ref),k=30,n_pc=zz.shape[1],same_source_fraction=float((a[ids]==a[:,None]).mean()),scope='one R12 source-structure diagnostic; primary remains q05'))
tab(OUT/'04_comparison/q05_vs_mean_source_structure.tsv',rows);js(OUT/'04_comparison/q05_vs_mean_manifest.json',dict(primary_changed=False,anchor_UIDs_unchanged=True,reference_UIDs_unchanged=True,ring_topology_reused=True,mean_q1_q05_q1_mean_L2=float(np.linalg.norm(qs[0]-f['q1'],axis=1).mean()),mean_q2_q05_q2_mean_L2=float(np.linalg.norm(qs[1]-f['q2'],axis=1).mean())));log('Mean source diagnostic complete')
