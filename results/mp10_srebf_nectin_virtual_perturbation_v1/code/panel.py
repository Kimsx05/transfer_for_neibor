from common import *
from scipy import sparse as sp
X=sp.load_npz('cache/selected_counts.npz'); meta=rd('cache/metadata.tsv.gz');features=rd('cache/features.tsv').feature
M=X[meta.main.to_numpy()].astype(np.float64); norm=M.multiply((1e4/np.asarray(M.sum(1)).ravel())[:,None]).tocsr();norm.data=np.log1p(norm.data)
mean=np.asarray(norm.mean(0)).ravel(); var=np.asarray(norm.power(2).mean(0)).ravel()-mean**2
pos=np.asarray(M.getnnz(0)); pat=sum((np.asarray(M[meta.loc[meta.main,'patient_id'].eq(p).to_numpy()].getnnz(0))>0).astype(int) for p in meta.patient_id.unique())
d=pd.DataFrame(dict(feature=features,mean=mean,variance=var,n_positive=pos,n_positive_patients=pat))
a=rd('qc/feature_map.tsv').dropna(subset=['raw_feature']).set_index('raw_feature')
d['gene']=d.feature.map(a.gene).fillna(d.feature);ok=d[(d.n_positive>=10)&(pat>=3)&(var>0)].copy()
ok['dispersion']=ok.variance/ok['mean'];bins=pd.qcut(ok['mean'],20,duplicates='drop')
ok['hvg_z']=ok.groupby(bins,observed=True).dispersion.transform(lambda x:(x-x.mean())/x.std());ok=ok.sort_values(['hvg_z','gene'],ascending=[False,True])
ok['hvg_top2000']=np.arange(len(ok))<2000
pres=ok[ok.gene.isin(a.gene)];chosen=pd.concat([pres,ok[~ok.gene.isin(pres.gene)].head(2000-len(pres))]).sort_values('gene')
chosen['reserved']=chosen.gene.isin(pres.gene);save(chosen,'tables/model_panel.tsv');save(d,'qc/panel_candidates.tsv.gz')
assert len(chosen)==2000 and chosen.gene.is_unique
A=X[:,chosen.index].toarray().astype('<i4');assert A.shape==(len(meta),2000)
A.tofile('cache/cells_by_genes_int32.bin')
# R's column-major matrix reads one cell at a time into a genes x cells matrix.
save(pd.DataFrame(dict(barcode=meta.barcode,panel_UMI=A.sum(axis=1))),'qc/orientation.tsv.gz')

