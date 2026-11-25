from base import *
ref=O/'shared/reference';x,genes,ids=read_export(ref);m=pd.read_csv(ref/'membership.tsv.gz',sep='\t').set_index('barcode').loc[ids];m=m[m.validate].copy();ix=pd.Index(ids).get_indexer(m.index);x=x[ix].copy()
rng=np.random.default_rng(123);depths={}
for s in XE:
 d=O/'xenium'/s;z=pd.read_csv(d/'identity.tsv.gz',sep='\t').set_index('id');c=pd.read_parquet(d/'cells.parquet').loc[z.index];depths[s]=c.loc[z.epi_eligible,'nGeneTranscripts'].to_numpy()
m['target_sample']=rng.choice(XE,len(m));m['target_transcripts']=[rng.choice(depths[s]) for s in m.target_sample];m['original_panel_depth']=np.asarray(x.sum(1)).ravel()
for i,n in enumerate(m.target_transcripts):
 a,b=x.indptr[i:i+2];counts=x.data[a:b].astype('int64')
 if counts.sum()>n:x.data[a:b]=rng.multivariate_hypergeometric(counts,int(n),method='marginals')
x.eliminate_zeros();m['thinned_panel_depth']=np.asarray(x.sum(1)).ravel();d=O/'shared/validation';d.mkdir(exist_ok=True);export(x,genes,m.index,d);put(m.reset_index(),d/'depth_design.tsv.gz')

