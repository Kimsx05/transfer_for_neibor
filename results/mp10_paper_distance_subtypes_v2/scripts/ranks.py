from run import *
from scipy.io import mmread
from scipy import sparse
import gzip
def rank(sec):
 r=M.set_index('section_id').loc[sec];d=pd.read_parquet(F/'spots'/f'{sec}.parquet').set_index('barcode')
 bc=Path(r.barcodes_path).read_text().splitlines();genes=Path(r.genes_path).read_text().splitlines();d=d.loc[bc];idx=np.where(d['query'])[0]
 if len(idx)<20:return sec
 with gzip.open(r.counts_path,'rb') if str(r.counts_path).endswith('.gz') else open(r.counts_path,'rb') as f:x=mmread(f).tocsr()[:,idx].astype(float)
 sums=np.asarray(x.sum(0)).ravel();x=x@sparse.diags(1e4/np.maximum(sums,1));x.data=np.log1p(x.data)
 g=d.iloc[idx];C=np.column_stack([np.ones(len(idx)),g.Epithelial_fraction,np.log1p(g.nCount_RNA),pd.get_dummies(g.lesion_location,drop_first=True,dtype=float)])
 y=g[RAW].to_numpy();res=y-C@np.linalg.lstsq(C,y,rcond=None)[0]
 mu=np.asarray(x.mean(1)).ravel();var=np.asarray(x.power(2).mean(1)).ravel()-mu**2
 score=np.asarray(x@res).ravel()/np.sqrt(np.maximum(var,1e-12)*len(idx)*np.sum(res**2))
 write(pd.DataFrame(dict(gene=genes,stat=score)).groupby('gene',as_index=False).stat.mean(),O/'tables'/f'{sec}_gene_ranks.tsv.gz')
 return sec
if __name__=='__main__':
 with ProcessPoolExecutor(max_workers=6) as p:
  for sec in p.map(rank,M.section_id):print(sec,flush=True)
