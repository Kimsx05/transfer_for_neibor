from run import *
import h5py
from scipy import sparse,io
p=ROOT/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad'
with h5py.File(p) as f:
 def obs(k):
  g=f['obs'][k]
  if isinstance(g,h5py.Dataset):return g.asstr()[:]
  if "values" in g:return g["values"].asstr()[:]
  codes=g['codes'][:];cats=g['categories'].asstr()[:];return np.where(codes>=0,cats[np.maximum(codes,0)],None)
 meta=pd.DataFrame({'cell':obs('_index'),'sample':obs('orig.ident'),'type':obs('reference_celltype_v2'),'dataset':obs('dataset')})
 sel=[];rng=np.random.default_rng(20261103)
 coverage=[]
 for (sample,t),g in meta[meta.type.isin(list(TARGETS)[:7]+['Epithelial'])].groupby(['sample','type']):
  idx=rng.choice(g.index,size=min(100,len(g)),replace=False);sel.extend(idx)
  coverage.append(dict(sample=sample,type=t,n_available=len(g),n_scored=len(idx)))
 sel=np.array(sorted(sel));x=f['X'];ptr=x['indptr'][:];blocks=[]
 for i in sel:
  a,b=ptr[i:i+2];blocks.append(sparse.csr_matrix((x['data'][a:b],x['indices'][a:b],np.array([0,b-a])),shape=(1,x.attrs['shape'][1])))
 mat=sparse.vstack(blocks).T.tocoo();genes=f['var']['_index']['values'].asstr()[:]
 io.mmwrite(O/'sc_counts.mtx',mat);pd.Series(genes).to_csv(O/'sc_genes.txt',index=False,header=False);write(meta.loc[sel],O/'tables/sc_cells.tsv');write(coverage,O/'tables/sc_coverage.tsv')
