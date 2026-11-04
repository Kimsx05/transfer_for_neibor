from core import *
import h5py
from scipy import io
# Reuse precisely the previous 19,053-cell stratified source-support sample.
meta=pd.read_csv(OLD/'tables/sc_cells.tsv',sep='\t');p=R/'subtype/v2/reference_training_v2/04_gene_harmonization/reference_v2_raw_counts_minimally_harmonized.h5ad'
with h5py.File(p) as f:
 g=f['obs']['_index'];cells=g['values'].asstr()[:] if isinstance(g,h5py.Group) else g.asstr()[:];lookup=pd.Series(np.arange(len(cells)),index=cells);sel=lookup.loc[meta.cell].to_numpy();x=f['X'];ptr=x['indptr'][:];rows=[]
 for i in sel:
  a,b=ptr[i:i+2];rows.append(sparse.csr_matrix((x['data'][a:b],x['indices'][a:b],[0,b-a]),shape=(1,x.attrs['shape'][1])))
 mat=sparse.vstack(rows).T;genes=f['var']['_index']['values'].asstr()[:];io.mmwrite(O/'scores/sc_temporary.mtx',mat);pd.Series(genes).to_csv(O/'scores/sc_genes.txt',index=False,header=False);save(meta,O/'tables/sc_cells.tsv')
