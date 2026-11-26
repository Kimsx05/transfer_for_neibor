from common import *
from scipy import sparse
import h5py
T=O/'.work'
def strings(g):
 if isinstance(g,h5py.Dataset): return g.asstr()[:]
 if 'values' in g:
  assert not np.any(g['mask'][:])
  return g['values'].asstr()[:]
 cats=g['categories'].asstr()[:]; codes=g['codes'][:]; assert min(codes)>=0
 return cats[codes]
def pack(name,m,genes,ids):
 # CSR cells x genes is CSC genes x cells for R
 m=m.tocsr();m.sum_duplicates();m.eliminate_zeros();m.sort_indices()
 assert len(ids)==m.shape[0] and len(genes)==m.shape[1] and pd.Index(ids).is_unique and pd.Index(genes).is_unique
 assert (m.data>=0).all() and np.array_equal(m.data,m.data.astype(np.int32))
 p=T/name;p.mkdir(exist_ok=True)
 for k,x in [('i',m.indices),('p',m.indptr),('x',m.data)]:x.astype('<i4').tofile(p/k)
 for k,x in [('genes',genes),('ids',ids)]: (p/k).write_text('\n'.join(x)+'\n')
def slice_h5(f,idx):
 # Read small consecutive blocks, in physical order.
 x=f['X'];ptr=x['indptr'][:]; idx=np.asarray(idx)
 out=[];groups=np.split(np.sort(idx),np.where(np.diff(np.sort(idx))>24)[0]+1)
 for g in groups:
  a,b=ptr[g[0]],ptr[g[-1]+1]
  block=sparse.csr_matrix((x['data'][a:b],x['indices'][a:b],ptr[g[0]:g[-1]+2]-a),shape=(g[-1]-g[0]+1,x.attrs['shape'][1]))
  out.append(block[g-g[0]])
 return sparse.vstack(out,format='csr')[np.argsort(np.argsort(idx))]
