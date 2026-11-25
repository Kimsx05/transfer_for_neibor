import os,sys,json,runpy
from pathlib import Path
import numpy as np,pandas as pd,h5py,pyarrow as pa
from scipy import sparse
O=Path(__file__).resolve().parents[1]; P=O.parents[1]; A=P/'highres_VisiumHD/analysis'
F=A/'hd_treg_mp10_spatial_all_v1'
HD=['Cervilla','D1_Johnson','UC_SCC']
XE=['A1R2','A1R4','A1R5','A1R6','A2R2']
pa.set_cpu_count(1);pa.set_io_thread_count(1)
def put(t,path): t.to_csv(path,sep='\t',index=False,na_rep='NA')
def tenx(path):
    with h5py.File(path) as f:
        m=f['matrix']; g=m['features']
        x=sparse.csr_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:][::-1]));x.eliminate_zeros()
        return x,m['barcodes'].asstr()[:],pd.DataFrame({k:g[k].asstr()[:] for k in ['name','id','feature_type']})
def export(x,genes,ids,d):
    # CSC slots of genes x objects are the CSR slots of its transpose.
    x=x.tocsr();x.sort_indices()
    for key,v in [('i',x.indices),('p',x.indptr),('x',x.data)]:
        np.asarray(v,dtype='<i4').tofile(d/(key+'.bin'))
    pd.Series(genes).to_csv(d/'genes.txt',index=False,header=False)
    pd.Series(ids).to_csv(d/'ids.txt',index=False,header=False)
def read_export(d):
    genes=np.atleast_1d(np.loadtxt(d/'genes.txt',dtype=str));ids=np.atleast_1d(np.loadtxt(d/'ids.txt',dtype=str))
    x=sparse.csr_matrix(tuple(np.fromfile(d/(k+'.bin'),dtype='<i4') for k in ['x','i','p']),shape=(len(ids),len(genes)))
    return x,genes,ids

