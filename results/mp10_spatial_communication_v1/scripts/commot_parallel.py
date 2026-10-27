"""Scheduling-only adapter for official COMMOT four complete COT components.
No LR partition, no parameter or solver change. Fixed original weighted sum order.
"""
from commot_loader import cotmod
import multiprocessing as mp,numpy as np,time
_SHARED=None
_NAMES=['cot_sparse','cot_row_sparse','cot_col_sparse','cot_blk_sparse']
def _solve(k):
 a,kw=_SHARED
 return getattr(cotmod,_NAMES[k])(*a,**kw[k])
def cot_combine_parallel(S,D,A,M,cutoff,eps_p=.1,eps_mu=None,eps_nu=None,rho=10,weights=(.25,)*4,nitermax=10000,stopthr=1e-8,verbose=False):
 global _SHARED
 expand=lambda v: v if isinstance(v,tuple) else (v,)*4
 ep=expand(eps_p);em=ep if eps_mu is None else expand(eps_mu);en=ep if eps_nu is None else expand(eps_nu);rh=expand(rho)
 kw=[dict(eps_p=ep[k],eps_mu=em[k],eps_nu=en[k],rho=rh[k],nitermax=nitermax,stopthr=stopthr,verbose=False) for k in range(4)]
 _SHARED=((S,D,A,M,cutoff),kw)
 with mp.get_context('fork').Pool(4) as pool:parts=pool.map(_solve,range(4))
 P={}
 for i in range(A.shape[0]):
  for j in range(A.shape[1]):
   if not np.isinf(A[i,j]):P[(i,j)]=float(weights[0])*parts[0][(i,j)]+float(weights[1])*parts[1][(i,j)]+float(weights[2])*parts[2][(i,j)]+float(weights[3])*parts[3][(i,j)]
 _SHARED=None
 return P
