from commot_loader import *
from commot_parallel import cot_combine_parallel
import pandas as pd,numpy as np,anndata as ad
from scipy import sparse
rng=np.random.default_rng(173);X=rng.integers(0,8,size=(12,6)).astype(float);X[:2,2:]=0;X[2:4,:2]=0
v=['L1','L2','R1','R2','R3','R4'];df=pd.DataFrame([['L1','R1','p'],['L1','R2_R3','q'],['L2','R1','p'],['L2','R4','r']]);xy=np.array([(100*(i%4),100*(i//4)) for i in range(12)],float)
a=ad.AnnData(sparse.csr_matrix(X),var=pd.DataFrame(index=v));a.obsm['spatial']=xy;b=a.copy();dis={(l,r):200.000001 for l in ['L1','L2'] for r in ['R1','R2_R3','R4']};dis[('L2','R4')]=100.000001
kwargs=dict(database_name='parity',df_ligrec=df,heteromeric=True,heteromeric_rule='min',dis_thr=dis,cot_nitermax=10000)
spatial.spatial_communication(a,**kwargs);spatial.cot_combine_sparse=cot_combine_parallel;spatial.spatial_communication(b,**kwargs)
rows=[]
for k in a.obsp:
 d=a.obsp[k]-b.obsp[k];err=float(np.max(np.abs(d.data))) if d.nnz else 0.;assert err<1e-12
 rows.append(dict(matrix=k,max_abs_error=err,status='PASS'))
pd.DataFrame(rows).to_csv(O/'commot/parallel_scheduling_parity.tsv',sep='\t',index=False);print('PARITY PASS',rows)
