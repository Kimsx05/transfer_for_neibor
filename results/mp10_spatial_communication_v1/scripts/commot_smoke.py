from commot_loader import *
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
X=np.array([[1,0,0,0],[2,0,0,0],[0,1,1,4],[0,3,4,1]],float)
a=ad.AnnData(sparse.csr_matrix(X),var=pd.DataFrame(index=['LigA','RecA','RecB','RecC']))
a.obsm['spatial']=np.array([[0,0],[1,0],[0,1],[1,2]],float)
df=pd.DataFrame([['LigA','RecA_RecB','PathA'],['LigA','RecC','PathA']])
dis={('LigA','RecA_RecB'):1.1,('LigA','RecC'):1.1}
spatial.spatial_communication(a,database_name='toy',df_ligrec=df,heteromeric=True,heteromeric_rule='min',dis_thr=dis)
T=a.obsp['commot-toy-total-total'].tocsr();assert np.isfinite(T.data).all();assert np.allclose(T[0,2],1.5,atol=1e-5)
Q=np.array([[.3,.7],[.8,.2],[.1,.9],[.5,.5]])
T.setdiag(0);T.eliminate_zeros();F=Q.T@T@Q;D=np.outer(Q.sum(0),Q.sum(0))-Q.T@Q
Fb=sum(T[i,j]*np.outer(Q[i],Q[j]) for i in range(4) for j in range(4) if i!=j)
Db=sum(np.outer(Q[i],Q[j]) for i in range(4) for j in range(4) if i!=j)
assert np.allclose(F,Fb) and np.allclose(D,Db)
db=database.ligand_receptor_database(species='human',signaling_type=None)
print('PASS official toy, full distance dictionary, cross-position F/D formula',T.sum(),len(db),db.iloc[:,3].value_counts().to_dict())
