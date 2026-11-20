from base6 import *
import anndata as ad
from scipy.optimize import minimize_scalar
from scipy.special import softmax
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold

def thin(x,depth,rng):
 z=x.copy();p=np.minimum(1,depth/np.maximum(1,np.asarray(z.sum(axis=1)).ravel()))
 z.data=rng.binomial(z.data.astype(int),np.repeat(p,np.diff(z.indptr)));z.eliminate_zeros();return z
def normalize(x):
 lib=np.asarray(x.sum(axis=1)).ravel();z=x.multiply(10000/np.maximum(1,lib)[:,None]).tocsr();z.data=np.log1p(z.data);return z
def reference():
 m=pd.read_parquet(V/'tables/reference_membership.parquet')
 donor=rd(A.parents[1]/'mp10_k32_reference_identifiability_v5/00_audit/provenance/080_sample_patient_mapping.tsv').drop_duplicates('sample').set_index('sample').patient_recorded
 m['group_id']=(m.dataset+'|'+m['sample']).map(donor)
 m['group_id']=m.group_id.mask(m.group_id.fillna('').str.lower().isin(['','unknown','na'])).fillna(m.dataset+'::sample::'+m['sample'])
 h=ad.read_h5ad(js(V/'models/reference_prepared.json')['reference'],backed='r')
 idx=h.obs_names.get_indexer(m.cell_id);assert min(idx)>=0
 sy=h.var_names.to_numpy();common=set.intersection(*[set(rd(V5/s/'genes.tsv').symbol) for s in sections]);ix=np.flatnonzero(np.isin(sy,list(common)));x=sparse.csr_matrix(h.X[idx,:])[:,ix].tocsr();h.file.close()
 assert np.all(x.data==x.data.astype(int));return m,x,sy[ix]
def broad(m):
 return np.select([m.label.eq('Epithelial'),m.label.eq('NK'),m.major_celltype.eq('T/NK'),m.major_celltype.eq('Fibroblast')],['Epithelial','NK','T','Stromal'],default=m.major_celltype.to_numpy(dtype='str'))
def infer(model,z):
 est,features,temp=model
 return softmax(est.decision_function(z[:,features])/temp,axis=1)
def fit_stage(z,zlo,y,train,cal,sy):
 overall=np.asarray(z[train].mean(0)).ravel();allowed=np.flatnonzero(~np.isin(sy,top+'TIGIT NECTIN2 CD96 NECTIN1 PVR CD226'.split()));feat=set()
 for cl in np.unique(y[train]):
  delta=np.asarray(z[train[y[train]==cl]].mean(0)).ravel()-overall;feat.update(allowed[np.argsort(delta[allowed],kind='stable')[-45:]])
 f=np.array(sorted(feat));model=LogisticRegression(C=.2,solver='lbfgs',max_iter=450,random_state=123)
 model.fit(sparse.vstack([z[train][:,f],zlo[train][:,f]]),np.tile(y[train],2))
 logits=model.decision_function(zlo[cal][:,f]);truth=pd.Index(model.classes_).get_indexer(y[cal]);assert truth.min()>=0
 res=minimize_scalar(lambda t:-np.log(softmax(logits/t,axis=1)[np.arange(len(truth)),truth].clip(1e-10)).mean(),bounds=(.5,3),method='bounded');return model,f,res.x
def predict(model,x):
 sy,ma,sub=model
 cols=[]
 for a in range(0,x.shape[0],8192):
  z=normalize(x[a:a+8192]);p=infer(ma,z);q=infer(sub,z);labs=sub[0].classes_;labels=ma[0].classes_;j=np.flatnonzero(labs!='Treg');other=j[q[:,j].argmax(1)]
  ix=np.arange(z.shape[0]);cols.append(pd.DataFrame(dict(major=labels[p.argmax(1)],major_support=p.max(1),p_T=p[:,list(labels).index('T')],p_epi=p[:,list(labels).index('Epithelial')],p_Treg_given_T=q[:,list(labs).index('Treg')],T_subtype=labs[q.argmax(1)],nonTreg_type=labs[other],nonTreg_support=q[ix,other])))
 return pd.concat(cols,ignore_index=True)
def to_ref(x,symbols,sy):
 i=pd.Index(sy).get_indexer(symbols);j=np.flatnonzero(i>=0);return (x@sparse.csr_matrix((np.ones(len(j),np.int32),(j,i[j])),shape=(len(symbols),len(sy)))).tocsr()
def train():
 m,x,sy=reference();y=broad(m);rng=np.random.default_rng(123)
 spl=StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=123)
 folds=np.full(len(m),-1)
 for k,(_,i) in enumerate(spl.split(np.zeros(len(m)),m.label,m.group_id)):folds[i]=k
 m['role']=np.select([folds==0,folds==1],['evaluation','calibration'],default='training');m.to_parquet(R/'reference_split.parquet',index=False)
 lo=normalize(thin(x,rng.choice([30,100,300,600],size=len(m)),rng));z=normalize(x);mods=[]
 for name,labels,sel in [('major',y,np.ones(len(m),bool)),('within_T',m.label.to_numpy(),y=='T')]:
  mod=fit_stage(z,lo,labels,np.flatnonzero((folds>=2)&sel),np.flatnonzero((folds==1)&sel),sy);mods.append(mod);pd.Series(sy[mod[1]],name='symbol').to_csv(R/'tables'/(name+'_features.tsv'),index=False)
  print(name,len(mod[1]),mod[2],mod[0].classes_,flush=True)
 joblib.dump((sy,*mods),R/'model.joblib')
if __name__=='__main__':train()
