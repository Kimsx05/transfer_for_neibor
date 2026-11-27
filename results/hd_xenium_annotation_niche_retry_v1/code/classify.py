from base_data import *
import celltypist
from scipy.special import expit
from sklearn.metrics import precision_recall_fscore_support,confusion_matrix
def score(model,x,genes):
 ix=pd.Index(genes).get_indexer(model.classifier.features)
 if np.any(ix<0):raise ValueError('unmapped model features')
 parts=[]
 for i in range(0,x.shape[0],2048):
  q=model.scaler.transform(x[i:i+2048,ix]);q.data=np.minimum(q.data,10);d=model.classifier.decision_function(q)
  if d.ndim==1:d=np.c_[-d,d]
  parts.append(d)
 return np.vstack(parts),model.classifier.classes_
def top(d,classes):
 ix=np.argsort(-d,axis=1,kind='stable');r=np.arange(len(ix));return pd.DataFrame(dict(first=classes[ix[:,0]],second=classes[ix[:,1]],prob=expit(d[r,ix[:,0]]).astype(float),margin=(d[r,ix[:,0]]-d[r,ix[:,1]]).astype(float)))
def xe():
 rng=np.random.default_rng(1729);x,g,ids=read_export(V/'shared/reference');ref=pd.read_csv(V/'shared/reference/membership.tsv.gz',sep='\t').set_index('barcode').loc[ids]
 identity_frame=ref.identity_ref.to_numpy(bool);x=x[identity_frame];ref=ref.loc[identity_frame].copy()
 ref['unit']=ref.patient_uid_recorded.fillna(ref.sample_uid);held=[]
 for ds,v in ref.groupby('dataset'):
  units=np.sort(v.unit.unique());held.extend(rng.choice(units,max(1,len(units)//5),replace=False))
 train=~ref.unit.isin(held);assert set(ref.unit[train]).isdisjoint(ref.unit[~train])
 cols=['unit','sample_uid','dataset','analysis_celltype_c2l_final','analysis_major_c2l_final'];ref[cols].assign(split=np.where(train,'train','held')).to_parquet(O/'xenium/reference_members.parquet')
 x=norm(x);metrics=[];con=[];models={};lev=ref.analysis_major_c2l_final.values;sub=ref.analysis_celltype_c2l_final.values
 for head in ['major','T/NK','Myeloid','Fibroblast','Endothelial','B']:
  ok=np.ones(len(ref),bool) if head=='major' else lev==head;y=lev if head=='major' else sub
  tr=ok&train.to_numpy();te=ok&~train.to_numpy()
  if len(np.unique(y[tr]))<2:continue
  o=O/'xenium'/('ct_identity_'+head.replace('/','_')+'.pkl')
  m=celltypist.models.Model.load(str(o)) if o.exists() else celltypist.train(X=x[tr],labels=y[tr],genes=g,with_mean=False,use_SGD=True,max_iter=150,n_jobs=1,check_expression=False,random_state=1729);m.write(str(o));models[head]=m
  sc,cl=score(m,x[te],g);q=top(sc,cl);pr,rec,_,num=precision_recall_fscore_support(y[te],q['first'],labels=cl,zero_division=0)
  metrics.extend(zip([head]*len(cl),cl,pr,rec,num));cm=confusion_matrix(y[te],q['first'],labels=cl);con.extend((head,c,b,cm[i,j]) for i,c in enumerate(cl) for j,b in enumerate(cl))
 write(pd.DataFrame(metrics,columns=['head','class','precision','recall','held_n']),'reference_evaluation.tsv');write(pd.DataFrame(con,columns=['head','truth','prediction','n']),'reference_confusion.tsv')
 for s in XE:
  a=dataset(s);ix=pd.Index(a.var_names).get_indexer(g);assert (ix>=0).all();xx=norm(a.X[:,ix]);t=top(*score(models['major'],xx,g)).add_prefix('ct_major_');res=t.copy()
  for head,m in models.items():
   if head=='major':continue
   ii=np.flatnonzero(t.ct_major_first==head)
   if not len(ii):continue
   add=top(*score(m,xx[ii],g)).add_prefix('ct_sub_');add.index=ii
   for col in add:res.loc[ii,col]=add[col]
  res.index=a.obs_names;paths(s).mkdir(exist_ok=True);res.to_parquet(paths(s)/'prediction.parquet')
def hd():
 for s in HD:
  if (paths(s)/'prediction.parquet').exists():continue
  a=dataset(s);mp=pd.read_csv(H/f'metadata/gene_mapping_{s}.tsv',sep='\t');gs=pd.Index(json.loads((H/'metadata/classifier_config.json').read_text())['gene_list']);ok=mp.reference_symbol.notna();ag=sparse.csr_matrix((np.ones(ok.sum(),'f4'),(mp.hd_feature_index[ok],gs.get_indexer(mp.reference_symbol[ok]))),shape=(len(mp),len(gs)));xx=norm(a.X@ag)
  m=celltypist.models.Model.load(str(H/'models/celltypist_major.pkl'));res=top(*score(m,xx,gs)).add_prefix('ct_major_');m=celltypist.models.Model.load(str(H/'models/celltypist_subtype.pkl'));d,cl=score(m,xx,gs)
  orig=pd.read_csv(V/'shared/taxonomy.tsv',sep='\t').set_index('original');
  for maj in res.ct_major_first.unique():
   ii=np.flatnonzero(res.ct_major_first==maj);jj=np.flatnonzero([c in orig.index and (orig.at[c,'major']==maj or (maj=='Fibroblast' and orig.at[c,'major']=='Mural')) for c in cl])
   if len(jj)<2:continue
   q=top(d[np.ix_(ii,jj)],cl[jj]).add_prefix('ct_sub_')
   for col in q:res.loc[ii,col]=q[col].values
   if maj=='Fibroblast':
    mural=q.ct_sub_first.map(orig.major).eq('Mural').values;res.loc[ii[mural],'ct_major_first']='Mural'
  res['ct_sub_first_original']=res.ct_sub_first;res['ct_sub_first']=res.ct_sub_first.map(orig.label)
  res.index=a.obs_names;res.to_parquet(paths(s)/'prediction.parquet')
if __name__=='__main__':(xe if sys.argv[1]=='xenium' else hd)()

