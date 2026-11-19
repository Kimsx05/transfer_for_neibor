from base import *
import anndata as ad,joblib,sklearn
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from datetime import datetime
kinds=['T8_from2','T16_from2','cell_nucleus','cell_conservative','existing_cell8']
def thin(x,depth,rng):
 z=x.copy();prob=np.minimum(1,depth/np.maximum(1,np.asarray(z.sum(1)).ravel()));z.data=rng.binomial(z.data.astype(int),np.repeat(prob,np.diff(z.indptr)));z.eliminate_zeros();return z
def lognorm(x):
 lib=np.asarray(x.sum(1)).ravel();z=x.multiply(10000/np.maximum(1,lib)[:,None]).tocsr();z.data=np.log1p(z.data);return z
def features(z,y,sy):
 features=set();overall=np.asarray(z.mean(0)).ravel()
 for lab in np.unique(y):
  mu=np.asarray(z[y==lab].mean(0)).ravel()-overall
  valid=np.flatnonzero(~np.isin(sy,exclude));features.update(valid[np.argsort(mu[valid],kind='stable')[-24:]])
 return np.array(sorted(features))
def train(z,y,ind,C):
 return LogisticRegression(C=C,solver='lbfgs',max_iter=400,tol=.001,random_state=123).fit(z[:,ind],y)
def project(raw,fn,model,threshold,sy):
 z=lognorm(raw);p=model.predict_proba(z[:,fn]);order=np.argsort(p,axis=1);top=order[:,-1];second=order[:,-2];ix=np.arange(len(top));lab=model.classes_[top];score=p[:,list(model.classes_).index('Treg')];t=np.asarray(raw[:,np.isin(sy,T)].sum(1)).ravel();positive=(lab=='Treg')&(score>=threshold)&(t>0);zero=raw.getnnz(1)==0
 d=pd.DataFrame(dict(predicted=lab,runner_up=model.classes_[second],margin=p[ix,top]-p[ix,second],Treg_support=score,model_candidate=(lab=='Treg')&~zero,supported_tier=positive&~zero))
 d.loc[zero,['predicted','runner_up','margin','Treg_support']]=None
 return d
def stats(meta,d,depth):
 records=[];truth=meta.label.eq('Treg').to_numpy();pred=d.supported_tier.to_numpy();tp=np.sum(truth&pred);fp=np.sum(~truth&pred)
 for label in ['ALL']+list(meta.label.unique()):
  m=np.ones(len(meta),bool)if label=='ALL' else meta.label.eq(label).to_numpy()
  records.append(dict(depth=depth,truth=label,n=m.sum(),samples=meta['sample'][m].nunique(),positive=(truth&m).sum(),negative=(~truth&m).sum(),TP=(pred&truth&m).sum(),FP=(pred&~truth&m).sum(),FN=(~pred&truth&m).sum(),raw_Treg=int(d.model_candidate[m].sum()),rejected_raw_candidates=int((d.model_candidate[m]&~d.supported_tier[m]).sum()),no_call=d.predicted[m].isna().sum(),precision=tp/max(tp+fp,1)if label=='ALL'and tp+fp else np.nan,recall=(pred&truth&m).sum()/(truth&m).sum()if (truth&m).any()else np.nan,FPR=(pred&~truth&m).sum()/(~truth&m).sum()if (~truth&m).any()else np.nan))
 return records
def main():
 meta=pd.read_parquet(V/'tables/reference_membership.parquet');ref=ad.read_h5ad(js(V/'models/reference_prepared.json')['reference'],backed='r');assert ref.var_names.is_unique
 ii=ref.obs_names.get_indexer(meta.cell_id);assert min(ii)>=0;raw=sparse.csr_matrix(ref.X[ii,:]);assert np.array_equal(raw.data,np.round(raw.data));rs=ref.var_names.to_numpy();ref.file.close()
 common=set(rs)
 for s in sections:common&=set(rd(R/s/'genes.tsv').symbol)
 keep=np.flatnonzero(np.isin(rs,list(common)));sy=rs[keep];x=raw[:,keep].tocsr();z=lognorm(x);y=meta.label.to_numpy();tr=np.flatnonzero(meta.partition.eq('train'));test=np.flatnonzero(meta.partition.eq('sample_holdout'));assert set(meta['sample'][tr]).isdisjoint(meta['sample'][test])
 folds=StratifiedGroupKFold(n_splits=3,shuffle=True,random_state=123);cv=[];oof={};foldmap=[]
 for fold,(it,iv)in enumerate(folds.split(tr,y[tr],meta['sample'][tr])):
  it,iv=tr[it],tr[iv];ind=features(z[it],y[it],sy);low=thin(x[iv],100,np.random.default_rng(123+fold));foldmap.extend((meta.cell_id[i],fold)for i in iv)
  for C in [.1,1.]:
   model=train(z[it],y[it],ind,C);d=project(low,ind,model,0,sy);d['true']=y[iv];d['T_present']=np.asarray(low[:,np.isin(sy,T)].sum(1)).ravel()>0;oof.setdefault(C,[]).append(d);cv.append(dict(fold=fold,C=C,macro_recall=balanced_accuracy_score(y[iv],d.predicted.fillna('unresolved')),features=len(ind)))
 cv=pd.DataFrame(cv);put(cv,'reference_internal_cv.tsv');put(pd.DataFrame(foldmap,columns=['cell_id','internal_fold']),'training_folds.tsv.gz');C=cv.groupby('C').macro_recall.mean().idxmax();d=pd.concat(oof[C]);cut=np.inf;truth=d.true.eq('Treg')
 for threshold in np.round(np.arange(.3,1,.05),2):
  pred=d.model_candidate&d.T_present&d.Treg_support.ge(threshold);tp=(pred&truth).sum();fp=(pred&~truth).sum()
  if tp>=10 and tp/(tp+fp)>=.8 and fp/(~truth).sum()<=.01:cut=threshold;break
 ind=features(z[tr],y[tr],sy);model=train(z[tr],y[tr],ind,C);joblib.dump(dict(model=model,features=sy[ind],normalization_symbols=sy,threshold=cut,C=C),R/'model.joblib')
 cfg=yaml.safe_load((R/'method_config.yaml').read_text());cfg.update(selected_C=float(C),threshold_value=float(cut)if np.isfinite(cut)else None,feature_names=sy[ind].tolist(),sklearn=sklearn.__version__,model_freeze_time=datetime.now().astimezone().isoformat());(R/'method_config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
 valid=[];conf=[];rng=np.random.default_rng(123)
 for dep in [25,100,500,'native']:
  xx=x[test]if dep=='native' else thin(x[test],dep,rng);d=project(xx,ind,model,cut,sy);m=meta.iloc[test].reset_index(drop=True);valid.extend(stats(m,d,dep));c=pd.crosstab(m.label,d.predicted.fillna('unresolved')).stack().rename('n').reset_index();c['depth']=dep;conf.append(c)
 put(pd.DataFrame(valid),'reference_validation.tsv');put(pd.concat(conf),'confusion.tsv')
 # Paired with frozen v2 parent profiles; new count draws only for selected original mixture conditions.
 par=pd.read_parquet(F/'tables/mixture_parent_membership.parquet');old=pd.read_parquet(F/'tables/expanded_mixtures.parquet');profile=old[(old.partition=='sample_holdout')&(old.background=='Epi70')&(old.rep<10)&(old.depth==500)].drop_duplicates(['profile_id','Treg_T_RNA_fraction']);mixrows=[];mats=[];cellidx=pd.Index(meta.cell_id)
 for q in profile.itertuples():
  b=par[par.profile_id==q.profile_id];prob=np.zeros(len(rs))
  for role,weight in [('Epithelial',.7),(q.competitor,.3-q.Treg_RNA_fraction),('Treg',q.Treg_RNA_fraction)]:
   vals=np.asarray(raw[cellidx.get_indexer(b[b.role.eq(role)].cell_id)].sum(0)).ravel();prob+=max(0,weight)*vals/vals.sum()
  mats.append(sparse.csr_matrix(rng.multinomial(500,prob/prob.sum())));mixrows.append(dict(profile_id=q.profile_id,competitor=q.competitor,Treg_RNA_fraction=q.Treg_RNA_fraction,depth=500))
 if mats:put(pd.concat([pd.DataFrame(mixrows),project(sparse.vstack(mats)[:,keep],ind,model,cut,sy)],axis=1),'holdout_mixture_predictions.tsv')
 for s in sections:
  genes=rd(R/s/'genes.tsv').symbol;gix=pd.Index(sy).get_indexer(genes);i=np.flatnonzero(gix>=0);mapping=sparse.csr_matrix((np.ones(len(i),int),(i,gix[i])),shape=(len(genes),len(sy)))
  for kind in kinds:
   units=pd.read_parquet(R/s/(kind+'.parquet'));rawx=sparse.load_npz(R/s/(kind+'.npz'));d=project((rawx@mapping).tocsr(),ind,model,cut,sy)
   units=pd.concat([units.reset_index(drop=True),d],axis=1);units['tier']=np.select([units.nCount.eq(0),units.supported_tier,units.model_candidate,units.FOXP3_raw_screen],['unmeasured_or_zero','reference_supported_Treg','uncertain_model_Treg','FOXP3_screen_only'],default='other_prediction');units['low_information']=units.nCount.lt(25);units.to_parquet(R/s/(kind+'_evidence.parquet'),index=False)
if __name__=='__main__':main()
