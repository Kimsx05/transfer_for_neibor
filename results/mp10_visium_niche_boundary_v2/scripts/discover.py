from lib import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ari,silhouette_score,adjusted_mutual_info_score
from scipy.signal import find_peaks
rng=np.random.default_rng(SEED)
x=read(O/'classification/immune_input.tsv.gz').set_index('section_id')
tum=M.index[M.source_context!='adjacent_normal_source'];t=x.loc[tum]
# Prespecified priority within biological families. Never read a score table.
candidates=['whole_T_density','penetration_T_200','CLQ_E_T_R1','Moran_T_R1','window500_E_T_fraction_IQR','Lee_E_T_R1','mixed_E_share','whole_I_fraction']
sel=[];dec=[]
for k in candidates:
    cor=t[sel].corrwith(t[k],method='spearman').abs() if sel else pd.Series(dtype=float)
    reason='missing>20%' if t[k].isna().mean()>.2 else 'redundant:'+cor.idxmax() if len(cor) and cor.max()>.85 else 'include'
    if reason=='include':sel.append(k)
    dec.append(dict(feature=k,decision=reason,coverage=t[k].notna().sum(),max_abs_r=cor.max()))
save(dec,'classification/feature_selection.tsv')
tx=x[sel].copy()
for k in ['whole_T_density','CLQ_E_T_R1','window500_E_T_fraction_IQR']:
    if k in tx:tx[k]=np.log1p(tx[k])
v=tx.loc[tum].dropna();scale=StandardScaler().fit(v)
pca=PCA().fit(scale.transform(v));pdim=min(5,np.searchsorted(pca.explained_variance_ratio_.cumsum(),.8)+1)
pc=pca.transform(scale.transform(v))[:,:pdim]
pcall=pd.DataFrame(index=x.index,columns=['PC'+str(i+1) for i in range(pdim)],dtype=float);valid=tx.notna().all(axis=1);pcall.loc[valid]=pca.transform(scale.transform(tx[valid]))[:,:pdim]
save(pcall.reset_index(),'classification/pca_coordinates.tsv');save(pd.DataFrame(pca.components_.T,index=sel,columns=['PC'+str(i+1) for i in range(len(sel))]).rename_axis('feature').reset_index(),'classification/pca_loadings.tsv')
save(dict(pc=np.arange(1,len(sel)+1),variance=pca.explained_variance_ratio_,retained=np.arange(len(sel))<pdim),'classification/pca_variance.tsv')
uid=M.loc[v.index,'analysis_unit_uid'];units=uid.dropna().unique()
def block_indices():return np.concatenate([np.flatnonzero(uid.eq(u).fillna(False)) for u in rng.choice(units,len(units),replace=True)])
def fit(y,k,seed=SEED,starts=4):return GaussianMixture(k,covariance_type='diag',reg_covar=.05,random_state=int(seed),n_init=starts,max_iter=300).fit(y)
results=[];predictions=[];reps=[];models={}
for k in range(1,6):
    mod=fit(pc,k,starts=20);models[k]=mod
    lab=mod.predict(pc);pr=mod.predict_proba(pc)
    init=[ari(lab,fit(pc,k,seed,1).predict(pc)) for seed in range(10)]
    indices=[block_indices() for j in range(80)]
    boot=[ari(lab[uid.notna()],fit(pc[idx],k,SEED+j,2).predict(pc[uid.notna()])) for j,idx in enumerate(indices)]
    km=KMeans(k,random_state=SEED,n_init=20).fit_predict(pc)
    for j,val in enumerate(boot):reps.append(dict(K=k,resample=j,ARI=val))
    ds=M.loc[v.index,'dataset_id']
    results.append(dict(K=k,BIC=mod.bic(pc),silhouette=silhouette_score(pc,lab) if k>1 and len(set(lab))>1 else np.nan,min_size=min(np.bincount(lab,minlength=k)),init_ARI=np.median(init),bootstrap_ARI=np.median(boot),bootstrap_q10=np.quantile(boot,.1),kmeans_ARI=ari(lab,km),mean_probability=pr.max(1).mean(),source_AMI=adjusted_mutual_info_score(ds,lab),**{'eta2_'+col:safe(sum(len(z)*(z.mean()-t.loc[v.index,col].mean())**2 for _,z in t.loc[v.index,col].groupby(lab)),sum((t.loc[v.index,col]-t.loc[v.index,col].mean())**2)) for col in ['observed_area','UMI_median','nGene_median']}))
    for sec,l,p in zip(v.index,lab,pr):predictions.append(dict(section_id=sec,K=k,cluster='C'+str(l+1),probability=p.max(),probabilities=json.dumps(p.tolist())))
model=pd.DataFrame(results);model['BIC_improvement_over_K1']=model.BIC.iloc[0]-model.BIC;save(model,'classification/model_selection.tsv');save(predictions,'classification/all_cluster_probabilities.tsv');save(reps,'classification/cluster_stability.tsv.gz')
k=int(model.loc[model.BIC.idxmin(),'K']);pred=pd.DataFrame(index=x.index);pred['data_driven_cluster']='not_evaluable';pred['cluster_probability']=np.nan
pr=models[k].predict_proba(pcall.loc[valid].to_numpy());pred.loc[valid,'data_driven_cluster']=['C'+str(a+1) for a in pr.argmax(axis=1)];pred.loc[valid,'cluster_probability']=pr.max(axis=1)
save(pred.reset_index(),'classification/natural_clusters.tsv')
def peaks(a,fac):
    kde=stats.gaussian_kde(a,bw_method=lambda z:z.scotts_factor()*fac)
    xx=np.linspace(np.min(a)-.1*np.ptp(a),np.max(a)+.1*np.ptp(a),240);h=kde(xx)
    ix=find_peaks(h,prominence=.06*h.max())[0]
    valleys=[xx[l+np.argmin(h[l:r+1])] for l,r in zip(ix[:-1],ix[1:])]
    return xx,h,ix,valleys
modes=[];density=[];univ=[]
for col in sel:
    val=tx.loc[tum,col].dropna()
    for j in range(1,5):
        mod=fit(val.to_numpy()[:,None],j);univ.append(dict(metric=col,K=j,BIC=mod.bic(val.to_numpy()[:,None])))
    for fac in [.8,1.,1.25]:
        xx,h,ix,valleys=peaks(val,fac)
        counts=[];troughs=[]
        for i in range(60):
            sample=tx.loc[v.index[block_indices()],col].dropna()
            _,_,ids,vs=peaks(sample,fac);counts.append(len(ids))
            if len(vs)==1:troughs.append(vs[0])
        modes.append(dict(metric=col,bandwidth_factor=fac,n_peaks=len(ix),peak_locations=json.dumps(xx[ix].tolist()),valleys=json.dumps(valleys),boot_same_n=np.mean(np.array(counts)==len(ix)),boot_two_peak_frequency=np.mean(np.array(counts)==2),valley_q10=np.quantile(troughs,.1) if troughs else np.nan,valley_q90=np.quantile(troughs,.9) if troughs else np.nan))
        density.extend(dict(metric=col,bandwidth_factor=fac,x=a,density=b) for a,b in zip(xx,h))
save(modes,'classification/multimodality.tsv');save(density,'plot_data/KDE.tsv.gz');save(univ,'classification/univariate_mixture.tsv')
save(t.corr(method='spearman').rename_axis('metric').reset_index(),'classification/correlations.tsv.gz')
save(t.isna().reset_index(),'classification/missingness.tsv.gz')
spec=dict(features=sel,mean=scale.mean_.tolist(),sd=scale.scale_.tolist(),log1p=[k for k in ['whole_T_density','CLQ_E_T_R1','window500_E_T_fraction_IQR'] if k in sel],n_PCs=int(pdim),K=k,seed=SEED,covariance='diagonal',reg_covar=.05,n_complete=len(v),selection='minimum BIC in single fixed PC space',no_outlier_deletion=True)
(O/'classification/natural_model.json').write_text(json.dumps(spec,indent=2));print(model.to_string(index=False))

