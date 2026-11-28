from common3 import *
from sklearn.metrics import silhouette_score,silhouette_samples,adjusted_rand_score as ari,adjusted_mutual_info_score
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from scipy.optimize import linear_sum_assignment
from scipy.signal import find_peaks
import warnings
rng=np.random.default_rng(seed);C=OUT/'classification'
families={'A':'amount','B':'supply_mixing','C':'supply_mixing','D':'heterogeneity','E':'organization','F':'amount','G':'amount','P':'periphery'}
def prep(a,train,cols=list('ABCDEFG')):
    tr=a.loc[train,cols].copy();full=a[cols].copy();chosen=[];aud=[]
    for k in cols:
        if k in ['B','C']:tr[k]=np.log1p(tr[k]);full[k]=np.log1p(full[k])
        vals=tr[k];rho=tr[chosen].corrwith(vals,method='spearman').abs() if chosen else pd.Series(dtype=float)
        reason='incomplete_optional' if vals.isna().any() else 'zero_variance' if vals.std(ddof=0)<1e-12 else 'redundant_'+rho.idxmax() if len(rho) and rho.max()>.8 else 'kept'
        if reason=='kept':chosen.append(k)
        aud.append(dict(feature=k,family=families[k],reason=reason,missing=vals.isna().sum(),zero=vals.eq(0).sum(),rho_max=rho.max(),mean=vals.mean(),sd=vals.std(ddof=0),transform='log1p' if k in ['B','C'] else 'identity'))
    if not chosen:raise ValueError('No complete nonconstant features')
    z=(full[chosen]-tr[chosen].mean())/tr[chosen].std(ddof=0)
    weight={k:sum(families[j]==families[k] for j in chosen)**-.5 for k in chosen}
    for v in aud:v['family_weight']=weight.get(v['feature'],np.nan)
    return z,z.mul(weight),aud
def pam(df,train,k):
    t=df.loc[train];ds=cdist(t,t);med=[int(ds.sum(0).argmin())];near=ds[:,med[0]]
    while len(med)<k:
        costs=np.minimum(ds,near[:,None]).sum(0);costs[med]=np.inf;m=int(np.argmin(costs));med.append(m);near=np.minimum(near,ds[:,m])
    for _ in range(200):
        old=ds[:,med].min(1).sum();best=old;swap=None
        for j in range(k):
            oth=np.delete(med,j);cost=np.minimum(ds,ds[:,oth].min(1)[:,None]).sum(0);cost[med]=np.inf;i=cost.argmin()
            if cost[i]<best-1e-9:best=cost[i];swap=(j,i)
        if swap is None:break
        med[swap[0]]=swap[1]
    # Deterministic names ordered by first retained feature at the medoid.
    med=sorted(med,key=lambda i:(t.iloc[i,0],str(t.index[i])));meds=list(t.index[med]);dist=cdist(df,t.iloc[med]);valid=np.isfinite(dist).all(axis=1);lab=pd.Series(-1,index=df.index);lab.loc[valid]=dist[valid].argmin(1)
    return lab,meds,dist
def mixture(z,train,k,cov='diag',reg=.001,rs=seed):
    x=z.loc[train].to_numpy();model=GaussianMixture(k,covariance_type=cov,reg_covar=reg,n_init=20,max_iter=500,random_state=rs).fit(x)
    if not model.converged_:raise ValueError('GMM did not converge')
    lab=pd.Series(-1,index=z.index);ok=z.notna().all(axis=1);lab.loc[ok]=model.predict(z.loc[ok].to_numpy());return model,lab
def align(ref,new,idx,k):
    v=ref.loc[idx];w=new.loc[idx];ok=(v>=0)&(w>=0);tab=np.zeros((k,k),int);np.add.at(tab,(v[ok].to_numpy(),w[ok].to_numpy()),1);rr,cc=linear_sum_assignment(-tab);return new.map(dict(zip(cc,rr))).fillna(-1).astype(int)
def redraw(idx):
    meta=M.loc[idx].copy();meta['id']=meta.unit.fillna('section:'+meta.index.to_series())
    tmp=meta.reset_index().groupby('id').agg(stratum=('dataset_id',lambda x:'|'.join(sorted(set(x)))),sections=('section_id',list))
    sampled=[]
    for _,g in tmp.groupby('stratum'):
        for u in rng.choice(g.index,len(g),replace=True):sampled.extend(tmp.at[u,'sections'])
    return sampled
def main():
    feats=read(O/'cache/immune_features.tsv');a=feats[(feats.grid==GRID)&(feats.axis==AXIS)].set_index('section_id')
    assert not any(any(k in v.lower() for k in ['mp10','p95','highlow']) for v in a.columns)
    tr=M.index[M.tumor & a.A.reindex(M.index).notna()];z,x,au=prep(a,tr)
    save(au,C/'feature_selection.tsv');save(z.reset_index(),C/'X_z.tsv');save(x.reset_index(),C/'X_pam.tsv')
    pc=PCA().fit(x.loc[tr]);codes=['PC'+str(j+1) for j in range(len(x.columns))];pcs=pd.DataFrame(pc.transform(x),index=x.index,columns=codes);save(pcs.reset_index(),C/'PCA_scores.tsv');save(pd.DataFrame(pc.components_.T,index=x.columns,columns=codes).rename_axis('feature').reset_index(),C/'PCA_loadings.tsv');save(dict(PC=codes,variance=pc.explained_variance_ratio_),C/'PCA_variance.tsv')
    refs={};pams=[];details=[];best={};gmods=[];gout=[];m3=None
    for k in [2,3,4]:
        labels,med,d=pam(x,tr,k);refs[('PAM',k)]=labels
        sil=silhouette_samples(x.loc[tr],labels.loc[tr]);sd=np.sort(d,axis=1);pams.append(dict(K=k,mean_silhouette=sil.mean(),cost=sd[x.index.isin(tr),0].sum(),medoids='|'.join(med)))
        for i,s in enumerate(x.index):details.append(dict(section_id=s,K=k,cluster='C'+str(labels[s]+1),medoid=med[labels[s]],is_medoid=s in med,nearest=sd[i,0],next_nearest=sd[i,1],margin=sd[i,1]-sd[i,0],silhouette=sil[tr.get_loc(s)] if s in tr else np.nan))
        if k==3:m3=med
    save(pams,C/'PAM_models.tsv');save(details,C/'PAM_sections.tsv')
    for cov in ['diag','tied']:
        for reg in [.001,.0001,.01]:
            for k in range(1,5):
                mod,lab=mixture(z,tr,k,cov,reg);gmods.append(dict(cov=cov,reg=reg,K=k,BIC=mod.bic(z.loc[tr].to_numpy()),converged=mod.converged_,n_iter=mod.n_iter_,n_init=20,min_cluster_size=lab.loc[tr].value_counts().min(),ARI_vs_PAM3=ari(refs['PAM',3].loc[tr],lab.loc[tr])))
                if (cov,reg)==('diag',.001):
                    refs['GMM',k]=lab
                    for s,pr in zip(z.index,mod.predict_proba(z.to_numpy())):gout.append(dict(section_id=s,K=k,cluster=lab[s]+1,probability=max(pr),posteriors=json.dumps(pr.tolist())))
    gm=pd.DataFrame(gmods);kbest=int(gm.query("cov=='diag' and reg==0.001").sort_values('BIC').iloc[0].K)
    save(gm,C/'GMM_models.tsv');save(gout,C/'GMM_sections.tsv')
    seedchecks=[]
    for k in sorted({3,kbest}):
        for rs in [1,22,103,509,2901]:
            _,lab=mixture(z,tr,k,rs=rs);seedchecks.append(dict(K=k,seed=rs,n_init=20,ARI=ari(refs['GMM',k].loc[tr],lab.loc[tr])))
    save(seedchecks,C/'GMM_seeds.tsv')
    uni=[]
    for key in z:
        v=z.loc[tr,key].to_numpy();xx=np.linspace(v.min()-.3,v.max()+.3,301)
        for k in [1,2,3]:
            mod,_=mixture(z[[key]],tr,k);uni.append(dict(feature=key,method='unit_variance_GMM',K=k,BIC=mod.bic(v[:,None])))
        for h in [.8,1,1.25]:
            kde=stats.gaussian_kde(v,bw_method=lambda z:z.scotts_factor()*h);yy=kde(xx);ii=find_peaks(yy,prominence=.06*max(yy))[0]
            uni.append(dict(feature=key,method='KDE',Scott_factor=h,n_peaks=len(ii),peaks_z=json.dumps(xx[ii].tolist())))
    save(uni,C/'univariate.tsv')
    boot=[];assignment=[];changes=[];adjrows=[]
    for rule,idx in [('primary',tr),('mapped_only',tr[M.loc[tr,'unit'].notna()])]:
        for j in range(100):
            ids=redraw(idx);zz,xx,au2=prep(a,ids)
            mods=[('PAM',kk) for kk in ([2,3,4] if rule=='primary' else [3])]
            if rule=='primary':mods += [('GMM',kk) for kk in sorted({3,kbest})]
            for method,kk in mods:
                try:
                    lab=pam(xx,ids,kk)[0] if method=='PAM' else mixture(zz,ids,kk,rs=seed+j)[1]
                    ref=refs[method,kk];al=align(ref,lab,idx,kk);ok=al.loc[idx]>=0
                    boot.append(dict(rule=rule,rep=j,method=method,K=kk,features=''.join(xx.columns),ARI=ari(ref.loc[idx][ok],al.loc[idx][ok]),valid=ok.sum(),status='ok'))
                    for s in a.index:assignment.append(dict(rule=rule,rep=j,method=method,K=kk,section_id=s,aligned=al[s],match=(al[s]==ref[s]) if al[s]>=0 else np.nan))
                except (ValueError,np.linalg.LinAlgError) as exc:boot.append(dict(rule=rule,rep=j,method=method,K=kk,status=str(exc)))
    save(boot,C/'bootstrap_runs.tsv.gz');aa=pd.DataFrame(assignment);save(aa,C/'bootstrap_labels.tsv.gz')
    ag=aa.groupby(['rule','method','K','section_id'],as_index=False).agg(stability=('match','mean'),n_effective=('match','count'));save(ag,C/'bootstrap_stability.tsv')
    ref=refs['PAM',3]
    def sensitivity(name,dat,train,cols=list('ABCDEFG')):
        zz,xx,info=prep(dat,train,cols);lab,md,_=pam(xx,train,3);al=align(ref,lab,train,3);allok=tr.intersection(al[al>=0].index)
        changes.append(dict(rule=name,n_train=len(train),n_projectable=len(allok),features=''.join(xx.columns),ARI=ari(ref.loc[allok],al.loc[allok]),n_changed=(ref.loc[allok]!=al.loc[allok]).sum(),medoids='|'.join(md)))
        for s in al.index:adjrows.append(dict(section_id=s,rule=name,cluster='C'+str(al[s]+1) if al[s]>=0 else None,in_training=s in train,changed=al[s]!=ref[s]))
    for grid,line in [('1000','T'),('x250','T'),('y250','T'),('xy250','T'),('500','I')]:
        sensitivity(grid+'_'+line,feats[(feats.grid==grid)&(feats.axis==line)].set_index('section_id'),tr)
    shared=tr[a.loc[tr,'P'].notna()];sensitivity('core_on_P_subset',a,shared);sensitivity('add_P',a,shared,list('ABCDEFGP'))
    # Explicitly add P as its own family on the identical core subset, even when correlated with B.
    zz,xx,info=prep(a,shared);pp=a.P;xx['P']=(pp-pp.loc[shared].mean())/pp.loc[shared].std(ddof=0)
    ll,mm,_=pam(xx,shared,3);al=align(ref,ll,shared,3)
    changes.append(dict(rule='core_plus_P_unpruned',n_train=len(shared),n_projectable=len(tr),features=''.join(xx.columns),ARI=ari(ref.loc[tr],al.loc[tr]),n_changed=(ref.loc[tr]!=al.loc[tr]).sum(),medoids='|'.join(mm)))
    for s in al.index:adjrows.append(dict(section_id=s,rule='core_plus_P_unpruned',cluster='C'+str(al[s]+1),in_training=s in shared,changed=al[s]!=ref[s]))
    sensitivity('verified_patients',a,tr[M.loc[tr,'patient_uid'].notna()])
    for ds in M.dataset_id.unique():sensitivity('omit_'+ds,a,tr[M.loc[tr,'dataset_id']!=ds])
    save(changes,C/'sensitivity_models.tsv');save(adjrows,C/'sensitivity_sections.tsv.gz')
    profile=a.join(ref.add(1).rename('cluster')).loc[tr];save(profile.groupby('cluster')[list('ABCDEFGP')+['E_area','win_n','UMI','nGene','islands','outer_coverage','whole_fraction','density_median']].median().reset_index(),C/'cluster_profiles.tsv')
    ds=M.loc[tr,'dataset_id'];cross=pd.crosstab(ds,ref.loc[tr]+1);save(cross.reset_index(),C/'cluster_source.tsv')
    qc=[]
    for col in ['E_area','valid_E_area','win_n','UMI','nGene','islands','median_island','whole_ab','mixed_E_share','outer_coverage']:
        val=profile[col];eta=sum(len(v)*(v.mean()-val.mean())**2 for _,v in val.groupby(profile.cluster))/sum((val-val.mean())**2)
        qc.append(dict(feature=col,eta_squared=eta))
    save(qc,C/'QC_confounding.tsv');old=read(V/'classification/section_classification.tsv')[['section_id','provisional_architecture_group','group_reason']].set_index('section_id')
    out=M.join(a[['E_area','valid_E_area','win_n','outer_coverage']]).join(old).join(ag.query("rule=='primary' and method=='PAM' and K==3").set_index('section_id')[['stability','n_effective']]);out['exploratory_cluster']=ref.map(lambda j:'C'+str(j+1));out['data_evaluable']=a.A.notna();out['medoid']=out.index.isin(m3);out['architecture_interpretation']='neutral exploratory cluster';out['interpretation_evidence']='Continuous summaries and local maps; cluster membership alone does not establish two coexisting local states.';save(out.reset_index(),C/'section_classification.tsv')
    save(pd.crosstab(out.loc[tr,'provisional_architecture_group'],out.loc[tr,'exploratory_cluster']).reset_index(),C/'old_to_new.tsv')
    print('PAM',pams,'GMM optimum',kbest,'source_AMI',adjusted_mutual_info_score(ds,ref.loc[tr]),flush=True)
if __name__=='__main__':main()

