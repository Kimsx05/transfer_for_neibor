from common3 import *
from scipy.linalg import block_diag
T=OUT/'tables';C=OUT/'classification'
def locked():
    label=C/'section_classification.tsv';lock=json.loads((OUT/'classification_lock.json').read_text())
    assert sha(label)==lock['hashes'][str(label)]
    return read(label).rename(columns={'exploratory_cluster':'group'})
def draws(units,ds,n=500):
    u=pd.DataFrame(dict(id=units,ds=ds)).drop_duplicates();ids=sorted(u.id.unique());strata=u.groupby('id').ds.agg(lambda x:'|'.join(sorted(set(x))))
    rng=np.random.default_rng(seed+17);weights=np.zeros((n+1,len(ids)));weights[0]=1
    for stratum,aa in strata.groupby(strata):
        ii=np.array([ids.index(j) for j in aa.index]);weights[1:,ii]=rng.multinomial(len(ii),[1/len(ii)]*len(ii),size=n)
    return pd.Index(ids),weights
def globtest(coefs):
    d=coefs[0];bs=coefs[1:];ok=np.isfinite(bs).all(axis=1);bs=bs[ok]
    if len(bs)<400:return np.nan,'fewer_than_400_estimable_draws'
    cov=np.atleast_2d(np.cov(bs.T))
    if np.linalg.matrix_rank(cov)<len(d):return np.nan,'rank_deficient_bootstrap_covariance'
    inv=np.linalg.inv(cov);stat=d@inv@d;dev=bs-d
    return (1+sum(np.einsum('bi,ij,bj->b',dev,inv,dev)>=stat))/(len(bs)+1),'ok'
def unit_bins(bins,labels):
    a=bins.merge(labels[['section_id','dataset_id','unit','group','stability','tumor']],on='section_id')
    ok=a[a.tumor&a.unit.notna()&a.n.ge(5)]
    keys=['boundary','domain','center','dataset_id','unit','group']
    u=ok.groupby(keys,as_index=False).agg(raw=('median','mean'),centered=('centered','mean'),E_fraction=('E_fraction','mean'),logUMI=('logUMI','mean'),points=('n','sum'),sections=('section_id','nunique'))
    return a,u
def runmodel(u,orig,b,domain,gs,variant,support,tests,curves):
    key=dict(boundary=b,domain=domain,contrast=':'.join(gs),variant=variant)
    avail=orig[orig.tumor&orig.unit.notna()&orig.group.isin(gs)][['unit','dataset_id','group']].drop_duplicates()
    sources=sorted(x for x,g in avail.groupby('dataset_id') if set(g.group)==set(gs))
    avail=avail[avail.dataset_id.isin(sources)];uu=u[u.dataset_id.isin(sources)&u.group.isin(gs)]
    denom={g:avail[avail.group==g].unit.nunique() for g in gs}
    support_bins=[];cont=True;prev=0
    for x in sorted(uu.center.unique(),key=abs):
        p=uu[uu.center==x];en=True;sg=[]
        for g in gs:
            z=p[p.group==g];n=z.unit.nunique();pct=divide(n,denom[g]);en=en and pct>=.7 and set(z.dataset_id)==set(sources)
            sg.append(dict(key,group=g,center=x,units=n,den=denom[g],fraction=pct,sources=z.dataset_id.nunique(),n_sections=z.sections.sum(),points=z.points.sum()))
        cont=cont and en and abs(x)-prev<=100;prev=abs(x)
        if cont:support_bins.append(x)
        support.extend(dict(v,connected_common=bool(cont),both_or_all_common_sources='|'.join(sources)) for v in sg)
    d=uu[uu.center.isin(support_bins)].copy()
    if variant=='fixed_units' and len(d):
        uq=d.groupby(['unit','group']).center.nunique();ii=uq[uq==len(support_bins)].index
        d=d[d.set_index(['unit','group']).index.isin(ii)]
        s2=[s for s,v in d.groupby('dataset_id') if set(v.group)==set(gs)];d=d[d.dataset_id.isin(s2)];sources=s2
    ns={g:d[d.group==g].unit.nunique() for g in gs}
    info=dict(key,n_sources=len(sources),sources='|'.join(sources),n_units=d.unit.nunique(),n_by_group=json.dumps(ns),n_bins=len(support_bins),xmin=min(support_bins,default=np.nan),xmax=max(support_bins,default=np.nan))
    reason='no_common_sources' if not sources else 'fewer_than_three_contiguous_bins' if len(support_bins)<3 else 'fewer_than_two_common_sources' if len(sources)<2 else 'fewer_than_five_mapped_units_in_arm' if min(ns.values())<5 else 'ok'
    def entries(why):
        for scale in ['raw','centered']:
            for term in ['function','shape']:tests.append(dict(info,scale=scale,term=term,p=np.nan,status=why,n_boot=0))
    if reason!='ok':entries(reason);return
    xs=np.array(sorted(support_bins));mx=xs.mean();sc=max(abs(xs-mx))
    def basis(x):t=(x-mx)/sc;return np.c_[np.ones(len(t)),t,t*t]
    bx=basis(d.center.to_numpy());X=np.hstack([bx*d.group.eq(g).to_numpy()[:,None] for g in gs]+[d.dataset_id.eq(s).to_numpy()[:,None] for s in sources[1:]])
    if variant=='adjusted':
        adj=d[['E_fraction','logUMI']].to_numpy();X=np.c_[X,(adj-adj.mean(0))/adj.std(0)]
    yy=d[['raw','centered']].to_numpy();cell=pd.factorize(pd.MultiIndex.from_frame(d[['center','group','dataset_id']]))[0]
    ids,ww=draws(d.unit,d.dataset_id);iunit=ids.get_indexer(d.unit)
    beta=np.full((len(ww),X.shape[1],2),np.nan);fail={}
    for j,counts in enumerate(ww):
        w=counts[iunit];den=np.bincount(cell,weights=w);wt=w/np.maximum(den[cell],1)
        # Equal total weight per occupied source × group × distance cell.
        xtw=X.T*wt;gram=xtw@X
        if np.linalg.matrix_rank(gram)<len(gram):fail['rank_deficient']=fail.get('rank_deficient',0)+1;continue
        beta[j]=np.linalg.solve(gram,xtw@yy)
    if not np.isfinite(beta[0]).all():entries('rank_deficient_original');return
    good=np.isfinite(beta[1:,0,0]);info.update(n_boot=int(good.sum()),bootstrap_failures=json.dumps(fail),min_source_arm_units=int(d.groupby(['dataset_id','group']).unit.nunique().min()))
    for e,scale in enumerate(['raw','centered']):
        cs=[beta[:,3*i:3*i+3,e] for i in range(len(gs))]
        delta=np.concatenate([c-cs[0] for c in cs[1:]],axis=1);idx=np.arange(delta.shape[1])%3!=0
        for term,dd in [('function',delta),('shape',delta[:,idx])]:
            pv,st=globtest(dd);tests.append(dict(info,scale=scale,term=term,p=pv,status=st))
        for (gi,g),(hi,h) in itertools.combinations(enumerate(gs),2):
            pr=(cs[gi]-cs[hi])@basis(xs).T;lo,up=np.nanquantile(pr[1:],[.025,.975],axis=0)
            for x,fit,l,r in zip(xs,pr[0],lo,up):curves.append(dict(info,scale=scale,g1=g,g2=h,center=x,diff=fit,lo=l,hi=r))
def descriptive(a,u):
    rows=[]
    for keys,z in u.groupby(['boundary','domain','center','group']):
        ds=z.groupby('dataset_id')[['raw','centered']].mean()
        for ep in ['raw','centered']:
            ci=stats.t.ppf(.975,len(ds)-1)*ds[ep].sem() if len(ds)>1 else np.nan;mean=ds[ep].mean()
            rows.append(dict(zip(['boundary','domain','center','group'],keys),scale=ep,mean=mean,lo=mean-ci,hi=mean+ci,n_units=z.unit.nunique(),n_sections=z.sections.sum(),n_sources=len(ds),points=z.points.sum(),sources='|'.join(ds.index)))
    save(rows,T/'descriptive_group_curves.tsv.gz');save(u,T/'unit_curves.tsv.gz')
def main():
    lab=locked();bins=read(O/'tables/section_bins.tsv.gz');a,u=unit_bins(bins,lab);descriptive(a,u)
    sensitivity=read(C/'sensitivity_sections.tsv.gz');supports=[];tests=[];curves=[]
    for variant in ['primary','adjusted','fixed_units','stable','1000_T','x250_T','y250_T','xy250_T','500_I']:
        aa,uu=a,u
        if variant=='stable':aa,uu=unit_bins(bins,lab[lab.stability>=.8])
        if variant in sensitivity.rule.unique():
            ll=lab.drop(columns='group').merge(sensitivity[sensitivity.rule==variant][['section_id','cluster']].rename(columns={'cluster':'group'}),on='section_id');aa,uu=unit_bins(bins,ll)
        for b in B:
            for dom in (['E_inner','nonE_outer']+(['EIN_EIS_outer'] if b==B[1] else []) if variant=='primary' else ['E_inner']):
                orig=aa[(aa.boundary==b)&(aa.domain==dom)];d=uu[(uu.boundary==b)&(uu.domain==dom)]
                for gs in ([G]+[list(g) for g in itertools.combinations(G,2)] if variant in ['primary','adjusted','stable'] else [G]):
                    runmodel(d,orig,b,dom,gs,variant,supports,tests,curves)
    res=pd.DataFrame(tests)
    res['family']=res.apply(lambda r:'_'.join([r.variant,'global' if r.contrast.count(':')==2 else 'pairwise','E_negative' if r.domain=='E_inner' else r.domain,r.term]),axis=1)
    res['q']=np.nan
    for _,ix in res.groupby('family').groups.items():res.loc[ix,'q']=bh(res.loc[ix,'p'],len(ix))
    save(res,T/'curve_tests.tsv');save(res.query("family=='primary_global_E_negative_function'"),T/'four_primary_tests.tsv');save(curves,T/'model_differences.tsv.gz');save(supports,T/'curve_support.tsv.gz')
if __name__=='__main__':main()

