from down import locked,draws,globtest
from common3 import *
C=OUT/'classification';T=OUT/'tables'
def contrasts(z,value,measure,family):
    rows=[]
    for g,h in itertools.combinations(G,2):
        a=z[z.group.isin([g,h])].dropna(subset=['unit',value]).groupby(['dataset_id','group','unit'],as_index=False)[value].mean()
        ds=[s for s,v in a.groupby('dataset_id') if v.group.nunique()==2];a=a[a.dataset_id.isin(ds)]
        info=dict(measure=measure,family=family,g1=g,g2=h,sources='|'.join(ds),n_sources=len(ds),n1=a[a.group==g].unit.nunique(),n2=a[a.group==h].unit.nunique())
        if not len(ds):rows.append(dict(info,scope='common_sources',status='no_common_source'));continue
        ids,wt=draws(a.unit,a.dataset_id);means=[]
        for name,v in a.groupby('dataset_id'):
            gr=[]
            for label in [g,h]:
                s=v[v.group==label];w=wt[:,ids.get_indexer(s.unit)];den=w.sum(1);gr.append(np.divide(w@s[value].to_numpy(),den,out=np.full(len(den),np.nan),where=den>0))
            means.append(gr[0]-gr[1]);rows.append(dict(info,scope=name,effect=means[-1][0],status='source_description'))
        x=np.mean(means,axis=0);bs=x[1:][np.isfinite(x[1:])];support=min(info['n1'],info['n2'])>=5 and len(ds)>=2
        # t interval and test over paired source contrasts, avoiding conditioning on nonempty arms in sparse unit resamples.
        vv=np.array(means)[:,0];se=stats.sem(vv) if len(vv)>1 else np.nan;half=stats.t.ppf(.975,len(ds)-1)*se if len(ds)>1 else np.nan
        pv=2*stats.t.sf(abs(x[0]/se),len(ds)-1) if support and se>0 else np.nan
        rows.append(dict(info,scope='common_sources',status='ok' if support else 'limited',effect=x[0],lo=x[0]-half,hi=x[0]+half,p=pv,finite_unit_resamples=len(bs),boot_lo=np.quantile(bs,.025) if len(bs)>=400 else np.nan,boot_hi=np.quantile(bs,.975) if len(bs)>=400 else np.nan))
    return rows
def main():
    cl=locked();lab=cl[cl.tumor][['section_id','dataset_id','unit','group','stability']];lev=read(V/'tables/all_section_levels.tsv')[['section_id','mean_MP10','median_MP10','highlow','old_highlow','n']]
    lev=lev.merge(cl[['section_id','group','tumor','dataset_id','unit','stability']],on='section_id');save(lev,T/'section_levels.tsv')
    tests=[];counts=[];e=read(V/'tables/section_effects.tsv.gz').merge(lab,on='section_id');save(e,T/'v2_effects_relabelled.tsv.gz')
    for value in ['mean_MP10','median_MP10']:tests+=contrasts(lev[lev.tumor],value,value,'overall_levels')
    h=lev[lev.tumor&lev.highlow.isin(['high','low'])].copy();h['is_high']=h.highlow.eq('high').astype(int);tests+=contrasts(h,'is_high','high_proportion_excluding_ties','highlow')
    save(lev[lev.tumor].groupby(['dataset_id','group','highlow'],dropna=False).agg(n_sections=('section_id','size'),n_units=('unit','nunique')).reset_index(),T/'highlow_counts.tsv')
    for keys,a in e[e.width.eq(200)&e.eligible&e.variant.isin(['native','strict','same_component','common_points'])].groupby(['boundary','variant']):
        for value in (['score','adjusted','P95_pp'] if keys[1]=='native' else ['score']):
            tests+=contrasts(a,value,'_'.join([*keys,value]),'D_'+keys[1]+'_'+value)
    pa=read(V/'tables/paired_boundaries.tsv.gz').merge(lab,on='section_id');save(pa,T/'paired_native_common.tsv.gz')
    for dom,a in pa[pa.eligible].groupby('domain'):
        tests+=contrasts(a,'paired_difference',dom+'_paired_boundary_difference','paired_geometry')
        for group,v in a.groupby('group'):
            u=v.groupby(['dataset_id','unit']).paired_difference.mean().groupby('dataset_id').mean();se=u.sem();half=se*stats.t.ppf(.975,len(u)-1) if len(u)>1 else np.nan
            counts.append(dict(domain=dom,group=group,n_sections=len(v),n_units=v.unit.nunique(),n_sources=len(u),effect=u.mean(),lo=u.mean()-half,hi=u.mean()+half))
    save(counts,T/'paired_geometry_summary.tsv')
    reg=read(V/'tables/region_niche.tsv.gz').merge(lab,on='section_id');save(reg,T/'region_niche_by_section.tsv.gz')
    for (typ,cat),a in reg.groupby(['type','category']):
        tests+=contrasts(a,'fraction',typ+'_'+str(cat),'unknown_QC' if cat=='Unknown' else typ)
    rr=pd.DataFrame(tests);rr['q']=np.nan
    for _,idx in rr[rr.scope=='common_sources'].groupby('family').groups.items():rr.loc[idx,'q']=bh(rr.loc[idx,'p'],len(idx))
    save(rr,T/'auxiliary_tests.tsv.gz')
    # Four pre-defined endpoints: overall mean, median, signed near-to-deep gradient over original E median bins of each mask.
    bs=read(O/'tables/section_bins.tsv.gz');sh=[]
    for (sec,b),a in bs.query("domain=='E_inner' and n>=5 and center>=-450").groupby(['section_id','boundary']):
        a=a.sort_values('center',ascending=False);a=a.iloc[:next((i for i,v in enumerate(a.center) if v!=-50-i*100),len(a))]
        if len(a)>=3:sh.append(dict(section_id=sec,boundary=b,slope=np.polyfit(a.center/100,a['median'],1)[0]))
    slope=pd.DataFrame(sh).pivot(index='section_id',columns='boundary',values='slope').add_prefix('slope_')
    z=lev[lev.tumor].set_index('section_id').join(slope).join(read(O/'cache/immune_features.tsv').query("axis==@AXIS and grid==@GRID").set_index('section_id')[list('ABCDEFG')])
    keep=read(C/'feature_selection.tsv').query("reason=='kept'").feature;out=[]
    for feat,end in itertools.product(keep,['mean_MP10','median_MP10']+list(slope.columns)):
        u=z.dropna(subset=[feat,end,'unit']).groupby(['dataset_id','unit'])[[feat,end]].mean().reset_index();rhos=[]
        for ds,a in u.groupby('dataset_id'):
            v,p=stats.spearmanr(a[feat],a[end]) if len(a)>=4 and a[feat].nunique()>1 and a[end].nunique()>1 else (np.nan,np.nan)
            if np.isfinite(v):rhos.append(v)
            out.append(dict(feature=feat,endpoint=end,scope=ds,n_units=len(a),rho=v,p=p))
        se=stats.sem(rhos) if len(rhos)>1 else np.nan;mu=np.mean(rhos);ci=se*stats.t.ppf(.975,len(rhos)-1) if len(rhos)>1 else np.nan;pv=2*stats.t.sf(abs(mu/se),len(rhos)-1) if len(rhos)>=3 and se>0 else np.nan
        out.append(dict(feature=feat,endpoint=end,scope='source_balanced',rho=mu,p=pv,lo=mu-ci,hi=mu+ci,n_sources=len(rhos)))
    res=pd.DataFrame(out);res['q']=np.nan;ix=res[res.scope=='source_balanced'].index;res.loc[ix,'q']=bh(res.loc[ix,'p']);save(res,T/'continuous_associations.tsv');save(z.reset_index(),T/'continuous_analysis_data.tsv')
if __name__=='__main__':main()

