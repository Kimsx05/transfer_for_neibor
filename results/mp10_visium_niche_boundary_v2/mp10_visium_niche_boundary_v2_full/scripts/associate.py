from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args()
out=O if arg.version=='primary' else O/'versions'/arg.version
c=read(out/'classification/section_classification.tsv').set_index('section_id')
mx=read(O/'tables/section_metrics_all.tsv.gz').pivot(index='section_id',columns='metric',values='value')
features=['F_entry','F_non','F_unresolved','window_share_entry','window_share_non','effective_area_coverage','H_supported']
mx=mx.join(c[features])
e=read(O/'tables/section_effects.tsv.gz').query("variant=='native' and width==200 and eligible").pivot(index='section_id',columns='boundary',values='score').add_prefix('D_')
lev=read(O/'tables/section_exposures.tsv.gz').query("rule=='all'").set_index('section_id')
y=e.join(lev[['mean_MP10','median_MP10']],how='outer')
full=mx.join(M[['dataset_id','analysis_unit_uid','source_context']]).join(y).join(c.provisional_architecture_group)
groups=read(out/'tables/unit_exposure_groups.tsv').query("rule=='all'")[['dataset_id','analysis_unit_uid','highlow','old_highlow','dataset_median']]
save(full.reset_index().merge(groups,how='left',on=['dataset_id','analysis_unit_uid']),out/'plot_data/section_all_metrics_levels.tsv.gz')
save(full[['dataset_id','analysis_unit_uid','source_context','provisional_architecture_group']].join(lev[['n','mean_MP10','median_MP10']]).reset_index().merge(groups,how='left',on=['dataset_id','analysis_unit_uid']),out/'tables/all_section_levels.tsv')
full=full[~full.source_context.eq('adjacent_normal_source')];dslist=full.dataset_id.unique()
rows=[];pts=[];dists=[]
for met in mx:
    for ds,a in full.groupby('dataset_id'):
        x=a[met].dropna()
        dists.append(dict(metric=met,source=ds,total=len(a),available=len(x),q10=x.quantile(.1),median=x.median(),q90=x.quantile(.9)))
    for end in y:
        z=full[[met,end,'dataset_id','analysis_unit_uid']].dropna().groupby(['dataset_id','analysis_unit_uid'],as_index=False)[[met,end]].mean().rename(columns={met:'x',end:'y'})
        pts.append(z.assign(metric=met,endpoint=end))
        vals=[]
        for ds in dslist:
            a=z[z.dataset_id==ds];n=len(a);rho=lo=hi=pv=np.nan
            ok=n>=4 and a.x.nunique()>1 and a.y.nunique()>1
            if ok:
                rho,pv=stats.spearmanr(a.x,a.y);rz=np.arctanh(np.clip(rho,-.999999,.999999));lo,hi=np.tanh(rz+np.array([-1,1])*1.96/np.sqrt(n-3));vals.append(rho)
            rows.append(dict(metric=met,endpoint=end,scope=ds,n_units=n,n_sources=int(ok),effect=rho,low=lo,high=hi,p=pv,status='ok' if ok else 'low_n_or_constant'))
        summary=unit_summary(vals)
        rows.append(dict(metric=met,endpoint=end,scope='overall',n_units=len(z[z.dataset_id.isin([r['scope'] for r in rows[-len(dslist):] if r['status']=='ok'])]),n_sources=len(vals),status='ok' if len(vals)>=3 else 'limited_source_support',**summary))
a=pd.DataFrame(rows);a['q']=np.nan
for _,idx in a[a.scope=='overall'].groupby('endpoint').groups.items():a.loc[idx,'q']=bh(a.loc[idx,'p'],len(idx))
save(a,out/'tables/all_metric_associations.tsv.gz');save(pd.concat(pts),out/'plot_data/association_points.tsv.gz');save(dists,out/'tables/metric_distributions_by_source.tsv.gz')
print(len(mx.columns),'metrics',len(rows),'correlations',flush=True)

