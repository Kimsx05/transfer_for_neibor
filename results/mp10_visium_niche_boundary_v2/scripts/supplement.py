from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args();out=O if arg.version=='primary' else O/'versions'/arg.version
def rd(f):return read(out/f)
def base(f):return read(O/f)
cl=rd('classification/section_classification.tsv')
w=rd('classification/window_states.tsv.gz')
loc=load('tables/window_component_metrics.parquet')
pq(loc.merge(w[['section_id','level','local_id','state_reason','local_state','neighbor_support']],how='left',on=['section_id','level','local_id'],validate='1:1'),out/'tables/window_component_metrics_with_state.parquet')
m=base('tables/section_metrics_all.tsv.gz')
statecols=['F_entry','F_non','F_unresolved','window_share_entry','window_share_non','effective_area_coverage','H_supported']
extra=[]
for _,r in cl.iterrows():
    a=w[(w.section_id==r.section_id)&(w.level=='window500')];cfg=json.loads((out/'classification/classification_lock.json').read_text())['config'];ok=a.E_n>=cfg['min_E_500'];nE=a.loc[ok,'E_n'].sum();nt=a.E_n.sum();n=ok.sum();i=a.loc[(a.local_state=='infiltrated')&a.neighbor_support,'E_n'].sum();l=a.loc[(a.local_state=='non_infiltrated')&a.neighbor_support,'E_n'].sum()
    for col,num,den,form in zip(statecols,[i,l,nE-i-l,r.n_entry,r.n_non,i+l,r.H_supported],[nE,nE,nE,n,n,nt,1],['entry E sites / E sites in eligible windows','limited entry E sites / eligible E sites','unresolved eligible E sites / eligible E sites','supported entry windows / eligible windows','supported limited-entry windows / eligible windows','supported E sites / all E sites','−sum p log2(p), among supported entry vs limited entry E area; counts in section_classification and windows']):
        extra.append(dict(section_id=r.section_id,metric=col,value=r[col],numerator=num,denominator=den,domain='all_E',scale='500 μm',formula=form,unit='bits' if col=='H_supported' else 'fraction',status='ok' if np.isfinite(r[col]) else 'undefined'))
extra=pd.DataFrame(extra);save(pd.concat([m,extra],ignore_index=True),out/'tables/section_metrics_with_states.tsv.gz')
ct=base('audit/metric_dictionary.tsv');more=extra.drop_duplicates('metric')[['metric','domain','scale','unit','formula']].assign(direction='larger named area/window share or higher pattern heterogeneity',missing_rule='zero eligible area/windows or no supported states for entropy',epsilon='none')
save(pd.concat([ct,more],ignore_index=True),out/'audit/metric_dictionary_with_states.tsv')
ct=rd('tables/group_contrasts.tsv.gz');t=ct.query("variant=='native' and width==200 and group_rule=='architecture' and scope!='overall' and endpoint=='score'").dropna(subset=['effect']);res=[]
for (b,g,h),v in t.groupby(['boundary','g1','g2']):
    for ds in v.scope:
        a=v[v.scope!=ds];res.append(dict(boundary=b,g1=g,g2=h,excluded_source=ds,n_sources=len(a),**unit_summary(a.effect)))
save(res,out/'tables/leave_one_source.tsv')
# Compare fractional architecture membership among high / low units with the same contributing origins.
u=rd('tables/highlow_architecture_units.tsv');gcols=G+['uncertain'];z=u.pivot_table(index=['dataset_id','analysis_unit_uid','highlow'],columns='architecture',values='weight',fill_value=0).reindex(columns=gcols,fill_value=0).reset_index()
reg=base('tables/region_niche.tsv.gz');levels=rd('tables/all_section_levels.tsv')
sub=pd.concat([read(O/'metrics'/f'{s}_subtypes.tsv') for s in M.index])
measures={'level_mean':levels[['section_id','mean_MP10']].rename(columns={'mean_MP10':'value'}),'level_median':levels[['section_id','median_MP10']].rename(columns={'median_MP10':'value'})}
for (typ,cat),v in reg.groupby(['type','category']):measures[typ+'_'+str(cat)]=v[['section_id','fraction']].rename(columns={'fraction':'value'})
for (dom,typ),v in sub.groupby(['domain','cell_type']):measures['subtype_'+dom+'_'+typ]=v[['section_id','fraction']].rename(columns={'fraction':'value'})
tabs={k:z[['dataset_id','analysis_unit_uid','highlow',k]].rename(columns={k:'value'}) for k in gcols}
for k,v in measures.items():
    t=v.merge(levels[['section_id','analysis_unit_uid','dataset_id','highlow']],on='section_id').dropna(subset=['highlow','analysis_unit_uid','value']);tabs[k]=t.groupby(['dataset_id','analysis_unit_uid','highlow'],as_index=False).value.mean()
rows=[]
for name,a in tabs.items():
    vals=[]
    for ds,v in a.groupby('dataset_id'):
        left=v[v.highlow=='high'].value;right=v[v.highlow=='low'].value;d=left.mean()-right.mean()
        if np.isfinite(d):vals.append(d)
        rows.append(dict(measure=name,scope=ds,n1=len(left),n2=len(right),effect=d))
    rows.append(dict(measure=name,scope='overall',n_sources=len(vals),**unit_summary(vals)))
t=pd.DataFrame(rows);t['status']='descriptive_or_estimable'
for grp in gcols:
    if z[grp].sum()==0:
        keep=(t.measure==grp)&(t.scope=='overall');t.loc[keep,['low','high','p']]=np.nan;t.loc[keep,'status']='no_observed_members_descriptive_zero_share'
t['family']=[('architecture' if k in gcols else k.split('_')[0]) for k in t.measure];t['q']=np.nan
for _,idx in t[t.scope=='overall'].groupby('family').groups.items():t.loc[idx,'q']=bh(t.loc[idx,'p'],len(idx))
save(t,out/'tables/highlow_levels_composition.tsv.gz')
# Representative corpus for an offline scientific appendix: scatter of every unit and source-specific OLS display only.
pts=rd('plot_data/association_points.tsv.gz');sections=rd('plot_data/section_all_metrics_levels.tsv.gz');sections=sections[~sections.source_context.eq('adjacent_normal_source')]
ids=sorted(pts.analysis_unit_uid.unique());sources=sorted(sections.dataset_id.unique());inputs={}
for (metric,end),v in pts.groupby(['metric','endpoint']):
    inputs.setdefault(metric,{})[end]=[[ids.index(r.analysis_unit_uid),sources.index(r.dataset_id),round(r.x,7),round(r.y,7)] for r in v.itertuples()]
for metric in inputs:
    inputs[metric]['section_distribution']=[[r.section_id,sources.index(r.dataset_id),round(float(getattr(r,metric)),7)] for r in sections[['section_id','dataset_id',metric]].dropna().itertuples()]
payload=json.dumps(dict(units=ids,sources=sources,metrics=inputs),separators=(',',':'),allow_nan=False)
(out/'metric_explorer.html').write_text((O/'scripts/explorer_template.html').read_text().replace('/*DATA*/',payload))
print('explorer',len(inputs),'metrics',len(pts),'unit-points')

