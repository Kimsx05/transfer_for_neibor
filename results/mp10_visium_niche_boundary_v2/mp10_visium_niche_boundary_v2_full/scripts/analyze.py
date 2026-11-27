from lib import *
from itertools import combinations
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--version',default='primary');arg=parser.parse_args()
OUT=O if arg.version=='primary' else O/'versions'/arg.version
cl=read(OUT/'classification/section_classification.tsv');lock=json.loads((OUT/'classification/classification_lock.json').read_text())
assert lock['class_sha256']==sha(OUT/'classification/section_classification.tsv')
def emit(rows,name):save(rows,OUT/name)
def table(name):return read(O/'tables'/f'{name}.tsv.gz')
meta=M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid']].reset_index(drop=True).merge(cl[['section_id','provisional_architecture_group','loose','strict','scale1000']],on='section_id').rename(columns={'provisional_architecture_group':'architecture'})
meta=meta[meta.source_context!='adjacent_normal_source'];uk=['dataset_id','analysis_unit_uid']
def eligible_meta(d):return d.merge(meta,on='section_id',validate='m:1').dropna(subset=['analysis_unit_uid'])
def paired_sources(d,col,group,levels):
    # Calculate source-wise contrasts of unit means. Shared units enter both groups;
    # cluster influence combines their two contributions before squaring.
    rows=[]
    for g1,g2 in combinations(levels,2):
        common=[]
        for ds,a in d.groupby('dataset_id'):
            z=a.dropna(subset=[col]).groupby(['analysis_unit_uid',group])[col].mean().unstack(group)
            v1=z[g1].dropna() if g1 in z else pd.Series(dtype=float);v2=z[g2].dropna() if g2 in z else pd.Series(dtype=float)
            diff=v1.mean()-v2.mean();n1=len(v1);n2=len(v2);uu=v1.index.union(v2.index);n=len(uu)
            se=low=high=p=np.nan
            if min(n1,n2)>=2 and n>=3:
                inf=((v1-v1.mean())/n1).reindex(uu,fill_value=0)-((v2-v2.mean())/n2).reindex(uu,fill_value=0)
                se=np.sqrt(sum(inf**2)*n/(n-1));crit=stats.t.ppf(.975,n-1)
                low=diff-crit*se;high=diff+crit*se;p=2*stats.t.sf(abs(diff/se),n-1) if se>0 else np.nan
            rows.append(dict(g1=g1,g2=g2,scope=ds,n1=n1,n2=n2,shared_units=len(v1.index.intersection(v2.index)),effect=diff,low=low,high=high,p=p,status='missing_arm' if not min(n1,n2) else 'limited_replication' if min(n1,n2)<2 else 'ok'))
            if np.isfinite(diff):common.append((ds,diff,n1,n2))
        z=[v[1] for v in common]
        rows.append(dict(g1=g1,g2=g2,scope='overall',n_sources=len(z),n1=sum(v[2] for v in common),n2=sum(v[3] for v in common),sources='|'.join(v[0] for v in common),status='no_common_source' if not z else 'fewer_than_three_sources' if len(z)<3 else 'ok',lodo_min=min([np.mean(np.delete(z,i)) for i in range(len(z))],default=np.nan) if len(z)>1 else np.nan,lodo_max=max([np.mean(np.delete(z,i)) for i in range(len(z))],default=np.nan) if len(z)>1 else np.nan,**unit_summary(z)))
    return rows
def by_group(d,col):
    u=d.groupby(uk+['architecture'],as_index=False)[col].mean().dropna()
    rows=[]
    for g in G+['uncertain','all']:
        a=u[u.architecture==g] if g!='all' else d.groupby(uk,as_index=False)[col].mean().dropna()
        ds=a.groupby('dataset_id')[col].mean()
        rows.extend(dict(group=g,scope=k,n_units=len(z),**unit_summary(z[col])) for k,z in a.groupby('dataset_id'))
        rows.append(dict(group=g,scope='overall',n_units=len(a),n_sources=len(ds),**unit_summary(ds)))
    return rows
ex=eligible_meta(table('section_exposures'))
u=ex.groupby(['rule']+uk,as_index=False).median_MP10.mean()
u['dataset_median']=u.groupby(['dataset_id','rule']).median_MP10.transform('median');u['highlow']=np.select([u.median_MP10.isna(),u.median_MP10>u.dataset_median,u.median_MP10<u.dataset_median],['unavailable','high','low'],default='tie')
old=read(O/'audit/v1_unit_exposures.tsv').query("rule=='all'")[uk+['group']].rename(columns={'group':'old_highlow'})
u=u.merge(old,on=uk,how='left');emit(u,'tables/unit_exposure_groups.tsv')
e=eligible_meta(table('section_effects'));e=e[e.eligible]
tests=[];summ=[];hl_tests=[]
for (b,v,w),d in e.groupby(['boundary','variant','width']):
    for col in ['score','adjusted','P95_pp']:
        if v=='border_spot' and col=='P95_pp' or v=='same_component' and col!='score':continue
        if v=='legacy_domain' and col=='P95_pp':continue
        vals=d[d.P95_eligible.eq(True)] if col=='P95_pp' else d
        key=dict(boundary=b,variant=v,width=w,endpoint=col)
        summ+= [dict(key,**x) for x in by_group(vals,col)]
        for rule in (['architecture','loose','strict','scale1000'] if (v,w,col)==('native',200,'score') else ['architecture']):
            rr=paired_sources(vals,col,rule,G)
            tests.extend(dict(key,group_rule=rule,**x) for x in rr)
    if (v,w)==('native',200):
        for rule in ['all','old','exclude_B1_'+b]:
            ref=u[u.rule==('all' if rule=='old' else rule)]
            a=d.merge(ref[uk+['highlow','old_highlow']],on=uk)
            hl_tests.extend(dict(boundary=b,exposure_rule=rule,**x) for x in paired_sources(a,'score','old_highlow' if rule=='old' else 'highlow',['high','low']))
test=pd.DataFrame(tests);test['q']=np.nan
# Six planned pairwise hypotheses per width/endpoint/variant, including empty arms.
for _,idx in test[test.scope=='overall'].groupby(['variant','width','endpoint','group_rule']).groups.items():test.loc[idx,'q']=bh(test.loc[idx,'p'],6)
emit(test,'tables/group_contrasts.tsv.gz');emit(test.query("variant=='native' and width==200 and endpoint=='score' and group_rule=='architecture' and scope=='overall'"),'tables/main_six_tests.tsv')
emit(summ,'tables/group_effects.tsv.gz');hl=pd.DataFrame(hl_tests);hl['q']=np.nan
for _,idx in hl[hl.scope=='overall'].groupby('exposure_rule').groups.items():hl.loc[idx,'q']=bh(hl.loc[idx,'p'],2)
emit(hl,'tables/highlow_boundary.tsv')
pair=eligible_meta(table('paired_boundaries'));out=[]
for dom,d in pair[pair.eligible].groupby('domain'):
    for col in ['D_all','D_without','paired_difference']:out.extend(dict(domain=dom,endpoint=col,**x) for x in by_group(d,col))
emit(out,'tables/paired_summary.tsv')
# Distinguish overall level, abundance within regions and compositional contrasts.
rows=[];levels=[];comp=eligible_meta(table('region_niche'));sub=eligible_meta(pd.concat([read(O/'metrics'/f'{sec}_subtypes.tsv') for sec in M.index]))
for key,d,col in [('overall_mean',ex[ex.rule=='all'],'mean_MP10'),('overall_median',ex[ex.rule=='all'],'median_MP10')]+[(k+'_'+str(n),v,'fraction') for (k,n),v in comp.groupby(['type','category'])]+[('within_all_E_'+n,v,'within_E') for n,v in comp[comp.type.eq('region')&comp.category.isin(['EIN','EIS','ENN','ENS'])].groupby('category')]+[('score_in_'+n,v,'score') for n,v in comp[comp.type.eq('region')&comp.category.isin(['EIN','EIS','ENN','ENS'])].groupby('category')]+[('subtype_'+dom+'_'+typ,v,'fraction') for (dom,typ),v in sub.groupby(['domain','cell_type'])]:
    rows.extend(dict(measure=key,**x) for x in paired_sources(d,col,'architecture',G))
    levels.extend(dict(measure=key,**x) for x in by_group(d,col))
out=pd.DataFrame(rows)
out['family']=out.measure.str.split('_').str[0];out.loc[out.measure=='region_Unknown','family']='Unknown_QC';out['q']=np.nan
for _,idx in out[out.scope=='overall'].groupby('family').groups.items():out.loc[idx,'q']=bh(out.loc[idx,'p'],len(idx))
emit(out,'tables/level_composition_tests.tsv.gz');emit(levels,'tables/level_composition_summary.tsv.gz')
a=ex[ex.rule=='all'].merge(u[u.rule=='all'][uk+['highlow','old_highlow']],on=uk);emit(a,'tables/section_levels_and_groups.tsv')
# A unit may contribute different section phenotypes. Allocate its one vote fractionally.
def frequency(g):
    z=g.groupby(uk+['highlow','architecture']).size().rename('s').reset_index()
    z['weight']=z.s/z.groupby(uk).s.transform('sum')
    return z
z=frequency(a);emit(z,'tables/highlow_architecture_units.tsv')
z=z.groupby(['dataset_id','highlow','architecture']).weight.sum().reset_index()
grid=pd.MultiIndex.from_product([a.dataset_id.unique(),['high','low','tie'],G+['uncertain','not_evaluable']],names=['dataset_id','highlow','architecture'])
z=z.set_index(grid.names).reindex(grid,fill_value=0).reset_index();z['den']=z.groupby(['dataset_id','highlow']).weight.transform('sum');z['fraction']=z.weight/z.den.replace(0,np.nan)
emit(z,'tables/highlow_architecture_proportions.tsv')
# Curves: smooth within section on consecutive measured centers, separately inside and EIN/EIS.
t=table('section_bins');t=t[t.n>=5].copy();parts=[]
for _,a in t.groupby(['section_id','boundary','side']):
    a=a.sort_values('center').copy()
    for col in ['score','centered']:
        h=a.set_index('center')[col]
        a[col+'_smooth']=[np.average([h.loc[i] for i in [x-100,x,x+100] if i in h],weights=[2 if i==x else 1 for i in [x-100,x,x+100] if i in h]) for x in a.center]
    parts.append(a)
t=pd.concat(parts);emit(t,'plot_data/section_curves.tsv.gz');t=eligible_meta(t);vals=['score','centered','score_smooth','centered_smooth','P95']
u2=t.groupby(uk+['boundary','center','architecture'],as_index=False)[vals].mean()
emit(u2,'plot_data/unit_curves.tsv.gz');curves=[];supports=[]
for (b,x),a in u2.groupby(['boundary','center']):
    for g in G+['uncertain']:
        aa=a[a.architecture==g];ds=aa.groupby('dataset_id')[vals].mean()
        scopes=t[(t.architecture==g)&(t.boundary==b)&(t.center==x)]
        for metric in vals:curves.append(dict(boundary=b,center=x,group=g,metric=metric,n_sections=scopes.section_id.nunique(),n_units=len(aa),n_sources=len(ds),sources='|'.join(ds.index),**unit_summary(ds[metric])))
    for g1,g2 in combinations(G,2):
        h=a[a.architecture==g1].groupby('dataset_id')[vals].mean();l=a[a.architecture==g2].groupby('dataset_id')[vals].mean();both=h.index.intersection(l.index)
        for col in ['score','centered']:supports.append(dict(boundary=b,center=x,g1=g1,g2=g2,metric=col,n_sources=len(both),sources='|'.join(both),n_units=a[a.dataset_id.isin(both)&a.architecture.isin([g1,g2])].analysis_unit_uid.nunique(),**unit_summary(h.loc[both,col]-l.loc[both,col])))
emit(curves,'plot_data/group_curves.tsv.gz');emit(supports,'tables/shared_source_curves.tsv.gz')

