from common3 import *
from down import locked
lab=locked()[['section_id','dataset_id','unit','tumor','group']]
def paired_bins(a,variant):
    a=a[(a.domain=='E_inner') & (a.n>=5)]
    p=a.pivot(index=['section_id','center'],columns='boundary',values=['median','n']).dropna()
    p.columns=['_'.join(c) for c in p.columns]
    p=p.reset_index().merge(lab,on='section_id',validate='many_to_one')
    p['paired_diff']=p['median_'+B[1]]-p['median_'+B[0]]
    p['point_domain']=variant
    return p
v=read(O/'cache/geometry_bins.tsv.gz')
p=pd.concat([paired_bins(read(O/'tables/section_bins.tsv.gz'),'native_matched_sections'),paired_bins(v[v.variant=='common_points'],'same_E_points')],ignore_index=True)
save(p,OUT/'tables/paired_curve_bins.tsv.gz')
summ=[]
for (var,c,g),a in p[p.tumor&p.unit.notna()].groupby(['point_domain','center','group']):
    u=a.groupby(['dataset_id','unit']).paired_diff.mean()
    s=u.groupby(level=0).mean()
    se=s.sem() if len(s)>1 else np.nan
    half=stats.t.ppf(.975,len(s)-1)*se if len(s)>1 else np.nan
    summ.append(dict(point_domain=var,group=g,center=c,diff=s.mean(),lo=s.mean()-half,hi=s.mean()+half,n_sections=len(a),n_units=a.unit.nunique(),n_sources=len(s),sources='|'.join(s.index)))
save(summ,OUT/'tables/paired_curve_summary.tsv.gz')

