from common3 import *
from down import locked,unit_bins,runmodel
T=OUT/'tables';lab=locked();z=read(T/'section_levels.tsv');z=z[z.tumor&z.unit.notna()]
rows=[]
for ep in ['mean_MP10','median_MP10']:
    u=z.groupby(['dataset_id','group','unit'])[ep].mean().groupby(['dataset_id','group']).mean().unstack().dropna()
    d=np.c_[u.C2-u.C1,u.C3-u.C1];n,k=d.shape;mu=d.mean(axis=0);cov=np.cov(d.T)
    stat=(n-k)/(k*(n-1))*n*mu@np.linalg.inv(cov)@mu if n>k and np.linalg.matrix_rank(cov)==k else np.nan
    rows.append(dict(endpoint=ep,n_sources=n,sources='|'.join(u.index),C2_minus_C1=mu[0],C3_minus_C1=mu[1],F=stat,df1=k,df2=n-k,p=stats.f.sf(stat,k,n-k) if np.isfinite(stat) else np.nan))
r=pd.DataFrame(rows);r['q']=bh(r.p);save(r,T/'overall_level_global.tsv')
bins=read(O/'tables/section_bins.tsv.gz');a,u=unit_bins(bins,lab)
su=[];ts=[];pr=[]
for b in B:
    aa=a[(a.boundary==b)&(a.domain=='E_inner')];uu=u[(u.boundary==b)&(u.domain=='E_inner')]
    sources=uu.groupby('dataset_id').group.nunique();sources=sources[sources==3].index
    for ds in sources:runmodel(uu[uu.dataset_id!=ds],aa[aa.dataset_id!=ds],b,'E_inner',G,'omit_'+ds,su,ts,pr)
r=pd.DataFrame(ts);r['q']=np.nan
for _,ix in r.groupby(['variant','term']).groups.items():r.loc[ix,'q']=bh(r.loc[ix,'p'])
save(r,T/'curve_leave_one_source.tsv');save(pr,T/'curve_leave_one_source_differences.tsv.gz')
support=read(T/'curve_support.tsv.gz');emp=[];wide=[]
for (b,dom,co),s in support.query("variant=='primary' and connected_common").groupby(['boundary','domain','contrast']):
    gs=co.split(':');sources=s.both_or_all_common_sources.iloc[0].split('|');z=u[(u.boundary==b)&(u.domain==dom)&u.center.isin(s.center)&u.dataset_id.isin(sources)&u.group.isin(gs)]
    means=z.groupby(['center','group','dataset_id'])[['raw','centered']].mean().reset_index()
    for x,a in means.groupby('center'):
        for ep in ['raw','centered']:
            h=a.pivot(index='dataset_id',columns='group',values=ep).reindex(columns=gs).dropna();assert len(h)==len(sources)
            for g in gs:wide.append(dict(boundary=b,domain=dom,contrast=co,center=x,scale=ep,group=g,mean=h[g].mean(),n_sources=len(h)))
            for g1,g2 in itertools.combinations(gs,2):
                v=h[g1]-h[g2];ci=v.sem()*stats.t.ppf(.975,len(v)-1) if len(v)>1 else np.nan
                emp.append(dict(boundary=b,domain=dom,contrast=co,center=x,scale=ep,g1=g1,g2=g2,effect=v.mean(),lo=v.mean()-ci,hi=v.mean()+ci,n_sources=len(v)))
save(emp,T/'common_source_bin_differences.tsv.gz');save(wide,T/'common_source_bin_means.tsv.gz')

