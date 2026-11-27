from lib import *
import argparse
from itertools import combinations
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');v=p.parse_args().version;out=O if v=='primary' else O/'versions'/v
m=read(out/'tables/all_section_levels.tsv')
sub=pd.concat([read(O/'metrics'/f'{s}_subtypes.tsv') for s in M.index])
sub=sub.merge(read(O/'audit/c2l_reference_usage.tsv')[['factor','major_lineage']],left_on='cell_type',right_on='factor',validate='m:1')
z=sub.groupby(['section_id','domain','major_lineage'],as_index=False).agg(n=('n','first'),abundance=('ab_sum','sum'),all30=('all30_sum','first'))
z['fraction']=z.abundance/z.all30.replace(0,np.nan);z=z.merge(m,on='section_id');save(z,out/'tables/major_lineage_by_section.tsv.gz');rows=[]
for (dom,k),a in z[z.source_context.ne('adjacent_normal_source')].dropna(subset=['analysis_unit_uid']).groupby(['domain','major_lineage']):
    for axis,groups in [('highlow',['high','low']),('provisional_architecture_group',G)]:
        u=a.groupby(['dataset_id','analysis_unit_uid',axis]).fraction.mean().unstack(axis)
        for g,h in combinations(groups,2):
            ds=u.reindex(columns=[g,h]).groupby('dataset_id').mean().dropna()
            for dataset,r in ds.iterrows():rows.append(dict(domain=dom,lineage=k,axis=axis,g1=g,g2=h,scope=dataset,effect=r[g]-r[h]))
            rows.append(dict(domain=dom,lineage=k,axis=axis,g1=g,g2=h,scope='overall',n_sources=len(ds),sources='|'.join(ds.index),status='estimable' if len(ds)>=3 else 'limited_or_no_common_sources',**unit_summary(ds[g]-ds[h])))
t=pd.DataFrame(rows);t['q']=np.nan
for _,ix in t[t.scope=='overall'].groupby(['domain','axis']).groups.items():t.loc[ix,'q']=bh(t.loc[ix,'p'],len(ix))
save(t,out/'tables/major_lineage_contrasts.tsv.gz')

