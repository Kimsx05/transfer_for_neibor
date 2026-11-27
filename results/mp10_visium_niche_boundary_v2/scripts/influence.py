from lib import *
import argparse
p=argparse.ArgumentParser();p.add_argument('--version',default='primary');arg=p.parse_args();out=O if arg.version=='primary' else O/'versions'/arg.version
m=read(out/'tables/all_section_levels.tsv').query("source_context!='adjacent_normal_source'").dropna(subset=['analysis_unit_uid'])
t=read(O/'tables/section_effects.tsv.gz').query("variant=='native' and width==200 and eligible").merge(m,on='section_id')
rows=[]
for b,a in t.groupby('boundary'):
    for key in ['analysis_unit_uid','dataset_id']:
        for id in a[key].unique():
            v=a[a[key]!=id]
            for col in ['score','adjusted']:
                zz=v.groupby(['dataset_id','analysis_unit_uid'])[col].mean().dropna();ds=zz.groupby('dataset_id').mean()
                rows.append(dict(boundary=b,excluded_kind=key,excluded_id=id,endpoint=col,n_units=len(zz),n_sources=len(ds),**unit_summary(ds)))
save(rows,out/'tables/all_section_influence.tsv')
# Primary sample-matched architecture contrasts for each boundary and domain.
z=read(O/'tables/paired_boundaries.tsv.gz').query('eligible').merge(m,on='section_id')
rows=[]
for domain,a in z.groupby('domain'):
    for col in ['D_all','D_without','paired_difference']:
        for g1,g2 in [(G[0],G[1]),(G[0],G[2]),(G[1],G[2])]:
            u=a.groupby(['dataset_id','analysis_unit_uid','provisional_architecture_group'])[col].mean().unstack()
            h=u.reindex(columns=[g1,g2]).groupby('dataset_id').mean().dropna()
            rows.append(dict(domain=domain,endpoint=col,g1=g1,g2=g2,n_sources=len(h),sources='|'.join(h.index),status='ok' if len(h)>=3 else 'sparse_or_missing_common_source',**unit_summary(h[g1]-h[g2])))
save(rows,out/'tables/matched_architecture_comparisons.tsv')

