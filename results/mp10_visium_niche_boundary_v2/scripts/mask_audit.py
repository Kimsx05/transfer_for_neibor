from lib import *
rows=[]
for s in M.index:
    d=load(f'spots/{s}.parquet');old=pd.read_parquet(ROOT/'results/mp10_p95_front_programs_v1/spots'/f'{s}.parquet',columns=['barcode','epi_anchor','lesion_location','query'])
    a=d[d.all_E].merge(old,on='barcode',validate='1:1');assert a.legacy_query.eq(a['query']).all()
    predicted=a.epi_anchor & a.score_before & a.lesion_location.ne('quality_excluded');assert predicted.eq(a.legacy_query).all()
    ct=a.groupby(['epi_anchor','score_before','lesion_location','legacy_query','new_epi_domain'],dropna=False).size().reset_index(name='n')
    rows.append(ct.assign(section_id=s))
save(pd.concat(rows),'audit/old_display_mask_reasons.tsv')

