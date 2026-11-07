from common import *
def main():
 mapping=read('canonical_programs.tsv');mapfull=dict(zip(mapping.full_program,mapping.canonical_program));mapno=dict(zip(mapping.noMP10_program,mapping.canonical_program));rows=[]
 for sec in T.section_id:
  pieces=[]
  for path in [OLD/'tables'/f'{sec}_rings.tsv.gz',S/'tables'/f'{sec}_new5_reference_rings.tsv.gz']:
   if path.exists():
    try:d=pd.read_csv(path,sep='\t')
    except pd.errors.EmptyDataError:continue
    if 'ring' in d:pieces.append(d[d.ring.eq('R0')&d.population.eq('Epi')])
  if not pieces:continue
  d=pd.concat(pieces);meta=T.set_index('section_id').loc[sec];rid=meta.analysis_unit_uid if pd.notna(meta.analysis_unit_uid) else 'section_only::'+sec
  for ver,lookup in [('full',mapfull),('noMP10',mapno)]:
   q=d[d.program.isin(lookup)].copy();q['feature']=q.program.map(lookup)
   for (pr,feature),z in q.groupby(['program','feature']):
    z=z.set_index('reference');p=z.loc['P95'] if 'P95' in z.index else None
    if p is not None:rows.append(dict(input_family='UCell_P95_R0_mean_'+ver,record_id=rid,dataset_id=meta.dataset_id,feature=feature,value=p['mean'],n_observed=p['n'],fraction_unlabelled=np.nan,transform='standard'))
    if {'P95_matched_base','control_base'}.issubset(z.index):rows.append(dict(input_family='UCell_P95_matched_R0_difference_'+ver,record_id=rid,dataset_id=meta.dataset_id,feature=feature,value=z.loc['P95_matched_base','mean']-z.loc['control_base','mean'],n_observed=min(z.loc['P95_matched_base','n'],z.loc['control_base','n']),fraction_unlabelled=np.nan,transform='standard'))
 new=pd.DataFrame(rows).groupby(['input_family','record_id','dataset_id','feature','transform'],dropna=False).agg(value=('value','mean'),n_observed=('n_observed','sum'),fraction_unlabelled=('fraction_unlabelled','mean')).reset_index();old=read('unit_background.tsv.gz');old=old[~old.input_family.str.startswith('UCell_P95')];save(pd.concat([old,new],ignore_index=True),'unit_background.tsv.gz');print('P95 state rows',len(new),flush=True)
if __name__=='__main__':main()
