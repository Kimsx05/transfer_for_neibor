# Regenerate regional measurements with already frozen parameters; no parameter selection.
from common import *
c=yaml.safe_load((R/'frozen_config.yaml').read_text());sig=c['signatures'];gate()
for s in SECTIONS:
 for kind in ['T32','T16','Epi64']:
  d=pd.read_parquet(R/s/f'{kind}_measure.parquet');rank=d.filter(like='rank_').rename(columns=lambda x:x[5:])
  for k,g in sig.items():d[k]=ucell(rank,g,c['UCell_maxRank'])
  d['MP10']=mp10(rank);d['Treg_main_detected']=(d[['raw_'+g for g in sig['Treg_main']]]>0).sum(1);d['Treg_main_UMI']=d[['raw_'+g for g in sig['Treg_main']]].sum(1);d['Treg_core_detected']=(d[['raw_'+g for g in ['FOXP3','IL2RA','CTLA4','IKZF2']]]>0).sum(1)
  if kind.startswith('T'):
   d['measurement_eligible']=(d.nFeature>=c['T_min_nFeature'])&(d.T_detected>=c['T_lineage_detected_min'])&(d.T_UMI>=c['T_lineage_UMI_min'])
   d['evidence']=np.where(~d.measurement_eligible,'insufficient',np.where((d.Treg_main_detected>=c['Treg_like_main_detected_min'])&(d.Treg_core_detected>=c['Treg_like_core_detected_min'])&(d.Treg_main>c['Treg_like_threshold_training_competitor_q95']),'Treg_like_multigene','T_lineage_or_activated'))
  else:
   d['measurement_eligible']=(d.nFeature>=c['Epi_min_nFeature'])&(d.MP10_detected>=c['Epi_min_MP10_detected'])&(d.Epi_detected>=2);d['evidence']=np.where(d.measurement_eligible,'Epi_enriched_scoreable','insufficient');q=d.loc[d.measurement_eligible,'MP10'].quantile([.5,.9,.95]);d['L50']=d.measurement_eligible&(d.MP10<=q.loc[.5])
   for name,qq in [('H50',.5),('H10',.9),('H05',.95)]:d[name]=d.measurement_eligible&(d.MP10>q.loc[qq])
  d['mixed_Epi_Fib']=(d.Epi_detected>=2)&(d.Fib_detected>=2);d.drop(columns=d.filter(like='rank_').columns).to_parquet(R/s/f'{kind}_regions.parquet',index=False)
print('RESTORED_FROZEN_SCORES')
