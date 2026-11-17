from common import *
from analysis_stats import fit_effect,bh
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());b=pd.read_parquet(R/'Cervilla/T32_neighborhoods.parquet');c=pd.read_parquet(R/'Cervilla/cell_region_neighborhood_comparison.parquet');rows=[]
for r in [50,200]:
 bb=b[b.radius_um.eq(r)&b.source_rule.eq('exclude_shared')&b.association_eligible].set_index('T_region_id');cc=c[c.radius_um.eq(r)&c.association_eligible].set_index('T_region_id');ids=bb.index.intersection(cc.index)
 for unit,d in [('bin_enriched',bb.loc[ids]),('cell_limited',cc.loc[ids])]:
  for block in [cfg['block_um'],cfg['larger_block_um']]:
   z=fit_effect(d,'Treg_main',block_um=block,covars=['T_lineage','T_nCount','Epi_nCount','observed_circle_fraction','T_background_balance']);rows.append(dict(section='Cervilla',radius_um=r,unit=unit,comparison='same_T_centers_same_covariates',**z))
put(bh(pd.DataFrame(rows),['comparison','block_um']),'cell_bin_common_effects.tsv')
