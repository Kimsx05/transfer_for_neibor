from common import *
from analysis_stats import fit_effect,bh
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());gate();results=[];coverage=[]
(R/'models/association_started.json').write_text(json.dumps({'started':datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),'configuration_frozen':cfg['frozen_at'],'source':'source-disjoint main; shared-bin inclusion diagnostic'},indent=2))
for s in SECTIONS:
 for grid in [32,16]:
  d=pd.read_parquet(R/s/f'T{grid}_neighborhoods.parquet')
  for radius in [50,200]:
   for source in ['exclude_shared','include_shared']:
    v=d[(d.radius_um==radius)&d.source_rule.eq(source)&d.association_eligible].copy()
    coverage.append(dict(section=s,T_grid_um=grid,radius_um=radius,source_rule=source,comparable_T_regions=len(v),n_Epi_neighbor_occurrences=int(v.Epi_n.sum()),Treg_like_regions=int(v.evidence.eq('Treg_like_multigene').sum()),median_Epi_neighbors=float(v.Epi_n.median()),median_observed_fraction=float(v.observed_circle_fraction.median())))
    if grid==32 and source=='exclude_shared':spec=[('Treg_main',False,cfg['block_um'],'main_unadjusted'),('Treg_main',True,cfg['block_um'],'main'),('Treg_main',True,cfg['larger_block_um'],'larger_blocks'),('Treg_literature',True,cfg['block_um'],'literature'),('T_lineage',True,cfg['block_um'],'T_lineage_control'),('Tfh',True,cfg['block_um'],'Tfh_control'),('Activated_T',True,cfg['block_um'],'activated_control')]
    else:spec=[('Treg_main',True,cfg['block_um'],'scale16'if grid==16 else 'shared_source_inclusion')]
    for pred,adj,block,family in spec:
     rr=fit_effect(v,pred,adjusted=adj,block_um=block);results.append(dict(section=s,T_grid_um=grid,radius_um=radius,source_rule=source,predictor=pred,adjusted=adj,family=family,**rr));print(s,grid,radius,source,family,rr['n'],rr['effect'],rr['P'],flush=True)
    if grid==32 and source=='exclude_shared':
     rr=fit_effect(v,'Treg_main',block_um=cfg['block_um'],covars=['T_lineage','T_nCount','Epi_nCount','Epi_area_um2','observed_circle_fraction','T_background_balance','Tfh','Activated_T']);results.append(dict(section=s,T_grid_um=grid,radius_um=radius,source_rule=source,predictor='Treg_main',adjusted=True,family='competing_state_adjustment',**rr))
 # Reverse: Epi scoring regions are observations; each neighboring T region contributes once.
 pairs=pd.read_parquet(R/s/'region_pairs.parquet');t=pd.read_parquet(R/s/'T32_regions.parquet').set_index('region_id');e=pd.read_parquet(R/s/'Epi64_regions.parquet').set_index('region_id')
 for radius in [50,200]:
  pp=pairs[(pairs.T_grid_um==32)&(pairs.radius_um==radius)&~pairs.shared_original_bins];rows=[]
  for eid,v in pp.groupby('Epi_region_id'):
   q=t.loc[v.T_region_id];er=e.loc[eid];rows.append(dict(Epi_region_id=eid,x_um=er.x_um,y_um=er.y_um,component=er.component,Epi_MP10=er.MP10,Treg_main=q.Treg_main.mean(),T_lineage=q.T_lineage.mean(),T_nCount=q.nCount.mean(),Epi_nCount=er.nCount,Epi_area_um2=er.support_area_um2,T_n=len(q)))
  v=pd.DataFrame(rows)
  if len(v):v.to_parquet(R/s/f'reverse_r{radius}.parquet',index=False)
  rr=fit_effect(v,'Epi_MP10',ycol='Treg_main',block_um=cfg['block_um'],covars=['T_lineage','T_nCount','Epi_nCount','Epi_area_um2']);results.append(dict(section=s,T_grid_um=32,radius_um=radius,source_rule='exclude_shared',predictor='Epi_MP10',adjusted=True,family='reverse_consistency',**rr))
# Families never select positive slices; every section/radius, including NA, is retained.
out=bh(pd.DataFrame(results),['family','source_rule']);put(out,'association_all_results.tsv');put(out[out.family.eq('main')],'association_primary.tsv');put(pd.DataFrame(coverage),'neighborhood_coverage.tsv')
print('ASSOCIATIONS_COMPLETE',flush=True)
