from common import *
from analysis_stats import fit_effect,bh
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());gate();s='Cervilla';e=pd.read_parquet(OLD/'tables/Cervilla_epithelial_regions.parquet');e=e[e.score_eligible].reset_index(drop=True);t=pd.read_parquet(R/s/'T32_regions.parquet');t=t[t.measurement_eligible].reset_index(drop=True);cm=pd.read_csv(OLD/'tables/epithelial_region_membership.tsv.gz',sep='\t');cm=cm[cm.region_id.isin(e.region_id)].copy();cm['cell_label']=cm.cell_id.str.rsplit('-',n=1).str[-1].astype(int)
# Link source2um assignments to original global8um coordinates, then to new T regions.
q=pd.read_parquet(A/'segmentation_assignment_v2/Cervilla/bin_to_cell.parquet',columns=['array_row','array_col','cell_label']);q=q[q.cell_label.isin(cm.cell_label)].merge(cm[['cell_label','region_id']],on='cell_label',validate='many_to_one');q['array_row8']=q.array_row//4;q['array_col8']=q.array_col//4
b=pd.read_parquet(R/s/'bin_evidence.parquet');tm=pd.read_parquet(R/s/'T32_members.parquet');lookup=b[['bin_id','array_row','array_col']].merge(tm,on='bin_id');lookup=lookup.rename(columns={'region_id':'T_region_id','array_row':'array_row8','array_col':'array_col8'});over=q[['array_row8','array_col8','region_id']].drop_duplicates().merge(lookup[['array_row8','array_col8','T_region_id']],on=['array_row8','array_col8']);shared=over.groupby('T_region_id').region_id.apply(set).to_dict();over[['region_id','T_region_id']].drop_duplicates().to_parquet(R/s/'cell_region_shared_source.parquet',index=False)
# Component numbers from old cell processing and new supplied grids need actual coordinate reconciliation.
bt=cKDTree(b[['x_um','y_um']]);dist,j=bt.query(e[['x_um','y_um']]);e['component_current']=b.component.iloc[j].to_numpy();e['component_match_distance_um']=dist;e['coordinate_supported']=dist<=16
rows=[];matches=[]
for rad in [50,200]:
 for comp,v in t.groupby('component'):
  ee=e[(e.component_current==comp)&e.coordinate_supported];tree=cKDTree(ee[['x_um','y_um']]);near=tree.query_ball_point(v[['x_um','y_um']],rad)
  primary=pd.read_parquet(R/s/'T32_neighborhoods.parquet');primary=primary[(primary.radius_um==rad)&primary.source_rule.eq('exclude_shared')].set_index('T_region_id')
  for row,js in zip(v.itertuples(),near):
   z=ee.iloc[js];z=z[~z.region_id.isin(shared.get(row.region_id,set()))];base=primary.loc[row.region_id].to_dict()
   if len(z):base.update(T_region_id=row.region_id,Epi_MP10=z.MP10.mean(),Epi_nCount=z.nCount_RNA.mean(),Epi_area_um2=np.nan,Epi_n=len(z),association_eligible=base['observed_circle_fraction']>=.5);rows.append(base)
  # No centroid-based bin-cell one-to-one labeling; the same source-disjoint T centers define the comparison.
 v=pd.DataFrame([z for z in rows if z['radius_um']==rad]);v=v[v.association_eligible]if len(v)else v
 rr=fit_effect(v,'Treg_main',block_um=cfg['block_um'],covars=['T_lineage','T_nCount','Epi_nCount','observed_circle_fraction','T_background_balance']);matches.append(dict(section=s,radius_um=rad,method='existing_Epi_cell_region_MP10',**rr))
pd.DataFrame(rows).to_parquet(R/s/'cell_region_neighborhood_comparison.parquet',index=False);put(bh(pd.DataFrame(matches),['method']),'cell_region_association_comparison.tsv')
# Coverage and common-center score agreement, not forced physical one-to-one region matches.
n=pd.read_parquet(R/s/'T32_neighborhoods.parquet');z=pd.DataFrame(rows);out=[]
for rad in [50,200]:
 a=n[(n.radius_um==rad)&n.source_rule.eq('exclude_shared')&n.association_eligible];v=z[(z.radius_um==rad)&z.association_eligible];m=a.merge(v,on='T_region_id',suffixes=('_bin','_cell'));out.append(dict(section=s,radius_um=rad,bin_T_centers=len(a),cell_T_centers=len(v),common_T_centers=len(m),neighborhood_MP10_rho=spearmanr(m.Epi_MP10_bin,m.Epi_MP10_cell).statistic if len(m)>5 else np.nan,cell_score_regions=len(e),note='Cell-only score source retained; common T-center environment comparison, no one-to-one cell/bin match'))
put(pd.DataFrame(out),'cell_bin_common_coverage.tsv');print('CELL_COMPARE_DONE',flush=True)
