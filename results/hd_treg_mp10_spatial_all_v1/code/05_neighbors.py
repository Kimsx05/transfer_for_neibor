from common import *
import sys
s=sys.argv[1];gate();cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());b=pd.read_parquet(R/s/'bin_evidence.parquet');e=pd.read_parquet(R/s/'Epi64_regions.parquet');e=e[e.measurement_eligible].reset_index(drop=True);em=pd.read_parquet(R/s/'Epi64_members.parquet');ei=pd.Series(np.arange(len(e)),index=e.region_id);em['ei']=em.region_id.map(ei);bin_e=em.dropna(subset='ei').set_index('bin_id').ei.astype(int);trees={c:cKDTree(v[['x_um','y_um']])for c,v in e.groupby('component')};eix={c:v.index.to_numpy()for c,v in e.groupby('component')};pool=b[b.eligible_pool];ptrees={c:cKDTree(v[['x_um','y_um']])for c,v in pool.groupby('component')};allpairs=[]
for kind in ['T32','T16']:
 t=pd.read_parquet(R/s/f'{kind}_regions.parquet');tm=pd.read_parquet(R/s/f'{kind}_members.parquet');tm['ei']=tm.bin_id.map(bin_e);shared=tm.dropna(subset='ei').groupby('region_id').ei.apply(lambda x:set(x.astype(int))).to_dict()
 out=[]
 for r in [50,200]:
  for c,v in t[t.measurement_eligible].groupby('component'):
   pos=v[['x_um','y_um']].to_numpy();obs=ptrees[c].query_ball_point(pos,r,return_length=True)*64/(np.pi*r*r);near=trees[c].query_ball_point(pos,r)if c in trees else [[]for _ in range(len(v))]
   for row,js,frac in zip(v.itertuples(),near,obs):
    ids=eix[c][js]if len(js) else np.array([],int);exclude=shared.get(row.region_id,set());independent=np.array([j for j in ids if j not in exclude],int)
    for src,kk in [('exclude_shared',independent),('include_shared',ids)]:
     ee=e.iloc[kk];n=len(ee);eligible=n>=1 and frac>=cfg['circle_observed_fraction_min'];out.append(dict(section=s,T_region_id=row.region_id,T_grid_um=int(kind[1:]),radius_um=r,source_rule=src,component=c,x_um=row.x_um,y_um=row.y_um,Treg_main=row.Treg_main,Treg_literature=row.Treg_literature,T_lineage=row.T_lineage,Tfh=row.Tfh,Activated_T=row.Activated_T,T_nCount=row.nCount,T_nFeature=row.nFeature,Treg_main_detected=row.Treg_main_detected,Treg_core_detected=row.Treg_core_detected,evidence=row.evidence,T_background_balance=(row.Epi_UMI-row.Fib_UMI)/(row.Epi_UMI+row.Fib_UMI+1),observed_circle_fraction=min(1.,frac),Epi_n=n,Epi_MP10=float(ee.MP10.mean())if n else np.nan,Epi_nCount=float(ee.nCount.mean())if n else np.nan,Epi_area_um2=int(ee.support_area_um2.sum())if n else 0,excluded_shared_n=len(ids)-len(independent),association_eligible=eligible,unresolved_reason='none'if eligible else 'no_independent_scoreable_Epi'if n==0 else 'insufficient_observed_input_circle',tissue_status=row.tissue_status))
    # One physical relation per T/Epi region pair and radius, not per-member score replication.
    for j in ids:allpairs.append((row.region_id,e.region_id.iloc[j],int(kind[1:]),r,bool(j in exclude),float(np.hypot(row.x_um-e.x_um.iloc[j],row.y_um-e.y_um.iloc[j]))))
 pd.DataFrame(out).to_parquet(R/s/f'{kind}_neighborhoods.parquet',index=False)
pd.DataFrame(allpairs,columns=['T_region_id','Epi_region_id','T_grid_um','radius_um','shared_original_bins','centroid_distance_um']).to_parquet(R/s/'region_pairs.parquet',index=False)
print('NEIGHBORS_DONE',s,len(allpairs),flush=True)
