from common import *
from scipy.sparse.csgraph import connected_components
eligibility=[];funnel=[];audit=[];strat=[];hot=[];des=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet');e=pd.read_parquet(R/s/'Epi_regions.parquet');tt=t[t.T_eligible].set_index('region_id');ee=e.set_index('region_id');pairs=pd.read_parquet(V/s/'region_pairs.parquet');pairs=pairs[pairs.T_grid_um.eq(32)&pairs.radius_um.eq(200)].copy();pairs['ti']=tt.index.get_indexer(pairs.T_region_id);pairs['ei']=ee.index.get_indexer(pairs.Epi_region_id);assert(pairs[['ti','ei']]>=0).all().all();dx=tt.x_um.to_numpy()[pairs.ti]-ee.x_um.to_numpy()[pairs.ei];dy=tt.y_um.to_numpy()[pairs.ti]-ee.y_um.to_numpy()[pairs.ei];assert np.allclose(np.hypot(dx,dy),pairs.centroid_distance_um);assert np.array_equal(tt.component.to_numpy()[pairs.ti],ee.component.to_numpy()[pairs.ei]);assert not pairs.duplicated(['T_region_id','Epi_region_id']).any()
 tm=pd.read_parquet(V/s/'T32_members.parquet');em=pd.read_parquet(V/s/'Epi64_members.parquet');both=tm.merge(em,on='bin_id',suffixes=('_T','_E'));shared=set(zip(both.region_id_T,both.region_id_E));actual=np.array([(a,b)in shared for a,b in zip(pairs.T_region_id,pairs.Epi_region_id)]);assert np.array_equal(actual,pairs.shared_original_bins)
 nn=pd.read_parquet(V/s/'T32_neighborhoods.parquet');nn=nn[nn.source_rule.eq('exclude_shared')];allout=[]
 for r in RADII:
  v=tt.reset_index().merge(nn[nn.radius_um.eq(r)][['T_region_id','observed_circle_fraction','Epi_n','Epi_MP10','Epi_nCount','Epi_area_um2','T_background_balance','association_eligible']],left_on='region_id',right_on='T_region_id',validate='one_to_one').drop(columns='T_region_id').rename(columns={'observed_circle_fraction':'coverage'});v['radius_um']=r;pool=pairs[pairs.centroid_distance_um<=r+1e-8];valid=pool[~pool.shared_original_bins];counts=valid.groupby('T_region_id').size();v['n_Epi_check']=v.region_id.map(counts).fillna(0).astype(int);assert np.array_equal(v.n_Epi_check,v.Epi_n);v['U_r']=v.coverage.ge(.5)&v.Epi_n.ge(1);assert v.U_r.equals(v.association_eligible)
  total=pool.groupby('T_region_id').size();v['raw_Epi_candidates']=v.region_id.map(total).fillna(0).astype(int);v['all_pairs_shared']=(v.raw_Epi_candidates>0)&(v.Epi_n==0)
  v['exclusion_reason']=np.select([v.coverage.lt(.5),v.raw_Epi_candidates.eq(0),v.all_pairs_shared],['coverage_below0.5','no_qualified_Epi_within_radius','all_pairs_shared_source'],default='eligible')
  for h in HIGHS:
   highpair=valid[ee[h].to_numpy()[valid.ei]];n=highpair.groupby('T_region_id').size();v['n_'+h]=v.region_id.map(n).fillna(0).astype(int);v['near_'+h]=pd.array(np.where(v.U_r,v['n_'+h].gt(0),None),dtype='boolean');v['fraction_'+h]=np.where(v.Epi_n>0,v['n_'+h]/v.Epi_n,np.nan)
  uv=v[v.U_r].copy();v['stratum']=-1;v.loc[uv.index,'stratum']=strata(uv,uv.Epi_n)
  for g,q in v[v.U_r].groupby('stratum'):
   for label in LABELS:strat.append(dict(section=s,radius_um=r,label=label,stratum=int(g),n=len(q),positive=int(q[label].sum()),mixed=0<q[label].sum()<len(q),component=int(q.component.iloc[0]),median_UMI=q.nCount.median(),median_T_UMI=q.T_UMI.median(),median_Epi_opportunity=q.Epi_n.median(),median_coverage=q.coverage.median()))
  eligibility.append(dict(section=s,radius_um=r,T_eligible=len(v),U_r=int(v.U_r.sum()),FOXP3_in_U=int(v.loc[v.U_r,'FOXP3_supported'].sum()),strict_in_U=int(v.loc[v.U_r,'strict'].sum()),known_tissue=s in SECTIONS[:3]))
  for reason,n in v.exclusion_reason.value_counts().items():funnel.append(dict(section=s,radius_um=r,reason=reason,n=int(n),denominator=len(v),count_type='first_exclusion'))
  for reason,mask in [('coverage_below0.5',v.coverage.lt(.5)),('no_qualified_Epi_within_radius',v.raw_Epi_candidates.eq(0)),('all_pairs_shared_source',v.all_pairs_shared)]:funnel.append(dict(section=s,radius_um=r,reason=reason,n=int(mask.sum()),denominator=len(v),count_type='rule_alone'))
  allout.append(v)
 pd.concat(allout,ignore_index=True).to_parquet(R/s/'neighborhoods.parquet',index=False)
 for label in LABELS:
  v=tt[tt[label]].reset_index();pp=pairlist(v,64);adj=sparse.coo_matrix((np.ones(len(pp)*2),(np.r_[pp[:,0],pp[:,1]],np.r_[pp[:,1],pp[:,0]])),shape=(len(v),len(v))).tocsr();nc,lab=connected_components(adj,directed=False)
  for i,row in v.iterrows():hot.append(dict(section=s,label=label,region_id=row.region_id,hotspot=int(lab[i]),x_um=row.x_um,y_um=row.y_um,definition='same-component positive centroids8connected by distance<=64um'))
 for group,mask in [('all_T_eligible',t.T_eligible),('FOXP3_supported',t.FOXP3_supported),('other_evaluable_T',t.T_eligible&~t.FOXP3_supported),('strict',t.strict),('FOXP3_only',t.FOXP3_supported&~t.strict),('strict_only',t.strict&~t.FOXP3_supported)]:
  q=t[mask]
  for col in ['nCount','nFeature','raw_FOXP3','raw_IL2RA','raw_CTLA4','raw_IKZF2','Treg_main','Tfh','Activated_T','T_lineage']:
   des.append(dict(section=s,group=group,feature=col,n=len(q),median=q[col].median(),q10=q[col].quantile(.1),q90=q[col].quantile(.9),positive_fraction=q[col].gt(0).mean()))
 audit.append(dict(section=s,pairs200=len(pairs),excluded_shared=int(actual.sum()),pair_unique=True,physical_distance_verified=True,same_component=True,source_mask_exact=True))
 print('GEOMETRY',s,flush=True)
put(pd.DataFrame(eligibility),'measurement_eligibility_by_section_radius.tsv');put(pd.DataFrame(funnel),'exclusion_funnel.tsv');put(pd.DataFrame(strat),'exchangeability_strata.tsv');put(pd.DataFrame(audit),'pair_geometry_validation.tsv');pd.DataFrame(hot).to_parquet(R/'tables/hotspot_membership.parquet',index=False);put(pd.DataFrame(des),'label_expression_distributions.tsv')
