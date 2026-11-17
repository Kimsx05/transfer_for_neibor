from common import *
from load_result import load_section,bin_membership
cfg=yaml.safe_load((R/'frozen_config.yaml').read_text());checks=[];coverage=[];rng=np.random.default_rng(123)
def check(name,cond,section='all',detail=''):
 assert cond,(section,name,detail);checks.append(dict(section=section,check=name,passed=bool(cond),detail=str(detail)))
# Artificial example: exact cumulative physical radius, component isolation and region source overlap.
xy=np.array([[0,0],[30,40],[51,0],[100,0]],float);near=cKDTree(xy).query_ball_point(xy[0],50);check('physical50_inclusive_not_knn',near==[0,1]);check('radius200_cumulative',set(near)<=set(cKDTree(xy).query_ball_point(xy[0],200)))
for s in SECTIONS:
 z=load_section(s,include_members=True,include_pairs=True);b=pd.read_parquet(R/s/'bin_evidence.parquet');info=z['inputs'];check('source_identity',Path(info['counts_path']).stat().st_size==info['counts_bytes']and Path(info['counts_path']).stat().st_mtime_ns==info['counts_mtime_ns'],s)
 check('unique_bin_IDs',b.bin_id.is_unique,s);check('18085_gene_background',info['input_genes']==18085,s)
 if s.startswith('MPBC'):check('unknown_original_tissue_not_filled',b.in_tissue_original.isna().all(),s)
 for kind in ['T32','T16','Epi64']:
  d=z[kind];m=z['members'][kind];check(kind+'_unique_members',m.bin_id.is_unique,s);check(kind+'_aligned_IDs',m.region_id.isin(d.region_id).all()and m.bin_id.isin(b.bin_id).all(),s)
  n=m.groupby('region_id').size().reindex(d.region_id).to_numpy();check(kind+'_membership_count',np.array_equal(n,d.n_bins),s)
  v=m.merge(b[['bin_id','x_um','y_um','nCount_RNA','in_tissue_original']],on='bin_id',validate='one_to_one');q=v.groupby('region_id').agg(x_um=('x_um','mean'),y_um=('y_um','mean'),nCount=('nCount_RNA','sum')).reindex(d.region_id);check(kind+'_centroid_and_counts',np.allclose(q[['x_um','y_um']],d[['x_um','y_um']])and np.array_equal(q.nCount,d.nCount),s)
  check(kind+'_no_known_outside_tissue',not(v.in_tissue_original==0).any(),s)
  audit=next(a for a in info['audit']if a['branch']==kind);check(kind+'_per_gene_conservation_record',audit['per_gene_conservation']and audit['included_UMI']+audit['excluded_UMI']==audit['input_UMI'],s)
  if kind=='Epi64':check('unscoreable_not_Low',not(d.L50&~d.measurement_eligible).any(),s)
  coverage.append(dict(section=s,branch=kind,all_regions=len(d),measurable=int(d.measurement_eligible.sum()),Treg_like=int(d.evidence.eq('Treg_like_multigene').sum()),observed_score_zeros=int((d.measurement_eligible&d['Treg_main' if kind.startswith('T')else 'MP10'].eq(0)).sum()),supplied_bins=len(b),mask_status='unknown'if s.startswith('MPBC')else 'original',included_UMI=audit['included_UMI'],excluded_UMI=audit['excluded_UMI']))
 e=z['Epi64'].set_index('region_id');pp=z['pairs'];check('unique_region_pairs',not pp.duplicated(['T_region_id','Epi_region_id','T_grid_um','radius_um']).any(),s)
 # Direct radius and raw-bin overlap on sampled pairs, plus full pair component/distance checks.
 for grid in [32,16]:
  t=z[f'T{grid}'].set_index('region_id');p=pp[pp.T_grid_um==grid];tv=t.loc[p.T_region_id];ev=e.loc[p.Epi_region_id];dist=np.hypot(tv.x_um.to_numpy()-ev.x_um.to_numpy(),tv.y_um.to_numpy()-ev.y_um.to_numpy());check(f'T{grid}_physical_pair_distance',np.all(dist<=p.radius_um.to_numpy()+1e-8)and np.allclose(dist,p.centroid_distance_um),s);check(f'T{grid}_component_isolation',np.array_equal(tv.component,ev.component),s)
  tm=z['members'][f'T{grid}'].groupby('region_id').bin_id.apply(set).to_dict();em=z['members']['Epi64'].groupby('region_id').bin_id.apply(set).to_dict()
  for row in p.sample(min(30,len(p)),random_state=123).itertuples():check(f'raw_overlap_{grid}_{row.Index}',bool(tm[row.T_region_id]&em[row.Epi_region_id])==row.shared_original_bins,s)
  d=pd.read_parquet(R/s/f'T{grid}_neighborhoods.parquet');check(f'T{grid}_not_replicated_stat_units',not d.duplicated(['T_region_id','radius_um','source_rule']).any(),s)
  for row in d[d.association_eligible].sample(min(15,int(d.association_eligible.sum())),random_state=123).itertuples():
   pairs=p[(p.T_region_id==row.T_region_id)&p.radius_um.eq(row.radius_um)];pairs=pairs[~pairs.shared_original_bins]if row.source_rule=='exclude_shared'else pairs;check(f'equal_region_mean_{grid}_{row.Index}',len(pairs)==row.Epi_n and np.isclose(e.loc[pairs.Epi_region_id,'MP10'].mean(),row.Epi_MP10),s)
  d50=d[d.radius_um==50].set_index(['T_region_id','source_rule']);d200=d[d.radius_um==200].set_index(['T_region_id','source_rule']);check(f'T{grid}_cumulative_neighbor_n',np.all(d200.loc[d50.index,'Epi_n'].to_numpy()>=d50.Epi_n.to_numpy()),s)
 check('Python_loader_membership',len(bin_membership(s))==len(b),s)
primary=pd.read_csv(R/'tables/association_primary.tsv',sep='\t');check('all_ten_primary_rows',len(primary)==10 and set(primary.section)==set(SECTIONS));check('source_disjoint_primary',primary.source_rule.eq('exclude_shared').all());check('NE_has_no_P',primary.loc[primary.status.str.startswith('descriptive'),'P'].isna().all())
started=json.loads((R/'models/association_started.json').read_text());check('parameter_freeze_precedes_association',cfg['frozen_at']<started['started']);put(pd.DataFrame(checks),'validation_checks.tsv');put(pd.DataFrame(coverage),'measurement_coverage.tsv');print('VALIDATED',len(checks),flush=True)
