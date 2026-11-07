from common import *
import json

def main():
 g=pd.read_csv(S/'tables/gene_sets.tsv',sep='\t');sets=g.groupby('gs_name').gene_symbol.agg(frozenset).to_dict();pan=g.groupby('gs_name').panel.first();mp=set('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split())
 groups={}
 for p in pan[pan.eq('main')].index:groups.setdefault(sets[p],[]).append(p)
 mapping=[];lookup={}
 for genes,names in groups.items():
  p='BUDDEN_CHOLESTEROL_SYNTHESIS_27' if 'BUDDEN_CHOLESTEROL_SYNTHESIS_27' in names else names[0];remain=genes-mp;matches=[a for a,s in sets.items() if s==remain and pan[a]=='sensitivity'];no=p if not genes&mp else next((a for a in matches if a.startswith(p.split('_27')[0])),matches[0] if matches else None)
  row=dict(canonical_program=p,full_program=p,noMP10_program=no,n_original=len(genes),n_shared=len(genes&mp),n_remaining=len(remain),shared_genes='|'.join(sorted(genes&mp)),aliases='|'.join(names),version='|'.join(sorted(g[g.gs_name.isin(names)].db_version.astype(str).unique())))
  mapping.append(row)
  for n in names:lookup[n]=(p,'full',n==p)
  if no!=p and no:lookup[no]=(p,'noMP10',True)
 lookup['PURAM_PEMT_TOP100']=('PURAM_PEMT_TOP30','Top100',True)
 # The accepted name can differ; identify by literal source list, never approximate members.
 for p in sets:
  if 'PURAM' in p and '100' in p:lookup[p]=('PURAM_PEMT_TOP30','Top100',True)
 mapping=pd.DataFrame(mapping).sort_values('canonical_program');save(mapping,'canonical_programs.tsv');save(M[META+['patient_id_status','P95_status','front_status','niche_status']], 'section_manifest.tsv')
 print('canonical',len(mapping),'overlap',mapping.n_shared.gt(0).sum(),'no overlap',mapping.n_shared.eq(0).sum(),flush=True)
 registry=[];pieces=[];coverage=[];missing=[];counter=0
 def ingest(d,module,statistic,value_col,method='UCell',range_type='original',candidate_col=None,n_col='n',n_other=None,se_col=None,blocks_col=None,lo_col=None,hi_col=None,program_col='program',source=''):
  nonlocal counter
  x=d.copy();x=x[x.section_id.isin(T.section_id)]
  if program_col and program_col in x:
   x[program_col]=x[program_col].astype(str);x=x[x[program_col].isin(lookup)];x=x[x[program_col].map(lambda p:lookup[p][2])];x['canonical_program']=x[program_col].map(lambda p:lookup[p][0]);x['signature']=x[program_col].map(lambda p:lookup[p][1]);x['actual_program']=x[program_col]
  else:x['canonical_program']='MP10' if module=='C4' else 'distance';x['signature']='not_applicable';x['actual_program']=x.canonical_program
  x['candidate_type']=x[candidate_col] if candidate_col else 'none'
  keys=['canonical_program','actual_program','signature','candidate_type']
  for vals,ss in x.groupby(keys,observed=True,dropna=False):
   counter+=1;eid=counter;cp,ap,sg,ca=vals;kind='r' if statistic in ['r','partial_r'] else 'identity'
   registry.append(dict(endpoint_id=eid,module=module,statistic=statistic,score_method=method,range_type=range_type,canonical_program=cp,actual_program=ap,signature=sg,candidate_type=ca,candidate_block=category(ca) if ca!='none' else 'none',kind=kind,core=cp in CORE,source=source))
   ss=ss.merge(T[META],on='section_id',suffixes=('','_manifest'),how='left',validate='m:1');nn=ss[n_col] if n_col in ss else pd.Series(np.nan,index=ss.index);oo=ss[n_other] if n_other and n_other in ss else nn;bb=ss[blocks_col] if blocks_col and blocks_col in ss else pd.Series(np.nan,index=ss.index)
   yy=pd.to_numeric(ss[value_col],errors='coerce');se=ss[se_col].to_numpy(float) if se_col and se_col in ss else np.full(len(ss),np.nan)
   if module=='C4' and statistic=='r' and lo_col:se=(np.arctanh(ss[hi_col].clip(-.999999,.999999))-np.arctanh(ss[lo_col].clip(-.999999,.999999))).to_numpy()/3.919928
   if module=='C2' and lo_col and hi_col:se=(ss[hi_col].to_numpy()-ss[lo_col].to_numpy())/3.919928
   out=pd.DataFrame(dict(endpoint_id=eid,section_id=ss.section_id.to_numpy(),dataset_id=ss.dataset_id.to_numpy(),analysis_unit_uid=ss.analysis_unit_uid.to_numpy(),patient_uid=ss.patient_uid.to_numpy(),value=yy.to_numpy(),n=nn.to_numpy(),n_other=oo.to_numpy(),blocks=bb.to_numpy(),spatial_se=se))
   if lo_col:out['source_ci_low']=ss[lo_col].to_numpy();out['source_ci_high']=ss[hi_col].to_numpy()
   if module=='C2':out['paired_r_P95']=ss.r_P95.to_numpy();out['paired_r_nonP95']=ss.r_nonP95.to_numpy()
   pieces.append(out)
 # Whole context table is a summary, not point-level parquet. Source-origin guards prevent new5/old80 duplication.
 c=pd.read_parquet(S/'tables/context85_effects.parquet');new=set(p for p in lookup if p.startswith('BUDDEN'))
 c=c[(c.origin.eq('reused_previous80')&~c.program.isin(new))|(c.origin.eq('new5_same_frozen_method')&c.program.isin(new))]
 for analysis in ['B1_R0','B1_R1','B1_R2','B3_Epi_distance','B3_tissue_distance']:
  for adj,st,col in [('raw','r','r'),('adjusted','partial_r','partial_r')]:
   q=c[c.analysis.eq(analysis)&c.adjustment.eq(adj)];assert not q.duplicated(['section_id','program']).any();ingest(q,analysis,st,col,range_type='0_1000' if analysis.startswith('B3') else 'graph_ring',blocks_col='blocks',source='context85_effects.parquet (guarded old80/new5 origins)')
 # Load aggregate tables in place, never point-level sources.
 p=pd.read_parquet(S/'tables/pearson_all.parquet');p=p[p.population.eq('P95')]
 for (method,ran),d in p.groupby(['score_method','range_type'],observed=True):
  if ran!='common_support_0_1000':d=d[d.program.map(lambda p:lookup.get(p,('',))[0]).isin(CORE)]
  ingest(d,'C1','r','pearson_r',method,ran,'candidate_type','n_spots',se_col='spatial_fisher_z_se',blocks_col='n_spatial_blocks',lo_col='spatial_ci_low',hi_col='spatial_ci_high',source='pearson_all.parquet')
  if ran=='common_support_0_1000':ingest(d,'C1','partial_r','partial_r_epi_depth',method,ran,'candidate_type','n_partial',blocks_col='n_spatial_blocks',source='pearson_all.parquet; no partial-r spatial SE')
 del p
 de=pd.read_parquet(S/'tables/delta_all.parquet')
 for (method,ran),d in de.groupby(['score_method','range_type'],observed=True):
  if ran!='common_support_0_1000':d=d[d.program.map(lambda p:lookup.get(p,('',))[0]).isin(CORE)]
  ingest(d,'C2','paired_delta_r','delta_r',method,ran,'candidate_type','n_P95','n_nonP95',blocks_col='n_blocks_P95',lo_col='spatial_ci_low',hi_col='spatial_ci_high',source='delta_all.parquet; same-draw joint spatial bootstrap')
 del de
 q=pd.read_csv(S/'tables/closeness_all.tsv.gz',sep='\t')
 q['zero_fraction_difference']=q.zero_fraction_P95-q.zero_fraction_nonP95
 for col in ['mean_difference_um','median_difference_um','matched_mean_difference_um','matched_median_pair_difference_um','zero_fraction_difference']:
  ingest(q,'C3',col,col,'distance_um' if 'um' in col else 'proportion','full_range','candidate_type','n_matched' if col.startswith('matched') else 'n_P95',program_col=None,lo_col='spatial_ci_low' if col=='mean_difference_um' else None,hi_col='spatial_ci_high' if col=='mean_difference_um' else None,source='closeness_all.tsv.gz; matched/median descriptive at section level')
 q=pd.read_csv(S/'tables/MP10_candidate_all.tsv.gz',sep='\t')
 for stat,col,n in [('r','pearson_r','n'),('partial_r','partial_r_epi_depth','n_partial')]:ingest(q,'C4',stat,col,'AUCell','full_range','candidate_type',n,blocks_col='n_blocks',lo_col='spatial_ci_low' if stat=='r' else None,hi_col='spatial_ci_high' if stat=='r' else None,program_col=None,source='MP10_candidate_all.tsv.gz')
 a=pd.concat(pieces,ignore_index=True);reg=pd.DataFrame(registry);save(reg,'endpoints.tsv.gz');del pieces
 a['record_id']=a.analysis_unit_uid.fillna('section_only::'+a.section_id);a['identity_level']=np.where(a.analysis_unit_uid.notna(),'analysis_unit','section_only');a['n']=pd.to_numeric(a.n,errors='coerce');a['value_transform']=a.value
 isr=a.endpoint_id.isin(reg.loc[reg.kind.eq('r'),'endpoint_id']);a.loc[isr,'value_transform']=transform(a.loc[isr,'value'],'r')
 a['variance']=a.spatial_se**2
 for col in ['paired_r_P95','paired_r_nonP95']:a[col]=np.arctanh(a[col].clip(-.999999,.999999))
 # Retain availability and source uncertainty; no source beta P is copied.
 unit=a.groupby(['endpoint_id','dataset_id','record_id','identity_level'],observed=True).agg(value_transform=('value_transform','mean'),n_sections=('section_id','nunique'),section_ids=('section_id',lambda s:'|'.join(s)),n_min=('n','min'),n_other_min=('n_other','min'),blocks_min=('blocks','min'),n_valid_sections=('value','count'),spatial_var_sum=('variance','sum'),n_spatial_se=('spatial_se','count'),paired_r_P95=('paired_r_P95','mean'),paired_r_nonP95=('paired_r_nonP95','mean'),source_ci_low=('source_ci_low','mean'),source_ci_high=('source_ci_high','mean')).reset_index()
 for col in ['paired_r_P95','paired_r_nonP95']:unit[col]=np.tanh(unit[col])
 unit['value']=unit.value_transform;mask=unit.endpoint_id.isin(reg.loc[reg.kind.eq('r'),'endpoint_id']);unit.loc[mask,'value']=np.tanh(unit.loc[mask,'value_transform']);unit['spatial_se']=np.where(unit.n_spatial_se.eq(unit.n_valid_sections)&unit.n_valid_sections.gt(0),np.sqrt(unit.spatial_var_sum)/unit.n_valid_sections,np.nan);unit=unit.drop(columns=['spatial_var_sum','value_transform']);save(unit,'unit_endpoints.tsv.gz')
 # Coverage includes endpoints unavailable in one or more registered sections.
 cov=a.merge(reg[['endpoint_id','module','statistic','score_method','range_type','signature','canonical_program','candidate_type']],on='endpoint_id',validate='m:1')
 for key,z in cov.groupby(['module','statistic','score_method','range_type','signature'],observed=True):
  v=z[np.isfinite(z.value)];coverage.append(dict(zip(['module','statistic','score_method','range_type','signature'],key))|dict(datasets='|'.join(sorted(v.dataset_id.unique())),n_datasets=v.dataset_id.nunique(),n_sections=v.section_id.nunique(),n_analysis_units=v.analysis_unit_uid.nunique(),n_programs=v.canonical_program.nunique(),n_candidates=v.candidate_type.nunique(),n_missing_rows=z.value.isna().sum(),missing_reason='see source status/reason; missing unit identity excluded from unit inference'))
 save(coverage,'input_coverage.tsv');save(cov[cov.value.isna()].groupby(['module','dataset_id','section_id']).size().reset_index(name='n_missing_endpoint_rows'),'endpoint_missing_coverage.tsv.gz')
 del a,cov
 # Formal R12 composition: keep denominators and labelled coverage; no replacement of unlabelled with zero.
 niche=pd.read_csv(S/'tables/niche_counts_all.tsv.gz',sep='\t');niche=niche[niche.section_id.isin(T.section_id)].copy();niche['record_id']=niche.analysis_unit_uid.fillna('section_only::'+niche.section_id)
 save(niche[['section_id','record_id','dataset_id','scope','population','niche','n_eligible','n_labelled','n_unlabelled','fraction_unlabelled','n_niche','fraction']],'niche_denominators.tsv.gz')
 summaries=[]
 for (scope,pop),z in niche.groupby(['scope','population']):
  for (ds,rid,k),ss in z.groupby(['dataset_id','record_id','niche']):
   summaries.append(dict(input_family='niche_'+scope+'_'+pop,record_id=rid,dataset_id=ds,feature=str(k),value=ss.loc[ss.n_labelled.gt(0),'fraction'].mean(),n_observed=ss.n_labelled.sum(),fraction_unlabelled=ss.fraction_unlabelled.mean(),transform='hellinger'))
 # Nine mutually exclusive major lineages recovered exactly by weighting existing disjoint region means by valid spot counts.
 comp=pd.read_csv(R/'mp10_paper_distance_full_v1/tables/region_composition.tsv',sep='\t');comp=comp[comp.axis.eq('lesion_location')&comp.section_id.isin(T.section_id)];lineages=['B','Endothelial','Epithelial','Fibroblast','Mast','Mural','Myeloid','Plasma','T_NK']
 for sec,z in comp.groupby('section_id'):
  meta=T.set_index('section_id').loc[sec];rid=meta.analysis_unit_uid if pd.notna(meta.analysis_unit_uid) else 'section_only::'+sec
  for l in lineages:
   good=z.n_valid.gt(0)&z[l+'_fraction'].notna();val=np.average(z.loc[good,l+'_fraction'],weights=z.loc[good,'n_valid']) if good.any() else np.nan;summaries.append(dict(input_family='cell_composition_whole_tissue',record_id=rid,dataset_id=meta.dataset_id,feature=l,value=val,n_observed=z.loc[good,'n_valid'].sum(),fraction_unlabelled=np.nan,transform='hellinger'))
 state=pd.read_csv(R/'mp10_paper_distance_full_v1/tables/region_programs.tsv',sep='\t');state=state[state.axis.eq('lesion_location')&state.region.eq('inside')&state.population.eq('epi_query')&state.section_id.isin(T.section_id)]
 if state.empty:
  print('state populations',pd.read_csv(R/'mp10_paper_distance_full_v1/tables/region_programs.tsv',sep='\t',usecols=['population']).population.unique(),flush=True)
 for rr in state.itertuples():
  meta=T.set_index('section_id').loc[rr.section_id];rid=meta.analysis_unit_uid if pd.notna(meta.analysis_unit_uid) else 'section_only::'+rr.section_id;cp=lookup.get(rr.program,(rr.program,))[0]
  summaries.append(dict(input_family='UCell_internal_Epi_median',record_id=rid,dataset_id=meta.dataset_id,feature=cp,value=rr.median,n_observed=rr.n_valid,fraction_unlabelled=np.nan,transform='standard'))
 su=pd.DataFrame(summaries).groupby(['input_family','record_id','dataset_id','feature','transform'],dropna=False).agg(value=('value','mean'),n_observed=('n_observed','sum'),fraction_unlabelled=('fraction_unlabelled','mean')).reset_index();save(su,'unit_background.tsv.gz')
 optional=[dict(input='14_Epi_MP_sample_activity',status='not_found_as_compatible_small_table',minimum_needed='section_id, frozen Epi population, MP1..MP14 mean/median/q25/q75/n; no full expression required'),dict(input='whole_tissue_major_lineage',status='available_exact_weighted_disjoint_region_means',minimum_needed='none; 9 major lineages, not 30-subtype composition'),dict(input='formal_R12_spatial_adjacency',status='no_compatible_whole_cohort_small_summary',minimum_needed='section_id, R12 pair adjacency counts, eligible edges, tissue-component-aware denominator'),dict(input='molecular_subtype_clinical',status='no_verified_crosswalk_linked_summary_used',minimum_needed='section/sample to verified subtype/pathology/clinical label and provenance'),dict(input='existing_dataset_combination_cache',status='named reanalysis/sweep directories not found; rebuilt from original effect summaries',minimum_needed='none')];save(optional,'optional_inputs.tsv')
 unitmeta=T[T.analysis_unit_uid.notna()].groupby(['analysis_unit_uid','dataset_id']).agg(section_ids=('section_id',lambda s:'|'.join(s)),n_Epi=('n_Epi','mean'),n_P95=('n_P95','mean'),niche_labelled=('n_niche_labelled_Epi','mean'),front_coverage=('fraction_front_exact_Epi','mean'),n_sections=('section_id','size'),patient_uid=('patient_uid','first')).reset_index().rename(columns={'analysis_unit_uid':'record_id'});unitmeta['niche_coverage']=unitmeta.niche_labelled/unitmeta.n_Epi.replace(0,np.nan);save(unitmeta,'unit_metadata.tsv')
 print('Prepared',len(reg),'endpoints',len(unit),'unit/section records',len(su),'background cells',flush=True)
if __name__=='__main__':main()
