from core import *
import json,shutil
CORE=['BUDDEN_TG_CATABOLISM_24','BUDDEN_CHOLESTEROL_SYNTHESIS_27','BUDDEN_S1P_RESPONSE_18','BUDDEN_TNF_EMT_RESPONSE_15','BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18','HALLMARK_MTORC1_SIGNALING__NO_MP10','HALLMARK_GLYCOLYSIS','HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_UNFOLDED_PROTEIN_RESPONSE','HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY','PURAM_PEMT_TOP30','FLORES_CRISPR_PRO_INVASION_10']
CORE=[p for p in CORE if p in P]
def adjust(a,p,col,keys):
 a[col]=a.groupby(keys,dropna=False,observed=True)[p].transform(lambda x:bh(x.to_numpy()))
 return a
def combine(kind):
 a=pd.concat([pd.read_parquet(O/'tables'/f'{s}_{kind}.parquet') for s in M.section_id],ignore_index=True)
 if kind=='pearson':p='p_pearson';prefix='q_pearson';keys=['score_method','program_family','population','range_type','source_context']
 elif kind=='delta':p='p_spatial_fisher_z';prefix='q_spatial';keys=['score_method','program_family','range_type','source_context']
 else:p='p_spatial' if kind=='closeness' else 'p_pearson';prefix='q_spatial' if kind=='closeness' else 'q_pearson';keys=['source_context']
 a['spatial_test_estimand']='P95_minus_nonP95_mean_distance_um' if kind=='closeness' else 'Pearson_or_Fisher_z_as_labelled'
 if kind=='pearson':
  a['partial_epi_depth_status']=np.where(a.n_partial.lt(5),'insufficient_complete_covariates_or_n',np.where(a.partial_r_epi_depth.notna(),'computed','degenerate_after_adjustment'));a['partial_front_status']=np.where(a.n_front_subset.lt(5),'insufficient_exact_front_subset',np.where(a.partial_r_front_subset_after.notna(),'computed_on_identical_front_subset','degenerate_after_front_adjustment'))
 if kind=='pearson':
  a.loc[a.n_partial.between(5,9),'partial_epi_depth_status']='sparse_n_lt10_descriptive_only';a.loc[a.n_front_subset.between(5,9),'partial_front_status']='sparse_n_lt10_front_adjustment_not_stable'
 adjust(a,p,prefix+'_within_section',['section_id']+keys);adjust(a,p,prefix+'_across_sections',keys)
 fam=a.groupby(keys,dropna=False,observed=True).agg(n_rows=(p,'size'),n_finite_tests=(p,'count'),n_sections=('section_id','nunique')).reset_index();save(fam,O/'tables'/f'{kind}_testing_families.tsv')
 for col in a.select_dtypes('object').columns:
  if a[col].nunique()<200:a[col]=a[col].astype('category')
 a.to_parquet(O/'tables'/f'{kind}_all.parquet',index=False)
 for sec,z in a.groupby('section_id',observed=True):z.to_parquet(O/'tables'/f'{sec}_{kind}_q.parquet',index=False)
 if kind in ['pearson','delta']:save(a[a.program.isin(CORE)],O/'tables'/f'{kind}_core.tsv.gz')
 else:save(a,O/'tables'/f'{kind}_all.tsv.gz')
 print(kind,len(a),flush=True);return a

def adjusted_delta(a,de):
 keys=['section_id','candidate_type','program','score_method','range_type'];vals=['partial_r_epi_depth','n_partial','partial_r_front_subset_before','partial_r_front_subset_after','n_front_subset'];h=a[a.population.eq('P95')][keys+vals];l=a[a.population.eq('nonP95')][keys+vals];z=h.merge(l,on=keys,suffixes=('_P95','_nonP95'),validate='1:1');z['delta_partial_r_epi_depth']=z.partial_r_epi_depth_P95-z.partial_r_epi_depth_nonP95;z['delta_partial_r_front_subset_before']=z.partial_r_front_subset_before_P95-z.partial_r_front_subset_before_nonP95;z['delta_partial_r_front_subset_after']=z.partial_r_front_subset_after_P95-z.partial_r_front_subset_after_nonP95;z['inference']='descriptive partial correlations; raw Pearson and spatial delta inference remain primary';z.to_parquet(O/'tables/partial_correlation_differences.parquet',index=False);save(z[z.program.isin(CORE)],O/'tables/partial_correlation_differences_core.tsv.gz')

def heterogeneity(a,de):
 z=a[a.score_method.eq('AddModuleScore_SCT_data')&a.range_type.eq('common_support_0_1000')&a.population.isin(['P95','nonP95'])&a.pearson_r.notna()].copy();z['fisher_z']=np.arctanh(z.pearson_r.clip(-.999999,.999999));keys=['source_context','dataset_id','candidate_type','program','program_family','population']
 unit=[]
 for level,idcol in [('analysis_unit','analysis_unit_uid'),('verified_patient','patient_uid')]:
  v=z[z[idcol].notna()].copy()
  if level=='verified_patient':v=v[v.eligible_for_patient_balanced_analysis.eq(True)]
  v['z_se_squared']=v.spatial_fisher_z_se**2
  u=v.groupby(keys+[idcol],observed=True,dropna=False).agg(n_sections=('section_id','size'),section_ids=('section_id',lambda x:'|'.join(x.astype(str))),fisher_z=('fisher_z','mean'),se_squared_sum=('z_se_squared','sum'),n_valid_se=('spatial_fisher_z_se','count'),n_spots_descriptive=('n_spots','sum'),minimum_n_P95=('n_spots','min'),distance_min_um=('distance_min_um','min'),distance_max_um=('distance_max_um','max')).reset_index().rename(columns={idcol:'unit_id'})
  u['level']=level;u['pearson_r']=np.tanh(u.fisher_z);u['spatial_z_se']=np.where(u.n_valid_se.eq(u.n_sections),np.sqrt(u.se_squared_sum)/u.n_sections,np.nan);unit.append(u)
 u=pd.concat(unit,ignore_index=True);u['spatial_ci_low']=np.tanh(u.fisher_z-1.96*u.spatial_z_se);u['spatial_ci_high']=np.tanh(u.fisher_z+1.96*u.spatial_z_se);save(u,O/'tables/unit_correlations.tsv.gz')
 ugroups={key:g for key,g in u[u.level.eq('analysis_unit')].groupby(keys,observed=True,dropna=False)}
 # Equal-weight descriptive effects; heterogeneity uses spatial SE, never n_spot-3.
 out=[]
 for k,g in z.groupby(keys,observed=True,dropna=False):
  rec=dict(zip(keys,k));ug=ugroups.get(k,u.iloc[:0])
  v=ug[np.isfinite(ug.spatial_z_se)&ug.spatial_z_se.gt(0)&ug.minimum_n_P95.ge(10)]
  qv,pv,i2=np.nan,np.nan,np.nan
  if len(v)>=3:
   w=1/v.spatial_z_se.to_numpy()**2;vals=v.fisher_z.to_numpy();weighted=np.sum(w*vals)/sum(w);qv=np.sum(w*(vals-weighted)**2);pv=stats.chi2.sf(qv,len(v)-1);i2=max(0,(qv-len(v)+1)/qv) if qv>0 else 0
  rec.update(n_sections=len(g),n_units=len(ug),n_units_with_spatial_interval=len(v),n_verified_patients=g.loc[g.eligible_for_patient_balanced_analysis.eq(True),'patient_uid'].nunique(),section_equal_r=np.tanh(g.fisher_z.mean()),unit_equal_r=np.tanh(ug.fisher_z.mean()) if len(ug) else np.nan,median_section_r=g.pearson_r.median(),n_negative=int((g.pearson_r<-.1).sum()),n_positive=int((g.pearson_r>.1).sum()),n_near_zero=int((g.pearson_r.abs()<=.1).sum()),r_min=g.pearson_r.min(),r_max=g.pearson_r.max(),Q_spatial=qv,p_heterogeneity=pv,I2_spatial=i2,range='common_support_0_1000',heterogeneity_status='estimated_conditional_spatial_SE' if len(v)>=3 else 'insufficient_units_with_spatial_SE');out.append(rec)
 h=pd.DataFrame(out);adjust(h,'p_heterogeneity','q_heterogeneity',['source_context','dataset_id','program_family','population']);save(h,O/'tables/dataset_heterogeneity.tsv.gz')
 summary=[]
 for k,g in h.groupby(['source_context','candidate_type','program','program_family','population'],observed=True):
  v=g.unit_equal_r.dropna();summary.append(dict(zip(['source_context','candidate_type','program','program_family','population'],k))|dict(n_datasets=len(g),n_datasets_with_units=len(v),dataset_equal_r=np.tanh(np.arctanh(v.clip(-.999999,.999999)).mean()) if len(v) else np.nan,n_sections=g.n_sections.sum(),n_units=g.n_units.sum(),n_verified_patients=g.n_verified_patients.sum(),n_positive_sections=g.n_positive.sum(),n_negative_sections=g.n_negative.sum(),n_near_zero_sections=g.n_near_zero.sum()))
 save(summary,O/'tables/dataset_equal_summary.tsv.gz')
 # Joint-bootstrap delta is descriptive at full range; prioritize formal common support.
 dz=de[de.score_method.eq('AddModuleScore_SCT_data')&de.range_type.eq('common_support')&de.delta_r.notna()].copy();dz['CI_excludes_zero']=(dz.spatial_ci_low>0)|(dz.spatial_ci_high<0)
 rep=dz.groupby(['dataset_id','candidate_type','program','program_family'],observed=True).agg(n_sections=('section_id','nunique'),n_units=('analysis_unit_uid','nunique'),median_delta_r=('delta_r','median'),n_negative=('delta_r',lambda x:(x<-.1).sum()),n_positive=('delta_r',lambda x:(x>.1).sum()),n_CI_supported=('CI_excludes_zero','sum'),minimum_n_P95=('n_P95','min')).reset_index();save(rep,O/'tables/delta_dataset_replication.tsv.gz')
 # Direct unit differences for the fixed core panel, preserving opposite directions.
 diff=[]
 for k,g in u[(u.level=='analysis_unit')&u.program.isin(CORE)&u.spatial_z_se.notna()].groupby(keys,observed=True):
  rows=list(g.itertuples())
  for i,aa in enumerate(rows):
   for bb in rows[i+1:]:
    se=np.hypot(aa.spatial_z_se,bb.spatial_z_se);df=aa.fisher_z-bb.fisher_z;diff.append(dict(zip(keys,k))|dict(unit_a=aa.unit_id,unit_b=bb.unit_id,r_a=aa.pearson_r,r_b=bb.pearson_r,difference_fisher_z=df,ci_low=df-1.96*se,ci_high=df+1.96*se,p_spatial_difference=2*stats.norm.sf(abs(df)/se) if se>0 else np.nan,interpretation='conditional_spatial_uncertainty_not_population_patient_validation'))
 dd=pd.DataFrame(diff)
 if len(dd):adjust(dd,'p_spatial_difference','q_spatial_difference',['dataset_id','program_family','population']);save(dd,O/'tables/core_unit_pair_differences.tsv.gz')

def optional_table(path):
 try:return pd.read_csv(path,sep='\t')
 except pd.errors.EmptyDataError:return pd.DataFrame()

def context_tables():
 from context_extra import extra
 for sec in M.section_id:extra(sec)
 old=pd.read_csv(OLD/'tables/effects.tsv.gz',sep='\t');old=old.merge(M[['section_id','dataset_id','analysis_unit_uid','patient_uid','source_context','eligible_for_patient_balanced_analysis']],on='section_id',validate='m:1');old['origin']='reused_previous80';new=pd.concat([pd.read_csv(O/'tables'/f'{s}_new5_context.tsv.gz',sep='\t') for s in M.section_id],ignore_index=True);new['origin']='new5_same_frozen_method';new=new.merge(M[['section_id','dataset_id','analysis_unit_uid','patient_uid','source_context','eligible_for_patient_balanced_analysis']],on='section_id',validate='m:1');a=pd.concat([old,new],ignore_index=True);a=a.rename(columns={c:c+'_legacy80' for c in a if c=='q' or c.startswith('q_')});a['program_family']=a.program.map(PAN).fillna('other_frozen_endpoint')
 a['test_p_context']=np.where(a.quadratic_beta.notna(),a.overall_p,a.p);adjust(a,'test_p_context','q_context85_within_section',['section_id','family','adjustment','program_family']);adjust(a,'test_p_context','q_context85_across_sections',['family','adjustment','program_family','source_context'])
 a.to_parquet(O/'tables/context85_effects.parquet',index=False);save(a[a.program.isin(CORE)],O/'tables/context_core_effects.tsv.gz')
 for name in ['niche_counts','front_only_niche','new5_bands','new5_curves','new5_reference_curves','new5_reference_rings']:
  frames=[optional_table(O/'tables'/f'{sec}_{name}.tsv.gz') for sec in M.section_id];v=pd.concat(frames,ignore_index=True)
  for sec,frame in zip(M.section_id,frames):
   if frame.empty and not len(frame.columns):save(pd.DataFrame(columns=v.columns),O/'tables'/f'{sec}_{name}.tsv.gz')
  normal_sections=set(M.loc[M.source_context.eq('adjacent_normal_source'),'section_id'])
  if name=='niche_counts':
   normal=v.section_id.isin(normal_sections);v=v[~normal|(v.scope.eq('all_Epi')&v.population.isin(['all','normal_Epi']))].copy();v.loc[v.section_id.isin(normal_sections),'population']='normal_Epi'
   for sec in normal_sections:save(v[v.section_id.eq(sec)],O/'tables'/f'{sec}_{name}.tsv.gz')
  if name=='front_only_niche':
   v.loc[v.section_id.isin(normal_sections),'status']='not_applicable_adjacent_normal'
  v=v.merge(M[['section_id','dataset_id','source_context','analysis_unit_uid','patient_uid']],on='section_id',validate='m:1')
  if name=='front_only_niche':adjust(v,'p_spatial','q_spatial_within_section',['section_id']);adjust(v,'p_spatial','q_spatial_across_sections',['source_context'])
  save(v,O/'tables'/f'{name}_all.tsv.gz')
 return a

def coverage():
 cov=[]
 for s in M.section_id:
  d=pd.read_parquet(O/'spots'/f'{s}.parquet');st=pd.read_csv(O/'tables'/f'{s}_status.tsv',sep='\t').iloc[0].to_dict();q=d['query'];tum=d.tumor_epi;hi=d.MP10_P95;f=tum&d.d_front.between(-200,0);nhi=int(hi.sum());st.update(n_tissue_valid=int(d.valid_tissue.sum()),n_front_exact_Epi=int((q&d.d_front.notna()).sum()),fraction_front_exact_Epi=float(d.loc[q,'d_front'].notna().mean()) if q.any() else np.nan,n_front_internal=int(f.sum()),n_front_P95=int((f&hi).sum()),MP10_P95_threshold=pd.read_csv(OLD/'tables'/f'{s}_coverage.tsv',sep='\t').iloc[0].P95_threshold,MP10_P95_actual_fraction=float(nhi/tum.sum()) if tum.any() else np.nan,n_niche_labelled_Epi=int((q&d.accepted_label.notna()).sum()),n_niche_unlabelled_Epi=int((q&d.accepted_label.isna()).sum()),n_pEMT_high_front=int(d.pEMT_high_front.sum()),SCT_AMS_status='complete_85' if st['n_AMS_valid']>0 else 'no_complete_score',P95_status='applicable' if nhi>=2 else 'single_P95_spot_not_inferential' if nhi==1 else 'not_applicable_normal' if s.startswith('GSE285715_N') else 'no_eligible_internal_Epi',front_status='available' if f.any() else 'no_eligible_front_Epi',niche_status='available' if (q&d.accepted_label.notna()).any() else 'no_labelled_eligible_Epi');cov.append(st)
 c=pd.DataFrame(cov).merge(M,on='section_id',suffixes=('','_manifest'),validate='1:1');save(c,O/'full_sections.tsv')
 for suffix,out in [('targets.tsv.gz','candidate_targets.tsv.gz'),('gene_coverage.tsv','program_gene_coverage.tsv.gz'),('software.tsv','software.tsv')]:
  root='scores' if suffix!='targets.tsv.gz' else 'tables';v=pd.concat([pd.read_csv(O/root/f'{s}_{suffix}',sep='\t') for s in M.section_id],ignore_index=True)
  if suffix=='gene_coverage.tsv':
   aliases={'CCN2':'CTGF','CAVIN3':'PRKCDBP','ACKR3':'CXCR7','GSDME':'DFNA5','P3H2':'LEPREL1'};sets=G.groupby('gs_name').gene_symbol.apply(list).to_dict();v['missing_SCT']=['|'.join(g for g in sets[r.program] if g not in str(r.used_symbols).split('|') and aliases.get(g,g) not in str(r.used_symbols).split('|')) for r in v.itertuples()]
   for sec,g in v.groupby('section_id'):save(g,O/'scores'/f'{sec}_gene_coverage.tsv')
  save(v,O/'tables'/out)
 save(pd.DataFrame([dict(section_id=s,spot_table=str(O/'spots'/f'{s}.parquet'),frozen_source=str(OLD/'spots'/f'{s}.parquet'),complete_abundance_source=str(F/'spots'/f'{s}.parquet'),scores=str(O/'scores'/f'{s}_SCT_AddModuleScore.tsv.gz')) for s in M.section_id]),O/'tables/input_index.tsv')
 return c
if __name__=='__main__':
 a=combine('pearson');de=combine('delta');combine('closeness');combine('MP10_candidate');adjusted_delta(a,de);heterogeneity(a,de);del a,de;context_tables();coverage();print('AGGREGATION COMPLETE',flush=True)
