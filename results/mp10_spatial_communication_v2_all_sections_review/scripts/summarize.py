#!/usr/bin/env python3
"""Read-only cache aggregation; no inference, resampling of spots, or candidate selection."""
from summary_common import *
import argparse,time
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,ListedColormap
p=argparse.ArgumentParser();p.add_argument('--allow-partial',action='store_true');p.add_argument('--no-figures',action='store_true');p.add_argument('--common-lr-input');a=p.parse_args();start=time.time();dest=O/'comparison';dest.mkdir(exist_ok=True);fdir=O/'figures/summary';fdir.mkdir(parents=True,exist_ok=True)
methods=['COMMOT','SpatialCellChat'];terminal=['completed','reused','failed','blocked'];successful=['completed','reused']
I=read(O/'input_review/section_status.tsv');mother=read(O/'input_review/candidate_sections.tsv');sections=mother.section_id.tolist();assert len(sections)==63 and len(set(sections))==63
if (O/'logs/task_status.tsv').exists():J=read(O/'logs/task_status.tsv')
else:J=pd.DataFrame(columns=['method','section_id','status'])
assert not J.duplicated(['method','section_id']).any()
S=pd.DataFrame([(m,s) for m in methods for s in sections],columns=['method','section_id']).merge(I,on='section_id',how='left').merge(J,on=['method','section_id'],how='left',suffixes=('','_job'));S['status']=S.status.fillna('queued');S['postprocess_marker']=[(job_dirs(m,s)[0]/'POSTPROCESS_COMPLETED.json').exists() for m,s in zip(S.method,S.section_id)]
complete=bool(S.status.isin(terminal).all() and S.loc[S.status.isin(successful),'postprocess_marker'].all())
if not complete and not a.allow_partial:raise RuntimeError('Full summary refused: queued/running or successful jobs without completed postprocessing. Pass --allow-partial for explicitly provisional outputs.')
S['delivery_status']='final' if complete else 'partial';save(S,dest/'section_method_status.tsv')
cols=[];figs=[];geoms=[]
for r in S.itertuples():
 out,source,var=job_dirs(r.method,r.section_id)
 if (out/'POSTPROCESS_COMPLETED.json').exists():
  c=read(out/'epi_state_type_contrasts.tsv');c['run_status']=r.status;c['result_available']=True;cols.append(c)
  figs.append(read(out/'figure_index.tsv'));geoms.append(read(out/'batch_geometry_diagnostics.tsv'))
 else:
  # Explicit missing contrasts preserve all sections and types; never zero-fill failed or pending jobs.
  info=I[I.section_id==r.section_id];info=info.iloc[0].to_dict() if len(info) else {'section_id':r.section_id}
  inp=O/'input_review/sections'/r.section_id
  cross=read(inp/'spatial_label_crosswalk.tsv') if (inp/'spatial_label_crosswalk.tsv').exists() else read(O/'input_review/spatial_label_crosswalk.tsv')
  cross=cross.dropna(subset=['analysis_fine_label']).drop_duplicates('analysis_fine_label');pairs=[('major',x,x) for x in LINEAGES]+[('fine',row.major,row.analysis_fine_label) for row in cross.itertuples() if row.analysis_fine_label not in EPI and row.major!='Epithelial']
  rows=[]
  for level,lineage,node in pairs:
   for direction in ['epi_to_TME','TME_to_epi']:
    rows.append({**info,'method':r.method,'network_level':level,'lineage':lineage,'node':node,'direction':direction,'run_status':r.status,'result_available':False,'computed':False,'NE_reason':'job_'+r.status if r.status in ['failed','blocked'] else 'results_pending','strength_low':np.nan,'strength_high':np.nan,'delta_strength':np.nan,'log2_HL':np.nan,'total_flow_low':np.nan,'total_flow_high':np.nan,'delta_flow':np.nan,'effect_direction':'NE','zero_category':'NE','score_coverage_status':'complete20' if info.get('n_mp10_genes_measured')==20 else 'score_coverage_incomplete'})
  cols.append(pd.DataFrame(rows))
C=pd.concat(cols,ignore_index=True)
# Current audited annotations supersede historical run snapshots; numerical inference remains untouched.
for k in ['dataset','patient_uid','patient_id_status','analysis_unit_uid','analysis_unit_type','tissue_condition','treatment_status','annotation_revision','tissue_evidence_level','tissue_evidence','tissue_evidence_file','clinical_annotation_source_sheet','clinical_annotation_source_rows','sample_disease_subtype','clinical_disease_background','sample_pathologic_stage','specimen_site','tumor_involvement_verified','annotation_limitations']:
 if k in I:C[k]=C.section_id.map(I.set_index('section_id')[k])
C['computed']=C.computed.map(truth);C['result_available']=C.result_available.map(truth)
for k in ['dataset','tissue_condition','treatment_status','score_coverage_status','analysis_unit_type','patient_id_status']:C[k]=C[k].map(unknown)
C['identity_stratum']=np.where(C.analysis_unit_type.str.startswith('SOURCE_VERIFIED_PATIENT')&C.patient_uid.notna(),'verified_patient','unverified_record')
C['unit_uid']=np.where(C.identity_stratum.eq('verified_patient'),C.patient_uid,C.analysis_unit_uid.where(C.analysis_unit_uid.notna(),C.section_id))
C['source_conflict_recorded']=C.analysis_unit_type.eq('SOURCE_VERIFIED_PATIENT_WITH_RECORDED_SOURCE_CONFLICT')
assert not C.duplicated(['method','section_id','network_level','node','direction']).any();assert not C.loc[~C.computed,['strength_low','strength_high','log2_HL']].notna().any().any()
save(C,dest/'all_section_contrasts.tsv')
I['identity_stratum']=np.where(I.analysis_unit_type.fillna('').str.startswith('SOURCE_VERIFIED_PATIENT')&I.patient_uid.notna(),'verified_patient','unverified_record');I['unit_uid']=np.where(I.identity_stratum.eq('verified_patient'),I.patient_uid,I.analysis_unit_uid.where(I.analysis_unit_uid.notna(),I.section_id));save(I,dest/'section_input_identity_coverage.tsv')
uid=I[I.identity_stratum.eq('verified_patient')].groupby('patient_uid',dropna=False).agg(n_sections=('section_id','nunique'),n_datasets=('dataset','nunique'),datasets=('dataset',lambda x:';'.join(sorted(set(x))))).reset_index();save(uid,dest/'patient_cross_dataset_identity_audit.tsv');assert (uid.n_datasets<=1).all(),'A verified UID occurs across datasets; avoid duplicate weighting and review metadata first.'
K=['method','dataset','tissue_condition','treatment_status','score_coverage_status','identity_stratum','network_level','lineage','node','direction'];directions=['higher','lower','equal','NE'];zeros=['both_positive','low_zero_only','high_zero_only','both_zero','NE']
def effect_direction(g):
 s=set(g.effect_direction)-{'NE'}
 return 'NE' if not s else 'mixed' if 'higher' in s and 'lower' in s else 'higher' if 'higher' in s else 'lower' if 'lower' in s else 'equal'
def stats_record(g):
 vals=g.log2_HL.dropna().to_numpy();n=len(vals)
 return dict(n_finite_logratio=n,median_log2_HL=float(np.median(vals)) if n else np.nan,q25_log2_HL=float(np.quantile(vals,.25)) if n else np.nan,q75_log2_HL=float(np.quantile(vals,.75)) if n else np.nan,min_log2_HL=float(np.min(vals)) if n else np.nan,max_log2_HL=float(np.max(vals)) if n else np.nan)
def section_counts(g):
 d={'n_sections':g.section_id.nunique(),'n_evaluable_sections':int(g.computed.sum()),'n_available_sections':int(g.result_available.sum())}
 for e in directions:d['section_'+e+'_n']=int(g.effect_direction.eq(e).sum());d['section_'+e+'_fraction_of_all']=d['section_'+e+'_n']/len(g);d['section_'+e+'_fraction_evaluable']=d['section_'+e+'_n']/d['n_evaluable_sections'] if d['n_evaluable_sections'] and e!='NE' else np.nan
 for z in zeros:d['section_zero_'+z+'_n']=int(g.zero_category.eq(z).sum())
 return d
units=[]
for key,g in C.groupby(K+['unit_uid'],dropna=False,sort=False):
 z={**dict(zip(K+['unit_uid'],key)),**section_counts(g),**stats_record(g),'unit_direction':effect_direction(g),'section_ids':';'.join(sorted(g.section_id)),'source_conflict_recorded':bool(g.source_conflict_recorded.any()),'patient_id_statuses':';'.join(sorted(set(g.patient_id_status)))}
 z['log2_HL']=z['median_log2_HL'];z['unit_evaluable']=bool(g.computed.any())
 for f in ['strength_low','strength_high','delta_strength','total_flow_low','total_flow_high','delta_flow','n_LR','n_spots','target_effective_Q_mass','target_expected_abundance']:
  if f in g:z['median_section_'+f]=g[f].median()
 units.append(z)
U=pd.DataFrame(units);save(U,dest/'patient_and_unverified_record_summary.tsv')
datasets=[]
for key,g in U.groupby(K,dropna=False,sort=False):
 z=dict(zip(K,key));cs=C
 for k,v in z.items():cs=cs[cs[k].eq(v)]
 z.update(section_counts(cs));z.update(stats_record(g));z['n_units']=len(g);z['n_evaluable_units']=int(g.unit_evaluable.sum());z['n_verified_patients']=len(g) if z['identity_stratum']=='verified_patient' else 0;z['n_unverified_records']=len(g) if z['identity_stratum']=='unverified_record' else 0;z['n_evaluable_verified_patients']=int(g.unit_evaluable.sum()) if z['identity_stratum']=='verified_patient' else 0
 for d in ['higher','lower','equal','mixed','NE']:
  z['unit_'+d+'_n']=int(g.unit_direction.eq(d).sum());z['unit_'+d+'_fraction_of_all']=z['unit_'+d+'_n']/len(g);z['unit_'+d+'_fraction_evaluable']=z['unit_'+d+'_n']/z['n_evaluable_units'] if z['n_evaluable_units'] and d!='NE' else np.nan
 z['heterogeneity_opposite_unit_directions']=bool(z['unit_higher_n'] and z['unit_lower_n']);datasets.append(z)
D=pd.DataFrame(datasets);save(D,dest/'dataset_stratified_summary.tsv')
DK=[x for x in K if x!='dataset'];equal=[];lodo=[]
def equal_record(g):
 finite=g.median_log2_HL.dropna();out={'n_datasets_total':len(g),'n_datasets_finite':len(finite),'datasets':';'.join(sorted(g.dataset)),'dataset_equal_mean_log2_HL':finite.mean(),'median_dataset_log2_HL':finite.median(),'min_dataset_log2_HL':finite.min(),'max_dataset_log2_HL':finite.max(),'dataset_opposite_directions':bool((finite>0).any() and (finite<0).any()),'weighting':'one equal weight per dataset with finite patient/record median; identity strata separate'}
 for d in ['higher','lower','equal','mixed']:out['dataset_equal_'+d+'_fraction_evaluable']=g['unit_'+d+'_fraction_evaluable'].mean()
 out['dataset_equal_NE_fraction_of_all']=g.unit_NE_fraction_of_all.mean();return out
for key,g in D.groupby(DK,dropna=False,sort=False):
 base=dict(zip(DK,key));full=equal_record(g);equal.append({**base,**full})
 for ds in g.dataset:
  sub=equal_record(g[g.dataset.ne(ds)]);lodo.append({**base,'omitted_dataset':ds,**sub,'full_dataset_equal_mean_log2_HL':full['dataset_equal_mean_log2_HL'],'change_from_full':sub['dataset_equal_mean_log2_HL']-full['dataset_equal_mean_log2_HL']})
E=pd.DataFrame(equal);L=pd.DataFrame(lodo);save(E,dest/'dataset_equal_summary.tsv');save(L,dest/'leave_one_dataset_out.tsv')
# Exploratory tests follow the committed prespecified threshold and family, never selected candidates.
tests=[];rng=np.random.default_rng(20261027)
for key,g in U[U.identity_stratum.eq('verified_patient')&U.score_coverage_status.eq('complete20')].groupby(K,dropna=False,sort=True):
 z=dict(zip(K,key));v=g.log2_HL.dropna().to_numpy();n=len(v);z.update(n_verified_patients=len(g),n_finite_patients=n,test_status='tested' if n>=8 else 'NE_insufficient_finite_verified_patients',p_value=np.nan,median_log2_HL=np.median(v) if n else np.nan,bootstrap_ci_low=np.nan,bootstrap_ci_high=np.nan)
 if n>=8:
  z['p_value']=1. if not np.any(v) else float(stats.wilcoxon(v,alternative='two-sided',zero_method='wilcox',method='auto').pvalue);boot=np.median(rng.choice(v,size=(2000,n),replace=True),axis=1);z['bootstrap_ci_low'],z['bootstrap_ci_high']=np.quantile(boot,[.025,.975])
 tests.append(z)
T=pd.DataFrame(tests);T['FDR_BH']=np.nan;family=[x for x in K if x not in ['node','lineage','identity_stratum']]
for _,g in T.groupby(family,dropna=False):
 vals=g.p_value.dropna().sort_values();n=len(vals)
 if n:T.loc[vals.index,'FDR_BH']=np.minimum.accumulate((vals.to_numpy()*n/np.arange(1,n+1))[::-1])[::-1].clip(0,1)
save(T,dest/'exploratory_patient_tests.tsv')
# Matched cross-method assessment: units differ; compare directions and ratios only.
join=[x for x in ['section_id','dataset','tissue_condition','treatment_status','score_coverage_status','identity_stratum','network_level','lineage','node','direction']]
mc=['computed','strength_low','strength_high','log2_HL','effect_direction','zero_category'];A=C[C.method.eq('COMMOT')][join+mc].merge(C[C.method.eq('SpatialCellChat')][join+mc],on=join,suffixes=('_COMMOT','_SpatialCellChat'),validate='one_to_one');A['both_evaluable']=A.computed_COMMOT&A.computed_SpatialCellChat;A['direction_agrees']=np.where(A.both_evaluable,A.effect_direction_COMMOT.eq(A.effect_direction_SpatialCellChat),pd.NA);A['finite_both']=A.log2_HL_COMMOT.notna()&A.log2_HL_SpatialCellChat.notna();save(A,dest/'method_concordance_section.tsv')
ams=[]
for key,g in A.groupby([k for k in join if k!='section_id'],dropna=False,sort=False):
 both=g[g.both_evaluable];fin=g[g.finite_both];rho=stats.spearmanr(fin.log2_HL_COMMOT,fin.log2_HL_SpatialCellChat).statistic if len(fin)>=3 and fin.log2_HL_COMMOT.nunique()>1 and fin.log2_HL_SpatialCellChat.nunique()>1 else np.nan
 ams.append({**dict(zip([k for k in join if k!='section_id'],key)),'n_sections':len(g),'n_both_evaluable':len(both),'n_finite_both':len(fin),'n_direction_agrees':int(both.direction_agrees.map(truth).sum()),'direction_agreement_fraction':both.direction_agrees.map(truth).mean() if len(both) else np.nan,'spearman_logratio':rho,'interpretation':'paired within-section descriptive agreement; algorithms not independent biological replicates'})
save(ams,dest/'method_concordance_summary.tsv')
G=pd.concat(geoms,ignore_index=True) if geoms else pd.DataFrame();save(G,dest/'batch_geometry_diagnostics.tsv');diagn=[]
if len(G):
 diagnostic_cols=[x for x in G if x not in C.columns and x not in ['variant']];B=C.merge(G[['method','section_id']+diagnostic_cols],on=['method','section_id'],how='left',validate='many_to_one')
 for variable in ['library_counts','n_genes_detected','epi_fraction_frozen']:
  h,l=B['high_'+variable+'_median'],B['low_'+variable+'_median'];B['log2_HL_'+variable]=np.log2((h/l).where((h>0)&(l>0)))
 factors=['n_spots','n_measured_genes','n_LR','fraction_reachable_pairs_105um','fraction_reachable_pairs_200um','n_epi_high','n_epi_low','target_effective_Q_mass','log2_HL_library_counts','log2_HL_n_genes_detected','log2_HL_epi_fraction_frozen']
 save(B[[c for c in ['method','section_id']+K[1:]+['log2_HL','computed','target_expected_abundance']+factors if c in B]],dest/'batch_geometry_target_diagnostics.tsv')
 for key,g in B.groupby(K,dropna=False,sort=False):
  for f in factors:
   if f not in g:continue
   good=g[['log2_HL',f]].replace([np.inf,-np.inf],np.nan).dropna();rho=stats.spearmanr(good.log2_HL,good[f]).statistic if len(good)>=4 and good.log2_HL.nunique()>1 and good[f].nunique()>1 else np.nan
   diagn.append({**dict(zip(K,key)),'diagnostic':f,'n_finite_sections':len(good),'spearman_rho':rho,'p_value':np.nan,'interpretation':'descriptive quality/composition/geometry association; not a batch-corrected effect or causal estimate'})
save(diagn,dest/'batch_geometry_descriptive_correlations.tsv')
# Preserve native-cache sensitivity from parent; never regenerate native inference here.
common=Path(a.common_lr_input) if a.common_lr_input else dest/'common_LR_sensitivity.tsv';common_ready=common.exists();sensitivity=[]
if not common_ready:
 for stale in ['common_LR_vs_primary.tsv','common_LR_sensitivity_summary.tsv']:
  if (dest/stale).exists():(dest/stale).unlink()
if common_ready:
 R=read(common);rkeys=['method','section_id','network_level','lineage','node','direction'];primary_cols=['log2_HL','computed','tissue_condition','treatment_status','score_coverage_status','dataset','identity_stratum'];R=R.drop(columns=[k for k in primary_cols if k in R]);R=R.merge(C[rkeys+primary_cols],on=rkeys,how='left',validate='one_to_one');R['delta_log2_HL_common_vs_primary']=R.log2_HL_common-R.log2_HL;R['finite_both']=R.log2_HL.notna()&R.log2_HL_common.notna();R['sign_agrees']=np.where(R.finite_both,np.sign(R.log2_HL).eq(np.sign(R.log2_HL_common)),pd.NA);save(R,dest/'common_LR_vs_primary.tsv')
 for key,g in R.groupby(K,dropna=False,sort=False):
  f=g[g.finite_both];sensitivity.append({**dict(zip(K,key)),'n_sections':len(g),'n_finite_both':len(f),'median_logratio_change':f.delta_log2_HL_common_vs_primary.median(),'fraction_sign_agrees':f.sign_agrees.map(truth).mean() if len(f) else np.nan,'median_retained_fraction':g.retained_fraction.median()})
 save(sensitivity,dest/'common_LR_sensitivity_summary.tsv')
# Scientific static figures. Every major/fine lineage is displayed, regardless of direction or p-value.
F=[]
def emit(fig,stem,kind,source,**meta):
 stem=fdir/stem
 if not complete and fig._suptitle:fig._suptitle.set_text('[PARTIAL] '+fig._suptitle.get_text())
 fig.savefig(str(stem)+'.pdf',transparent=True,bbox_inches='tight');fig.savefig(str(stem)+'.png',transparent=True,dpi=300,bbox_inches='tight');plt.close(fig);F.append({'figure_type':kind,'source_data':str(source),'pdf':str(stem)+'.pdf','png':str(stem)+'.png','relative_pdf':str(stem.relative_to(O))+'.pdf','relative_png':str(stem.relative_to(O))+'.png',**meta})
def cmap_values(ax,mat,title,cmap='RdBu_r',lim=None):
 vals=mat.to_numpy(dtype=float);finite=np.abs(vals[np.isfinite(vals)]);lim=lim or (max(float(np.quantile(finite,.98)),.25) if len(finite) else 1.);cm=plt.get_cmap(cmap).copy();cm.set_bad('#e5e5e5');im=ax.imshow(vals,aspect='auto',cmap=cm,vmin=-lim,vmax=lim);ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=65,ha='right',fontsize=7);ax.set_yticks(range(len(mat)),mat.index,fontsize=7);ax.set_title(title,color='black');return im
if not a.no_figures:
 # Complete mother-list coverage, including failed, blocked and pending jobs.
 cov=S.pivot(index='section_id',columns='method',values='status').reindex(sections);codes={'queued':0,'running':1,'completed':2,'reused':3,'failed':4,'blocked':5};v=cov.apply(lambda s:s.map(codes));fig,ax=plt.subplots(figsize=(6,16));im=ax.imshow(v,cmap=ListedColormap(['#eeeeee','#fee08b','#66c2a5','#3288bd','#d53e4f','#9e0142']),vmin=-.5,vmax=5.5,aspect='auto');ax.set_yticks(range(63),sections,fontsize=7);ax.set_xticks([0,1],cov.columns);fig.colorbar(im,ax=ax,ticks=list(codes.values())).ax.set_yticklabels(list(codes));ax.set_title(('FINAL' if complete else 'PARTIAL')+' | all 63 sections, two algorithms');fig.tight_layout();emit(fig,'all_sections_method_coverage','coverage',dest/'section_method_status.tsv',method='both',network_level='all',lineage='all')
 groups=[('major','all_major')]+[('fine',lin) for lin in LINEAGES]
 for level,lineage in groups:
  group=C[(C.network_level==level)&(C.lineage.eq(lineage) if level=='fine' else True)]
  names=list(group.node.drop_duplicates())
  for method in methods:
   gm=group[group.method.eq(method)]
   fig,axs=plt.subplots(1,2,figsize=(max(12,len(names)*1.4),max(12,len(sections)*.22+2)),squeeze=False)
   # Rows ordered by dataset and tissue; incomplete-score and identity uncertainty explicitly tagged.
   ordtab=gm.drop_duplicates('section_id').sort_values(['dataset','tissue_condition','treatment_status','score_coverage_status','section_id']);order=ordtab.section_id.tolist();labels=[f'{r.section_id} | {r.tissue_condition}'+(' | score<20' if r.score_coverage_status!='complete20' else '')+(' | unverified' if r.identity_stratum=='unverified_record' else '') for r in ordtab.itertuples()]
   finite=np.abs(gm.log2_HL.dropna().to_numpy());lim=max(np.quantile(finite,.98),.25) if len(finite) else 1
   for ax,d in zip(axs[0],['epi_to_TME','TME_to_epi']):
    mat=gm[gm.direction.eq(d)].pivot(index='section_id',columns='node',values='log2_HL').reindex(index=order,columns=names);mat.index=labels;im=cmap_values(ax,mat,d,lim=lim)
   axs[0,1].set_yticklabels([]);fig.suptitle(f'{method} | {level}: {lineage}\nGray: NE or zero category; exact categories retained in tables. Clip color at 98th percentile.',fontsize=12);fig.subplots_adjust(top=.90,bottom=.11,left=.34,right=.86,wspace=.12);cax=fig.add_axes([.89,.27,.018,.45]);fig.colorbar(im,cax=cax,label='within-section log2(H/L), both positive only');emit(fig,safe(f'{method}_{level}_{lineage}_section_effects'),'section_effect_heatmap',dest/'all_section_contrasts.tsv',method=method,network_level=level,lineage=lineage)
   # Separate tissue, score and identity panels prevent unknown/normal/tumor pooling.
   gd=D[(D.method==method)&(D.network_level==level)&(D.lineage.eq(lineage) if level=='fine' else True)]
   for sk,ds in gd.groupby(['tissue_condition','treatment_status','score_coverage_status','identity_stratum'],dropna=False,sort=False):
    stratum=' | '.join(sk);fig,axs=plt.subplots(2,2,figsize=(max(11,len(names)*1.2),max(7,ds.dataset.nunique()*.8+4)),squeeze=False)
    for ii,d in enumerate(['epi_to_TME','TME_to_epi']):
     z=ds[ds.direction==d];mat=z.pivot(index='dataset',columns='node',values='median_log2_HL').reindex(columns=names);im=cmap_values(axs[ii,0],mat,d+' | patient/record median log2(H/L)');fig.colorbar(im,ax=axs[ii,0],shrink=.5)
     mat=z.pivot(index='dataset',columns='node',values='n_evaluable_units').reindex(columns=names);ax=axs[ii,1];im=ax.imshow(mat.to_numpy(float),aspect='auto',cmap='Blues',vmin=0);ax.set_xticks(range(len(names)),names,rotation=65,ha='right',fontsize=7);ax.set_yticks(range(len(mat)),mat.index,fontsize=7);ax.set_title(d+' | evaluable patient/record units');fig.colorbar(im,ax=ax,shrink=.5)
     for y in range(len(mat)):
      for x in range(len(names)):ax.text(x,y,str(int(mat.iloc[y,x])) if pd.notna(mat.iloc[y,x]) else 'NE',ha='center',va='center',fontsize=7,color='black')
    fig.suptitle(f'{method} | {level}: {lineage}\n{stratum}',fontsize=12);fig.tight_layout(rect=[0,0,1,.93]);emit(fig,safe(f'{method}_{level}_{lineage}_{stratum}_dataset'),'dataset_effect_and_coverage',dest/'dataset_stratified_summary.tsv',method=method,network_level=level,lineage=lineage,stratum=stratum)
   # Unit distributions retain within-dataset heterogeneity and zero/NE frequencies alongside finite effects.
   gu=U[(U.method==method)&(U.network_level==level)&(U.lineage.eq(lineage) if level=='fine' else True)]
   for sk,g in gu.groupby(['tissue_condition','treatment_status','score_coverage_status','identity_stratum'],dropna=False,sort=False):
    stratum=' | '.join(sk);fig,axs=plt.subplots(1,2,figsize=(14,max(5,len(names)*.5+2)))
    dsorder=sorted(g.dataset.unique());colors=plt.get_cmap('tab10')(np.linspace(0,1,max(len(dsorder),2)))
    for ax,d in zip(axs,['epi_to_TME','TME_to_epi']):
     for j,ds in enumerate(dsorder):
      z=g[(g.direction==d)&(g.dataset==ds)];offset=(j-(len(dsorder)-1)/2)*.075
      for ni,n in enumerate(names):
       vv=z[z.node==n].log2_HL.dropna().to_numpy()
       if len(vv):ax.scatter(vv,np.full(len(vv),ni+offset),s=9,alpha=.55,color=colors[j]);ax.scatter(np.median(vv),ni+offset,s=32,marker='|',color=colors[j])
      ax.scatter([],[],color=colors[j],s=15,label=ds)
     ax.axvline(0,c='black',lw=.7);ax.set_yticks(range(len(names)),names);ax.invert_yaxis();ax.set_title(d);ax.set_xlabel('patient/record median within-section log2(H/L)');ax.grid(axis='x',alpha=.15)
    axs[1].legend(fontsize=7,loc='upper left',bbox_to_anchor=(1,1));fig.suptitle(f'{method} | {level}: {lineage}\n{stratum}; finite values only, zero/NE counts in tables');fig.tight_layout(rect=[0,0,1,.92]);emit(fig,safe(f'{method}_{level}_{lineage}_{stratum}_distribution'),'unit_effect_distribution',dest/'patient_and_unverified_record_summary.tsv',method=method,network_level=level,lineage=lineage,stratum=stratum)
  ag=A[(A.network_level==level)&(A.lineage.eq(lineage) if level=='fine' else True)];fig,axs=plt.subplots(1,2,figsize=(12,max(6,len(names)*.45+2)))
  for ax,d in zip(axs,['epi_to_TME','TME_to_epi']):
   z=ag[(ag.direction==d)&ag.finite_both]
   for n in names:
    q=z[z.node==n];ax.scatter(q.log2_HL_COMMOT,q.log2_HL_SpatialCellChat,s=18,label=n,alpha=.65)
   ax.axhline(0,color='#777',lw=.6);ax.axvline(0,color='#777',lw=.6);ax.set_xlabel('COMMOT log2(H/L)');ax.set_ylabel('SpatialCellChat log2(H/L)');ax.set_title(d);ax.grid(alpha=.15)
  axs[1].legend(fontsize=7,bbox_to_anchor=(1,1),loc='upper left');fig.suptitle(f'{level}: {lineage} | matched section method concordance\nDescriptive pairs; tissue, identity and score strata retained in tables. Distinct native strength units.');fig.tight_layout(rect=[0,0,1,.9]);emit(fig,safe(f'{level}_{lineage}_method_concordance'),'method_concordance',dest/'method_concordance_section.tsv',method='both',network_level=level,lineage=lineage)
 if len(G):
  fig,axs=plt.subplots(2,3,figsize=(14,8));first=G.drop_duplicates('section_id')
  for ax,x,y in zip(axs.ravel(),['n_spots','n_spots','n_epi_high','low_library_counts_median','low_n_genes_detected_median','low_epi_fraction_frozen_median'],['fraction_reachable_pairs_105um','fraction_reachable_pairs_200um','n_epi_low','high_library_counts_median','high_n_genes_detected_median','high_epi_fraction_frozen_median']):
   for ds,g in first.groupby('dataset'):ax.scatter(g[x],g[y],s=17,label=ds,alpha=.7)
   ax.set_xlabel(x);ax.set_ylabel(y);ax.grid(alpha=.15)
  axs[0,0].legend(fontsize=6);fig.suptitle('Geometry, H/L support and local RNA/composition diagnostics; no batch correction');fig.tight_layout(rect=[0,0,1,.94]);emit(fig,'geometry_and_quality_diagnostics','batch_geometry',dest/'batch_geometry_diagnostics.tsv',method='both',network_level='all',lineage='all')
 if len(L):
  z=L[L.network_level.eq('major')];strata=['method','tissue_condition','treatment_status','score_coverage_status','identity_stratum']
  for sk,g in z.groupby(strata,dropna=False,sort=False):
   fig,axs=plt.subplots(1,2,figsize=(14,6));ss=' | '.join(sk)
   for ax,d in zip(axs,['epi_to_TME','TME_to_epi']):
    q=g[g.direction==d]
    for i,n in enumerate(LINEAGES):
     r=q[q.node==n];ax.scatter(r.dataset_equal_mean_log2_HL,np.full(len(r),i),s=18,alpha=.6,label='LODO' if i==0 else None)
     if len(r):ax.scatter(r.full_dataset_equal_mean_log2_HL.iloc[0],i,marker='D',s=35,color='black',label='All datasets' if i==0 else None)
    ax.set_yticks(range(8),LINEAGES);ax.axvline(0,color='#777',lw=.6);ax.invert_yaxis();ax.set_title(d);ax.set_xlabel('equal-weight mean of dataset medians, log2(H/L)');ax.legend(fontsize=8)
   fig.suptitle('Leave-one-dataset-out | '+ss);fig.tight_layout(rect=[0,0,1,.93]);emit(fig,safe('LODO_'+ss),'leave_one_dataset_out',dest/'leave_one_dataset_out.tsv',method=sk[0],network_level='major',lineage='all',stratum=ss)
if F:save(F,fdir/'figure_index.tsv')
if (fdir/'figure_index.tsv').exists():figs.append(read(fdir/'figure_index.tsv'))
FI=pd.concat(figs,ignore_index=True) if figs else pd.DataFrame();save(FI,O/'figure_index.tsv')
if len(FI):assert FI.pdf.map(lambda p:Path(p).exists()).all() and FI.png.map(lambda p:Path(p).exists()).all()
summary={'delivery_status':'final' if complete else 'partial','n_expected_sections':63,'n_expected_jobs':126,'n_prepared_sections':len(I),'n_success_jobs':int(S.status.isin(successful).sum()),'n_terminal_jobs':int(S.status.isin(terminal).sum()),'n_postprocessed_jobs':int(S.postprocess_marker.sum()),'n_failed_or_blocked_jobs':int(S.status.isin(['failed','blocked']).sum()),'n_evaluable_section_method_type_direction':int(C.computed.sum()),'n_figures':len(FI),'n_evaluable_section_method_jobs':int(C[C.computed][['method','section_id']].drop_duplicates().shape[0]),'n_exploratory_tests':int(T.p_value.notna().sum()),'common_LR_sensitivity_available':common_ready,'common_LR_input':str(common),'common_LR_is_QA_override':bool(a.common_lr_input),'elapsed_seconds':time.time()-start,'completion_rule':'126 terminal jobs; all completed/reused jobs have POSTPROCESS_COMPLETED; failures retain specific reasons and recovery entry'}
atom_json(summary,dest/'SUMMARY_STATUS.json');print(json.dumps(summary,ensure_ascii=False),flush=True)
