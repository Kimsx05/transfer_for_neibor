from common import *
import sys, re, warnings
from collections import defaultdict
J=Path('qc/jobs')
def job(j,part):
 f=J/f'{j}_{part}.tsv.gz'
 return rd(f) if f.exists() else pd.DataFrame()
K=['SREBF1','SREBF2']; endpoints=targets[2:]
meta=rd('cache/metadata.tsv.gz');mainmeta=meta[meta.main]; elig=rd('tables/eligible_metadata.tsv.gz')
panel=rd('tables/model_panel.tsv').set_index('gene');rawcov=rd('qc/raw_coverage.tsv')
members=rd('tables/program_members.tsv');progs=programs();controls=job('main','controls')
sign=lambda x:'up' if x>0 else 'down' if x<0 else None
def sample(j):
 if j=='alt':return meta[meta.alt]
 if j=='pilot':return meta[meta.pilot]
 if j.startswith('state_'):return mainmeta[mainmeta.state==j[6:]]
 if j.startswith('drop_'):return mainmeta[mainmeta.dataset!=j[5:]]
 if j.startswith('sub'):return mainmeta[mainmeta.patient_id.isin(job(j,'membership').patient_id)]
 return mainmeta
def design(j):
 s=sample(j)
 return dict(analysis_id=j,dataset_id='ALL' if not j.startswith('drop_') else 'without_'+j[5:],method='scTenifoldKnk',method_version='2.0.0',state=j[6:] if j.startswith('state_') else 'All assigned clusters',
 n_patients=s.patient_id.nunique(),n_cells=len(s),metric_type='heat_distance',effect_scale='WT heat-kernel regulatory distance (t=10)',uncertainty_type=None,uncertainty_low=None,uncertainty_high=None,p_value=None,q_value=None,prior_id='none_de_novo',config_id='epi-v1-seed123/scTenifold-v2',direction_available=False,predicted_direction=None,assessment_status='NE',ne_reason=None,effect_value=None,direction_score=None)
def gene_result(base,k,g,eff,cov):
 row={**base,'perturbed_gene':k,'endpoint_gene':g}
 if cov.empty:reason='state cells <80 or patients <4'
 elif k not in cov.index or g not in cov.index:reason='not in frozen modeling panel or absent feature'
 elif not cov.loc[k,'in_network'] or not cov.loc[g,'in_network']:reason='local detection/patient coverage or variance insufficient'
 elif cov.loc[k,'outdegree']<=0:reason='no WT outgoing edge for KO'
 elif (k,g) not in eff.index:reason='KO output unavailable'
 elif not eff.loc[(k,g),'reachable']:reason='no fitted directed WT path'
 else:reason=None
 if reason:return {**row,'ne_reason':reason}
 x=eff.loc[(k,g)];score=x.directionScore
 return {**row,'assessment_status':'EVALUABLE','effect_value':x.distance,'rank':x['rank'],'rank_denominator':int(cov.in_network.sum())-1,'predicted_direction':sign(score),'direction_available':sign(score) is not None,'direction_score':score,'depth_score':x.depth_direction_score,'pc1_score':x.pc1_direction_score,'p_value':x['p.value'],'q_value':x['p.adj']}
rows=[]; robust=[];coverage=[]; psummary=[]
frozen=pd.read_csv(P/'resources/mp10_regulatory_priors/collectri.tsv',sep='\t')
mainjobs=['main']+['state_'+x for x in sorted(mainmeta.state.unique())]
other=['alt']+['drop_'+x for x in sorted(mainmeta.dataset.unique())]+['sub%02d'%i for i in range(1,21)]
for j in mainjobs+other:
 base0=design(j);eff=job(j,'effects');cov=job(j,'coverage')
 if not cov.empty:cov=cov.set_index('gene')
 if not eff.empty:eff=eff.set_index(['perturbed_gene','gene'])
 if j in other:
  for k in K:robust.append(gene_result(base0,k,'NECTIN2',eff,cov))
  continue
 ks=K+sorted(set(controls.gene)) if j=='main' else K
 for k in ks:
  for g in endpoints if k in K else ['NECTIN2']:rows.append(gene_result(base0,k,g,eff,cov))
  if k not in K:continue
  for program,full in progs.items():
   for mask in ['Full','M20']:
    gs=full if mask=='Full' else full-progs['MP10']
    allowed=gs-{k};ok=[g for g in sorted(allowed) if not cov.empty and g in cov.index and cov.loc[g,'in_network'] and (k,g) in eff.index and eff.loc[(k,g),'reachable']] if not eff.empty else []
    values=[float(eff.loc[(k,g),'distance']) for g in ok]
    evaluable=len(ok)>=3 and len(ok)>=.5*len(gs)
    psummary.append({**base0,'perturbed_gene':k,'endpoint_gene':program+'_'+mask,'program':program,'mask':mask,'metric_type':'program_median_heat_distance','effect_scale':'median regulatory distance across covered non-KO members','effect_value':np.median(values) if evaluable else None,'assessment_status':'EVALUABLE' if evaluable else 'NE','ne_reason':None if evaluable else ('no members after MP10 mask' if not gs else 'state/coverage/path gate insufficient'),'covered_members':','.join(ok),'n_covered':len(ok),'n_frozen':len(gs)})
 for g in sorted(set(targets).union(*progs.values())):
  label=base0['state'] if j!='main' else 'ALL'
  x=rawcov[(rawcov.state==label)&(rawcov.gene==g)].iloc[0].to_dict()
  for k in K:
   y={'source_outdegree':cov.loc[k,'outdegree'] if not cov.empty and k in cov.index else None,'endpoint_indegree':cov.loc[g,'indegree'] if not cov.empty and g in cov.index else None}
   present=not cov.empty and g in cov.index and bool(cov.loc[g,'in_network'])
   reason=gene_result(base0,k,g,eff,cov).get('ne_reason')
   coverage.append({**x,'analysis_id':j,'perturbed_gene':k,'raw_state':label,'prior_id':'none_de_novo','prior_used':False,'external_collectri_target':bool(((frozen.source==k)&(frozen.target==g)).any()),'prior_coverage_note':'No prior required by scTenifoldKnk; external CollecTRI listed for provenance only','in_hvg_top2000':bool(panel.loc[g,'hvg_top2000']) if g in panel.index else False,'reserved_in_panel':bool(g in panel.index and panel.loc[g,'reserved']),'in_model_panel':g in panel.index,'in_fitted_network':present,'fitted_n_positive':cov.loc[g,'n_positive'] if not cov.empty and g in cov.index else None,'fitted_positive_patients':cov.loc[g,'n_positive_patients'] if not cov.empty and g in cov.index else None,**y,'path_available':None if not present or eff.empty or (k,g) not in eff.index else bool(eff.loc[(k,g),'reachable']),'assessment_status':'NE' if reason else 'EVALUABLE','ne_reason':reason})
R=pd.DataFrame(robust);R['analysis_type']=R.analysis_id.map(lambda k:'patient_subsample' if k.startswith('sub') else 'LODO' if k.startswith('drop_') else 'alternative_seed')
S=pd.DataFrame(rows);main=S[(S.analysis_id=='main')&(S.endpoint_gene=='NECTIN2')].set_index('perturbed_gene')
RR=[]
for (k,mode),a in R.groupby(['perturbed_gene','analysis_type']):
 ok=a[a.assessment_status=='EVALUABLE'];d=ok.dropna(subset=['predicted_direction'])
 ref=main.loc[k,'predicted_direction']
 RR.append(dict(perturbed_gene=k,analysis_type=mode,n_planned=len(a),n_evaluable=len(ok),n_direction=len(d),n_up=sum(d.predicted_direction=='up'),n_down=sum(d.predicted_direction=='down'),sign_agreement=(d.predicted_direction==ref).mean() if len(d) else None,distance_min=ok.effect_value.min(),distance_max=ok.effect_value.max(),direction_score_min=ok.direction_score.min(),direction_score_max=ok.direction_score.max(),rank_min=ok['rank'].min() if len(ok) else None,rank_max=ok['rank'].max() if len(ok) else None,n_q05=(ok.q_value<.05).sum(),depth_sign_agreement=(np.sign(d.depth_score)==np.sign(d.direction_score)).mean() if len(d) else None,pc1_sign_agreement=(np.sign(d.pc1_score)==np.sign(d.direction_score)).mean() if len(d) else None))
 for ix in a.index:R.loc[ix,'direction_matches_main']=bool(R.loc[ix,'predicted_direction']==ref) if R.loc[ix,'direction_available'] else None
 if mode=='patient_subsample' and len(ok):
  i=(S.analysis_id=='main')&(S.endpoint_gene=='NECTIN2')&(S.perturbed_gene==k)&(S.assessment_status=='EVALUABLE')
  S.loc[i,['uncertainty_type','uncertainty_low','uncertainty_high']]=['patient_subsample_refit_min_max',ok.effect_value.min(),ok.effect_value.max()]
save(R,'tables/virtual_perturbation_robustness.tsv');save(pd.DataFrame(RR),'tables/refit_ranges.tsv')
PP=pd.DataFrame(psummary);save(PP,'tables/program_response.tsv')
C=pd.DataFrame(coverage);save(C,'tables/virtual_perturbation_coverage.tsv')
# Expand every patient-by-stratum combination so that sparse strata remain explicit.
pr=job('main','patient');designs=[]
for p,grp in mainmeta.groupby('patient_id'):
 for st in ['All','H05','H10','H50','L50']:
  cc=len(grp) if st=='All' else int(grp[st].sum())
  for k in K:designs.append(dict(patient_id=p,dataset=grp.dataset.iloc[0],perturbed_gene=k,stratum=st,n_cells=cc))
pr=pd.DataFrame(designs).merge(pr,on=['patient_id','dataset','perturbed_gene','stratum','n_cells'],how='left',validate='one_to_one')
pr.loc[pr.n_cells<10,'reason']='fewer than 10 cells';pr['assessment_status']=np.where(pr.effect_value.isna(),'NE','EXPLORATORY')
save(pr,'tables/patient_conditioned.tsv')
pat=[]
for (k,st),a in pr.groupby(['perturbed_gene','stratum']):
 v=a.dropna(subset=['effect_value']);base0=design('main')
 pat.append({**base0,'analysis_id':'patient_median_'+st,'state':st,'perturbed_gene':k,'endpoint_gene':'NECTIN2','metric_type':'median_patient_conditioned_distance','effect_scale':'median across patients on fixed WT heat kernel','effect_value':v.effect_value.median() if len(v) else None,'n_patients':len(v),'n_cells':v.n_cells.sum(),'assessment_status':'EXPLORATORY' if len(v) else 'NE','ne_reason':None if len(v) else 'no eligible patients','uncertainty_type':'observed_patient_min_max' if len(v) else None,'uncertainty_low':v.effect_value.min(),'uncertainty_high':v.effect_value.max()})
S=pd.concat([S,PP,pd.DataFrame(pat)],ignore_index=True);save(S,'tables/virtual_perturbation_summary.tsv')
# Direct plotting data, all axes and counts kept in native units.
coldata=C[(C.analysis_id=='main')&(C.perturbed_gene=='SREBF1')&C.gene.isin(targets)].copy()
save(coldata,'plot_data/F01_coverage.tsv')
save(pr[pr.stratum=='All'],'plot_data/F02_patients.tsv')
save(S[(S.analysis_id.str.startswith('state_'))&(S.endpoint_gene=='NECTIN2')],'plot_data/F03_states.tsv')
co=controls.rename(columns={'gene':'perturbed_gene'}).merge(main.reset_index(),on='perturbed_gene',validate='many_to_one')
mt=main.loc[K].reset_index();mt['matched_to']=mt.perturbed_gene;mt['role']='SREBF';co['role']='Control'
co=pd.concat([co,mt],ignore_index=True);save(co,'plot_data/F04_controls.tsv');save(controls,'tables/matched_controls.tsv')
save(R,'plot_data/F05_refits.tsv');save(S[(S.analysis_id=='main')&S.perturbed_gene.isin(K)&S.endpoint_gene.isin(endpoints)],'plot_data/F06_secondary.tsv')
save(PP[PP.analysis_id=='main'],'plot_data/F07_programs.tsv')
save(pd.DataFrame(pat).query("state!='All'"),'plot_data/F08_strata.tsv')
save(main.loc[K].reset_index().melt(id_vars=['perturbed_gene'],value_vars=['direction_score','depth_score','pc1_score'],var_name='adjustment',value_name='score'),'plot_data/F09_direction.tsv')
jobs=['pilot',*mainjobs,*other]
resources=[];stages=[]
for j in jobs:
 s=Path('logs/usage_'+j+'.txt').read_text()
 key=lambda field:re.search(field+r':\s*([^\n]+)',s)[1].strip()
 tm=key(r'Elapsed \(wall clock\) time \(h:mm:ss or m:ss\)');w=sum(float(x)*60**i for i,x in enumerate(reversed(tm.split(':'))))
 resources.append(dict(analysis_id=j,wall_s=w,cpu_s=float(key(r'User time \(seconds\)'))+float(key(r'System time \(seconds\)')),peak_RSS_KiB=int(key(r'Maximum resident set size \(kbytes\)')),output_bytes=sum(p.stat().st_size for p in J.glob(j+'_*'))+(Path('cache/jobs/main_WT.rds').stat().st_size if j=='main' else 0),**job(j,'status').iloc[0].to_dict()))
 tf=Path('logs/time_'+j+'.tsv')
 if tf.exists():stages.append(rd(tf).assign(analysis_id=j))
save(pd.DataFrame(resources),'tables/resource_usage.tsv');save(pd.concat([rd('logs/stages.tsv').assign(analysis_id='extraction'),*stages]),'tables/stage_resources.tsv')
plan=rd('tables/resource_plan.tsv');t=rd('logs/time_pilot.tsv').wall_s.sum();plan['pilot_stage_wall_s']=t;plan['estimated_wall_s']=t*(2000/plan.pilot_genes)**3*plan.formal_jobs/3;save(plan,'tables/resource_plan.tsv')

