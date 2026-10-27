from common import *
checks=[]
def check(name,ok,details=''):
 checks.append(dict(check=name,passed=bool(ok),details=details));assert ok,(name,details)
s=read('samples_v2.tsv');m=read('cells_v2.tsv.gz');p=read('proportions_v2.tsv.gz');v=read('score_validation.tsv');r=read('proportion_models_v2.tsv');memb=read('model_sample_membership.tsv')
check('unique barcode',m.barcode.is_unique);check('counts overlap identical',v.count_max_abs_error.eq(0).all());check('score overlap identical',v.score_abs_error.max()<1e-12);check('all epi finite coverage',m.loc[m.major.eq('Epithelial'),'MP10_cell'].notna().all());check('TME state never assigned',m.loc[m.major.ne('Epithelial'),'epi_MP10_state_within_specimen'].isna().all());check('contamination exclusion retained',m.annotation_status_v2.ne('excluded_before_round2').all());check('reference flag not filter',m.spatial_reference_include.eq(False).any());check('proportions equal v1 raw counts',p.equals(pd.read_csv(V/'tables/proportions.tsv.gz',sep='\t')))
for keys,a in r.groupby(['analysis_variant','comparison','scope','annotation_level','denominator_type']):
 pp=a.P.dropna();ix=np.argsort(pp.to_numpy());out=np.empty(len(pp));out[ix]=np.minimum.accumulate((pp.to_numpy()[ix]*len(pp)/np.arange(1,len(pp)+1))[::-1])[::-1].clip(0,1)
 check('BH family '+str(keys),np.allclose(a.loc[pp.index,'q'],out))
base=['scope','annotation_level','denominator_type','parent_family']
for keys,a in memb.groupby(base):
 for var in a.analysis_variant.unique():
  b=a[a.analysis_variant.eq(var)];x=set(b.loc[b.comparison.eq('MP10_stage_complete_unadjusted'),'sample_uid']);y=set(b.loc[b.comparison.eq('MP10_stage_adjusted'),'sample_uid']);check('stage matched '+str(keys)+' '+var,x==y)
 for co in a.comparison.unique():
  b=a[a.comparison.eq(co)];x=set(b.loc[b.analysis_variant.eq('parent30_asin'),'sample_uid']);y=set(b.loc[b.analysis_variant.eq('parent30_logit'),'sample_uid']);check('transform matched '+str(keys)+' '+co,x==y)
if (T/'network_status.tsv').exists():
 n=read('network_status.tsv');plan=read('network_plan.tsv');check('all planned networks complete',len(n)==len(plan) and n.status.isin(['completed','reused']).all());con=read('network_cell_contributions.tsv.gz');check('sampling cap100',con.n_sampled.le(100).all());state=con[con.network_id.str.contains('epi_state_within_specimen')&con.node.isin(['Epi_MP10_low','Epi_MP10_high'])];ct=state.pivot(index=['network_id','sample_uid'],columns='node',values='n_sampled');check('equal epi sampling per specimen',ct.Epi_MP10_low.eq(ct.Epi_MP10_high).all());check('sampling original counts',con.n_sampled.le(con.n_raw).all())
 c=read('B_Plasma_LR_contrasts.tsv.gz');check('raw difference never threshold masked',np.allclose(c.raw_delta,c.prob_high-c.prob_low,equal_nan=True));f=read('B_Plasma_LR_all_statuses.tsv.gz');check('uncomputed probabilities NA',f.loc[f.status.isin(['node_insufficient','gene_unmeasured_or_unmatched','gene_measurement_unresolved','LR_screened_out']),'raw_probability'].isna().all())
 circle=read('circle_nodes_v2.tsv');a=circle[circle.comparison_type.eq('epi_state_within_specimen')];check('both epi nodes in all state panels',a.groupby(['axis_id','scope','condition']).node.apply(lambda x:{'Epi_MP10_low','Epi_MP10_high'}.issubset(x)).all())
save(checks,'validation_checks.tsv');print('PASS',len(checks),'checks')
