exec(open(__file__.replace('07_validate.py','01_prepare.py')).read().split('# Archive contents')[0])
import platform,subprocess,datetime
from PIL import Image
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,pass_=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError((name,detail))
s=pd.read_csv(B/'02_cell_lines/sample_scores.tsv',sep='\t');check('60_GSM_unique',len(s)==60 and s.GSM.is_unique and s.matrix_id.is_unique);check('5x4x3_design',s.groupby(['cell_line','state']).size().eq(3).all() and s.cell_line.nunique()==5 and s.state.nunique()==4);check('RNA_scores_finite',np.isfinite(s.MP10).all())
d=pd.read_csv(B/'02_cell_lines/line_deltas.tsv',sep='\t');check('all_15_paired_line_deltas',len(d)==15 and np.allclose(d.delta,d.resistant_mean-d.parental));su=pd.read_csv(B/'02_cell_lines/overall_effects.tsv',sep='\t');check('sign_test_n5_lower_bound',su.P.min()>=.0625 and su.n_lines.eq(5).all())
a=pd.read_csv(B/'03_clinical/patient_inclusion_endpoint_audit.tsv',sep='\t');check('188_GSM_unique',len(a)==188 and a.GSM.is_unique);check('within_cohort_patient_unique',not a.duplicated(['cohort','patient_id']).any());check('endpoint_complements',np.allclose(a.nonresponse,1-a.good,equal_nan=True) and np.allclose(a.non_pCR,1-a.pCR,equal_nan=True));check('pCR_implies_good',a.loc[a.pCR.eq(1),'good'].eq(1).all());h=a[a.cohort=='GSE212810'];check('39_AMC_35_cisplatin_4_GCb',len(h)==39 and h.strict_inclusion.sum()==35 and h.regimen.eq('GCb_carboplatin').sum()==4);check('GCb_excluded_strict',not h.loc[h.regimen=='GCb_carboplatin','strict_inclusion'].any());check('author_21R_18NR',(h.author_good==1).sum()==21 and (h.author_good==0).sum()==18);check('known_NR_not_lost_to_N_missing',h.loc[h.ypT_code.ge(2)&h.ypT_code.ne(777),'good'].eq(0).all());check('ypT0Nx_unknown',h.loc[h.ypT_code.eq(777)&h.ypN_code.eq(999),['good','pCR']].isna().all().all());check('noninvasive_Nx_not_pCR',h.loc[h.ypT_code.eq(0)&h.ypN_code.eq(999),'pCR'].eq(0).all());g=a[a.cohort=='GSE169455'];check('125_NAC_24_induction',g.therapy_intent.eq('neoadjuvant').sum()==125 and g.therapy_intent.eq('induction').sum()==24)
sc=pd.read_csv(B/'03_clinical/score_standardization.tsv',sep='\t')
for co,z in a.groupby('cohort'):
 ss=sc[sc.cohort==co].iloc[0];check('fixed_z_'+co,np.allclose(z.MP10_z,(z.MP10-ss['mean'])/ss.SD))
r=pd.read_csv(B/'03_clinical/model_results.tsv',sep='\t');use=pd.read_csv(B/'03_clinical/model_patient_usage.tsv',sep='\t');check('model_uses_each_patient_once',not use.duplicated(['cohort','scope','endpoint','adjustment','patient_id']).any());check('finite_estimates_for_all_completed_models',np.isfinite(r.loc[r.converged.fillna(False),['OR','CI_low','CI_high','P']]).all().all());check('CIs_contain_point',((r.CI_low<=r.OR)&(r.OR<=r.CI_high)|r.OR.isna()).all());check('prespecified_Firth_applied',r.method.eq('Firth_Jeffreys_profile_likelihood').sum()==1)
for _,row in r[r.converged.fillna(False)].iterrows():
 z=use[(use.cohort==row.cohort)&(use.scope==row.scope)&(use.endpoint==row.endpoint)&(use.adjustment==row.adjustment)];check('model_count_'+str(row.name),len(z)==row.n and z.response.sum()==row.events)
e=pd.read_csv(B/'01_epi/sample_effects.tsv',sep='\t');q=json.loads((B/'01_epi/execution_QC.json').read_text());check('MP10_exact_replay',q['MP10_reproduction_max_abs_error']<1e-12);check('no_lowN_effects',e.loc[e.n_used<50,'rho'].isna().all());check('missing_cycle_not_imputed',e.loc[e.adjustment!='unadjusted','rho'].isna().all());p=pd.read_csv(B/'01_epi/patient_summary.tsv',sep='\t');check('no_false_patient_inference',p.n_patients.eq(0).all() and p[['rho','P','CI_low','CI_high']].isna().all().all());check('no_fake_cell_P','P' not in e.columns)
pre=json.loads((B/'00_audit/pre_replay_hashes.json').read_text());rr=[]
for f,old in pre.items():
 sha=hashlib.sha256((B/f).read_bytes()).hexdigest();rr.append(dict(path=f,before=old,after=sha,identical=old==sha))
save(pd.DataFrame(rr),'00_audit/reproducibility.tsv');check('core_tables_identical_on_full_replay',all(x['identical'] for x in rr))
figs=[]
for f in sorted((B/'figures').glob('*.png')):
 im=Image.open(f);dpi=im.info.get('dpi',(0,0));check('PNG_300dpi_'+f.stem,abs(dpi[0]-300)<1);check('PNG_white_'+f.stem,im.convert('RGB').getpixel((0,0))==(255,255,255));check('paired_PDF_'+f.stem,f.with_suffix('.pdf').exists());figs.append(dict(file=f.name,width=im.width,height=im.height,dpi=dpi[0]))
save(pd.DataFrame(figs),'00_audit/figure_validation.tsv')
resource=[]
for ld in sorted((B/'logs').glob('replay_*')):
 if not ld.is_dir():continue
 for f in ld.glob('*.log'):
  text=f.read_text();m=re.search(r'Maximum resident set size \(kbytes\): (\d+)',text);exit_=re.search(r'Exit status: (\d+)',text)
  if m:resource.append(dict(step=f.stem,peak_RSS_KiB=int(m[1]),exit_status=int(exit_[1]) if exit_ else None,log=str(f.relative_to(B))))
save(pd.DataFrame(resource),'00_audit/resource_usage.tsv');check('resource_cap',max(z['peak_RSS_KiB'] for z in resource)<12*1024**2);check('all_replay_steps_zero_exit',all(z['exit_status']==0 for z in resource))
save(pd.DataFrame(checks).rename(columns={'pass_':'pass'}),'00_audit/validation_checks.tsv');(B/'00_audit/python_environment.json').write_text(json.dumps(dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,platform=platform.platform()),indent=2));print('VALIDATED',len(checks),'checks; all core table hashes identical')
