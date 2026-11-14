import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,hashlib
import pandas as pd,numpy as np
from scipy.stats import t
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v1';O=P/'results/mp10_immune_endothelial_mechanisms_v1';checks=[]
def rd(n):return pd.read_csv(B/'tables'/n,sep='\t')
def check(name,value,detail=''):
 checks.append(dict(check=name,passed=bool(value),detail=str(detail)))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
r=rd('model_results.tsv.gz');reg=rd('model_registry.tsv');m=rd('model_membership.tsv.gz');cy=rd('treg_cycle_composition.tsv');labels=rd('treg_cell_cycle_labels.tsv.gz');s=rd('sample_level_scores.tsv.gz');q=pd.read_csv(O/'tables/specimen_exposure_eligibility.tsv',sep='\t')
check('registry_complete_unique',set(r.model_id)==set(reg.model_id) and r.model_id.is_unique and reg.model_id.is_unique,len(r))
for fam,n in [('epi_continuous_primary',10),('bridge_primary',18),('bridge_noncycling',9),('treg_secondary',45)]:
 d=r[r.FDR_family==fam];check('planned_'+fam,len(d)==n and d.planned_n.eq(n).all(),f'{len(d)} registered / {sum(d.P.notna())} testable')
for fam,d in r.groupby('FDR_family'):
 p=d.P.to_numpy();idx=np.flatnonzero(np.isfinite(p));qq=np.full(len(d),np.nan);order=idx[np.argsort(p[idx])]
 if len(order):qq[order]=np.minimum(1,np.minimum.accumulate((p[order]*d.planned_n.iloc[0]/np.arange(1,len(order)+1))[::-1])[::-1])
 check('BH_'+fam,np.allclose(qq,d.q,equal_nan=True,atol=1e-12) and d.planned_n.nunique()==1 and len(d)==d.planned_n.iloc[0])
check('NE_not_fake_P1',r.loc[r.status=='not_evaluable',['beta','P','q','CI_low','CI_high']].isna().all().all())
check('completed_df_and_units',((r.loc[r.status=='completed','residual_df']>=3)&(r.loc[r.status=='completed','n_recorded_patient_units']>=4)).all())
check('no_repeated_patient_in_model',not m.duplicated(['model_id','patient_uid_recorded']).any())
check('unique_specimen_in_model',not m.duplicated(['model_id','sample_uid']).any())
check('a_b_not_both',all(not (d.sample_uid.str.contains('T4463209a').any() and d.sample_uid.str.contains('T4463209b').any()) for _,d in m.groupby('model_id')))
check('phase_all_counted',(cy.n_cycling+cy.n_noncycling+cy.n_invalid).eq(cy.n_Treg).all() and cy.n_Treg.sum()==len(labels))
check('phase_barcodes_unique',labels.barcode.is_unique)
check('no_invalid_silently_noncycling',labels.loc[~labels.valid_phase,'stratum'].eq('invalid').all())
check('fraction_valid_zero_semantics',np.allclose(cy.loc[cy.n_invalid.eq(0),'cycling_fraction'],cy.loc[cy.n_invalid.eq(0),'n_cycling']/cy.loc[cy.n_invalid.eq(0),'n_Treg']))
check('phase_G1_rule',(labels.loc[labels.Phase.eq('G1'),'S.Score']<=1e-12).all() and (labels.loc[labels.Phase.eq('G1'),'G2M.Score']<=1e-12).all())
a=pd.read_csv(O/'tables/analysis_cell_metadata.tsv.gz',sep='\t',usecols=['barcode','analysis_celltype_c2l_final','sample_uid']);check('exact_original_Treg_membership',set(a.loc[a.analysis_celltype_c2l_final.eq('Treg'),'barcode'])==set(labels.barcode))
check('classification_complete',cy.n_invalid.eq(0).all())
for (fam,metric),d in s.groupby(['family','metric']):
 ref=d[d.profile_type.eq(fam)&d.main_patient.eq(True)&d.n_cells.ge(50)]
 check('fixed_score_reference_'+fam+'_'+metric,d.center.nunique()==1 and d.scale.nunique()==1 and abs(ref.score.mean())<1e-10 and abs(ref.score.std()-1)<1e-10)
check('no_duplicate_scores',not s.duplicated(['profile_id','metric']).any())
for _,row in r[r.FDR_family.str.startswith('bridge_noncycling')].iterrows():
 ctrl=r[r.FDR_family.str.startswith('bridge_matched_all')&r.x.eq(row.x)&r.y.eq(row.y)&r.setting.eq(row.setting)]
 check('matched_control_registered_'+row.model_id,len(ctrl)==1)
 if len(ctrl)!=1:continue
 da=m[m.model_id.eq(row.model_id)].sort_values('sample_uid');db=m[m.model_id.eq(ctrl.model_id.iloc[0])].sort_values('sample_uid')
 check('matched_members_'+row.model_id,da.sample_uid.tolist()==db.sample_uid.tolist())
 if len(da)!=len(db):continue
 check('matched_exposure_'+row.model_id,np.allclose(da.x,db.x,equal_nan=True))
 check('matched_dataset_'+row.model_id,da.dataset.tolist()==db.dataset.tolist())
 if row.covariate!='none':check('matched_covariate_'+row.model_id,np.allclose(da[row.covariate],db[row.covariate],equal_nan=True))
# Independent matrix implementation of all estimable HC3 fits.
for _,row in r[r.status.eq('completed')&r.method.eq('HC3')].iterrows():
 d=m[m.model_id.eq(row.model_id)].sort_values('sample_uid');n=len(d);X=np.column_stack([np.ones(n),pd.get_dummies(d.dataset,drop_first=True,dtype=float).to_numpy(), *([d[row.covariate].to_numpy()] if row.covariate!='none' else []),d.x.to_numpy()]);y=d.y.to_numpy();inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;resid=y-X@beta;h=np.sum((X@inv)*X,axis=1);meat=sum(np.outer(X[i],X[i])*(resid[i]/(1-h[i]))**2 for i in range(n));se=np.sqrt((inv@meat@inv)[-1,-1]);df=n-X.shape[1];pv=2*t.sf(abs(beta[-1]/se),df)
 check('HC3_independent_'+row.model_id,np.allclose([beta[-1],se,pv],[row.beta,row.se,row.P],rtol=1e-6,atol=1e-8))
# Independent IRLS check of the three principal proportion coefficients.
for id in ['F00012','F00039','F00066']:
 row=r[r.model_id.eq(id)].iloc[0];d=m[m.model_id.eq(id)].merge(cy[['sample_uid','n_cycling','n_noncycling']],on='sample_uid').sort_values('sample_uid');n=len(d);X=np.column_stack([np.ones(n),pd.get_dummies(d.dataset,drop_first=True,dtype=float).to_numpy(),d.x]);N=(d.n_cycling+d.n_noncycling).to_numpy();yy=d.n_cycling.to_numpy();beta=np.zeros(X.shape[1])
 for _ in range(100):
  pp=1/(1+np.exp(-(X@beta)));ww=N*pp*(1-pp);delta=np.linalg.solve(X.T@(ww[:,None]*X),X.T@(yy-N*pp));beta+=delta
  if abs(delta).max()<1e-12:break
 pp=1/(1+np.exp(-(X@beta)));ww=N*pp*(1-pp);df=n-X.shape[1];disp=np.sum((yy-N*pp)**2/ww)/df;se=np.sqrt(np.linalg.inv(X.T@(ww[:,None]*X))[-1,-1]*disp)
 check('quasi_independent_'+id,np.allclose([beta[-1],se,disp],[row.beta,row.se,row.dispersion],rtol=1e-5,atol=1e-8))
d=r[r.status.eq('completed')];crit=t.ppf(.975,d.residual_df);check('CI_and_P_same_t_df',np.allclose(d.CI_low,d.beta-crit*d.se) and np.allclose(d.CI_high,d.beta+crit*d.se) and np.allclose(d.P,2*t.sf(abs(d.beta/d.se),d.residual_df)))
repair=rd('gsea_numeric_repair_comparison.tsv.gz');check('repaired_all_original_NA',repair.loc[repair.pval_old.isna(),'pval_repaired'].notna().all());check('repair_mask_scope_complete',set(repair['mask'])=={'Full','M20-masked'})
check('old_input_files_unchanged',all(sha(pathlib.Path(x.path))==x.sha256 for x in rd('input_reuse_audit.tsv').itertuples()))
# The original large object is verified, not copied. One SHA pass per final audit.
orig=rd('original_source_input_manifest.tsv');obj=orig[orig.path.str.endswith('final_annotation_v2_pre_reference.qs')].iloc[0];check('original_raw_object_sha256',sha(pathlib.Path(obj.path))==obj.sha256)
check('M0_original_registry_preserved',sha(B/'tables/model_registry_v1_frozen.tsv')==json.loads((B/'config/analysis_plan.yaml').read_text())['registry_sha256'])
for n in ['README.md','METHODS.md','RESULTS.md','EVIDENCE_REVIEW.md','HYPOTHESIS_v1.md','HYPOTHESIS_UPDATE.md','HD_VALIDATION_PENDING.md','FIGURE_INDEX.md']:
 check('required_document_'+n,(B/n).exists())
for fam,tf in [('SCENIC','SREBF2')]:check('SREBF2_SCENIC_NE_retained','SREBF2' in (B/'RESULTS.md').read_text() and 'NE' in (B/'RESULTS.md').read_text())
ov=rd('gene_overlap.tsv');check('SCENIC_SREBF2_not_zero_filled',ov.loc[ov.set1.eq('SCENIC_SREBF2')|ov.set2.eq('SCENIC_SREBF2'),'n_shared'].isna().all())
check('alias_ACAT2_not_SOAT2',set(rd('gene_alias_audit.tsv').gene)=={'NECTIN2','PVRL2','ACAT2','SOAT2'})
th=pd.read_csv(B/'logs/thread_audit.tsv',sep='\t');check('all_runtime_threads_one',th[['BLAS','OMP','data_table']].eq(1).all().all())
put=pd.DataFrame(checks);put.to_csv(B/'qc/validation.tsv',sep='\t',index=False);print('Validation',sum(put.passed),'/',len(put));print(put[~put.passed].to_string(index=False));assert put.passed.all()
