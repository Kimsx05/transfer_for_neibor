import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,hashlib,pandas as pd,numpy as np
from scipy import stats
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';T=B/'tables';checks=[]
def ck(name,ok,detail=''):checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
 return h.hexdigest()
def bh(p,n):
 p=np.asarray(p,float);ids=np.flatnonzero(np.isfinite(p));ids=ids[np.argsort(p[ids])];q=np.full(len(p),np.nan)
 if len(ids):q[ids]=np.minimum(1,np.minimum.accumulate((p[ids]*n/np.arange(1,len(ids)+1))[::-1])[::-1])
 return q
reg=pd.read_csv(T/'model_registry.tsv',sep='\t');score=pd.read_csv(T/'score_model_results.tsv.gz',sep='\t');gs=pd.read_csv(T/'gsea_model_results.tsv.gz',sep='\t');scenic=pd.read_csv(T/'recovered_SCENIC_model_results.tsv',sep='\t');res=pd.concat([score,gs,scenic],ignore_index=True);mem=pd.read_csv(T/'model_membership.tsv.gz',sep='\t');scmem=pd.read_csv(T/'recovered_SCENIC_model_membership.tsv.gz',sep='\t')
ck('registry_exact_result_ids',set(res.model_id)==set(reg.model_id));ck('unique_model_ids',res.model_id.is_unique);ck('NE_P_q_missing',res.loc[res.status.ne('completed'),['P','q']].isna().all().all());ck('H10_H50_not_v2_primary',not reg.loc[reg.state.isin(['H10','H50']),'primary'].any())
for fam,d in res.groupby('FDR_family'):
 n=reg.loc[reg.FDR_family.eq(fam),'planned_n'].iloc[0];ck('BH_'+fam,np.allclose(bh(d.P,n),d.q,equal_nan=True,rtol=1e-9,atol=1e-12));ck('planned_'+fam,len(d)==n)
for mid,d in mem.groupby('model_id'):
 r=score[score.model_id.eq(mid)].iloc[0];ck('patient_unique_'+mid,d.patient_uid_recorded.is_unique);ck('cells_'+mid,(d.n_sender>=r.min_cells).all() and(d.n_receiver>=r.min_cells).all());ck('n_'+mid,len(d)==r.n_specimens)
 if r.status!='completed':continue
 paired=r.kind!='cross';cols=[np.ones(len(d))]
 if not paired and d.dataset.nunique()>1:cols.extend(pd.get_dummies(d.dataset,dtype=float).iloc[:,1:].to_numpy().T)
 if r.setting in ['quality_RNA','composition_cluster0','quality_sender_RNA','quality_receiver_RNA']:cols.append(d.covariate.to_numpy())
 if not paired:cols.append(d.x.to_numpy())
 X=np.column_stack(cols);y=d.y.to_numpy();n,k=X.shape;df=n-np.linalg.matrix_rank(X);ix=0 if paired else k-1
 if r.y=='cycling_fraction':
  nt=(d.n_cycling+d.n_noncycling).to_numpy();beta=np.zeros(k)
  for iteration in range(300):
   eta=np.clip(X@beta,-30,30);mu=1/(1+np.exp(-eta));w=nt*mu*(1-mu);I=X.T@(X*w[:,None]);grad=X.T@(nt*(y-mu));step=np.linalg.solve(I,grad);beta+=step
   if np.max(abs(step))<1e-10:break
  mu=1/(1+np.exp(-np.clip(X@beta,-30,30)));w=nt*mu*(1-mu);disp=np.sum(nt*(y-mu)**2/(mu*(1-mu)))/df;V=np.linalg.inv(X.T@(X*w[:,None]))*disp;b=beta[ix];se=np.sqrt(V[ix,ix]);tol=2e-5
 else:
  V0=np.linalg.inv(X.T@X);coef=V0@X.T@y;err=y-X@coef;h=np.einsum('ij,jk,ik->i',X,V0,X);V=V0@(X.T@(X*(err**2/(1-h)**2)[:,None]))@V0;b=coef[ix];se=np.sqrt(V[ix,ix]);tol=1e-8
 p=2*stats.t.sf(abs(b/se),df);ci=b+np.array([-1,1])*stats.t.ppf(.975,df)*se
 ck('numeric_'+mid,np.allclose([b,se,*ci,p],[r.beta,r.SE,r.CI_low,r.CI_high,r.P],rtol=tol,atol=tol),f'Python independent {r.method}')
 if r.kind=='direct_D':ck('direct_D_identity_'+mid,np.allclose(d.y,d.delta_sterol-d.delta_FA))
for mid,d in scmem.groupby('model_id'):
 r=scenic[scenic.model_id.eq(mid)].iloc[0];ck('SCENIC_patient_'+mid,d.patient_uid_recorded.is_unique)
 if r.status=='completed':
  n=len(d);b=d.delta.mean();se=np.sqrt(np.sum((d.delta-b)**2)/(1-1/n)**2/n**2);ck('SCENIC_numeric_'+mid,np.allclose([b,se,2*stats.t.sf(abs(b/se),n-1)],[r.beta,r.SE,r.P],atol=1e-10))
a=pd.read_csv(T/'scenic_module_ledger.tsv',sep='\t');b=pd.read_csv(T/'scenic_corrected_module_ledger.tsv',sep='\t');ck('original_module_unique',a.module_id.is_unique);ck('same_module_denominator',a.set_index('module_id').n_original.equals(b.set_index('module_id').n_original));ck('all_original_SREBF2_mapping_fail',(~a.loc[a.TF.eq('SREBF2'),'mapping_pass']).all());ck('all_original_SREBF2_NES_annotation_present',(a.loc[a.TF.eq('SREBF2'),'n_enriched_annotated']>0).all());ck('SREBF2_top50_rescue_41of51',len(b[b.TF.eq('SREBF2')&b.context.eq('activating;top50')&b.n_mapped_corrected.eq(41)&b.n_original.eq(51)&b.status.eq('retained')])==1)
for mask in ['Full','M20']:
 d=pd.read_csv(T/f'AUCell_original_numerical_check_{mask}.tsv',sep='\t');ck('original_AUCell_numeric_'+mask,(d.max_absolute_error<1e-10).all() and d.n_checked.min()==3000)
coverage=pd.read_csv(T/'signature_coverage.tsv',sep='\t');ck('fixed_signature_coverage_rule',np.array_equal(coverage.evaluable,(coverage.n_covered>=5)&(coverage.n_covered/coverage.n_expected>=.5)));ck('two_primary_routes_NE',len(res[res.FDR_family.eq('cross_primary_main')&res.status.eq('not_evaluable')])==2)
f=pd.read_csv(T/'mTORC1_EC_filter_audit.tsv',sep='\t');ck('EC_filter176_98',f.old_design_keep.sum()==176 and f.new_fixed_keep.sum()==98);ck('EC_counterfactual102_157',f.old_samples_group_keep.sum()==102 and f.new_samples_design_keep.sum()==157)
for r in pd.read_csv(T/'input_manifest.tsv',sep='\t').itertuples():ck('input_unchanged_'+pathlib.Path(r.path).name,sha(r.path)==r.sha256)
frozen=json.loads((B/'config/M0_FROZEN.json').read_text())
for p,h in frozen.items():
 actual=B/p
 if p=='tables/model_registry.tsv':actual=T/'model_registry_M0_original.tsv'
 ck('M0_original_'+p,sha(actual)==h)
thread=pd.read_csv(B/'logs/thread_audit.tsv',sep='\t');ck('R_threads_one',(thread[['BLAS','OMP','data_table']]==1).all().all())
x=pd.DataFrame(checks);x.to_csv(B/'qc/validation.tsv',sep='\t',index=False);print('VALIDATION',x.passed.sum(),'/',len(x));print(x[~x.passed].to_string(index=False));assert x.passed.all()
