import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,pandas as pd,numpy as np,json
B=pathlib.Path.cwd()/'results/mp10_hypothesis_followup_v2';T=B/'tables';checks=[]
def ck(n,v):checks.append(dict(check=n,passed=bool(v)))
r=pd.read_csv(T/'model_results.tsv.gz',sep='\t');reg=pd.read_csv(T/'model_registry_complete.tsv',sep='\t');ck('1328_registered_and_reported',len(r)==1328 and set(r.model_id)==set(reg.model_id));ck('NE_not_zero',r.loc[r.status!='completed',['P','q']].isna().all().all());ck('1073_testable_255_NE',sum(r.status=='completed')==1073 and sum(r.status=='not_evaluable')==255)
for fam,d in r.groupby('FDR_family'):
 n=int(d.planned_n.iloc[0]);p=d.P.dropna().sort_values();q=np.minimum(1,np.minimum.accumulate((p.values*n/np.arange(1,len(p)+1))[::-1])[::-1]);ck('BH_'+fam,np.allclose(q,d.loc[p.index,'q'],atol=1e-10));ck('family_registered_n_'+fam,len(d)==n)
for fn in ['AUCell_original_numerical_check_Full.tsv','AUCell_original_numerical_check_M20.tsv']:
 d=pd.read_csv(T/fn,sep='\t');ck(fn,(d.max_absolute_error==0).all())
for mask in ['Full','M20']:
 d=pd.read_csv(T/f'SREBF_corrected_cell_AUCell_{mask}.tsv.gz',sep='\t');ck('AUC116448_'+mask,len(d)==116448 and d.barcode.is_unique)
fig=pd.read_csv(T/'figure_index.tsv',sep='\t');ck('18_figures',len(fig)==18)
for f in fig.itertuples():
 d=pd.read_csv(B/'plot_data'/f.plot_data,sep='\t');ck('figure_has_data_'+f.figure,len(d)>0)
 if 'model_id' in d:
  src=r.set_index('model_id');d=d[d.model_id.isin(src.index)]
  for k in ['beta','CI_low','CI_high','q','P']:
   if k in d:ck('figure_matches_'+f.figure+'_'+k,np.allclose(d[k],src.loc[d.model_id,k],equal_nan=True,atol=1e-10))
rd=pd.read_csv(B/'qc/render_audit.tsv',sep='\t');ck('all62_vector_png_checks',len(rd)==62 and rd.passed.all());ck('raw_counts_checks170',pd.read_csv(T/'counts_membership_validation.tsv',sep='\t').passed.all());ck('independent_validation4654',pd.read_csv(B/'qc/validation.tsv',sep='\t').passed.all())
pd.DataFrame(checks).to_csv(B/'qc/final_delivery_validation.tsv',sep='\t',index=False);print('FINAL',sum(x['passed'] for x in checks),'/',len(checks));assert all(x['passed'] for x in checks)
