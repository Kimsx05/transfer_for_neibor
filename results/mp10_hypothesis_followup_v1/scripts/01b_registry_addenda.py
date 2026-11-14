"""Reconstruct the documented conditional registry additions; default verifies, --rebuild writes.
No effect statistics are read. The initial M0 registry remains immutable.
"""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,pandas as pd
B=pathlib.Path('results/mp10_hypothesis_followup_v1');base=pd.read_csv(B/'tables/model_registry_v1_frozen.tsv',sep='\t');cl=pd.read_csv(B/'tables/clinical_field_audit.tsv',sep='\t');assert cl.Type.isin(['MIBC','NMIBC']).all() and sum(cl.groupby('dataset').Type.nunique()>1)>=2
clinical=base[(base.setting=='main')&base.FDR_family.isin(['bridge_primary','bridge_noncycling','bridge_matched_all'])].copy();clinical['setting']='clinical_type';clinical['covariate']='clinical_type';clinical['role']='sensitivity';clinical['FDR_family']+='__clinical_type';clinical['model_id']=['F%05d'%i for i in range(len(base)+1,len(base)+len(clinical)+1)]
control=base[(base.setting=='main')&((base.FDR_family=='bridge_matched_all')|((base.FDR_family=='treg_decomposition')&(base.cohort=='NC_matched')))].copy();control['setting']='receiver_MT';control['covariate']='receiver_MT';control['role']='sensitivity';control['FDR_family']+='__receiver_MT';control.loc[control.y.eq('G2M_CHECKPOINT_Full'),'FDR_family']='treg_matched_G2M_control__receiver_MT';control['planned_n']=control.groupby('FDR_family').model_id.transform('size');control['model_id']=['F%05d'%i for i in range(len(base)+len(clinical)+1,len(base)+len(clinical)+len(control)+1)]
allr=pd.concat([base,clinical,control],ignore_index=True)
if '--rebuild' in sys.argv:
 allr.to_csv(B/'tables/model_registry.tsv',sep='\t',index=False);clinical.to_csv(B/'tables/model_registry_clinical_addendum.tsv',sep='\t',index=False);control.to_csv(B/'tables/model_registry_matched_control_addendum.tsv',sep='\t',index=False)
else:pd.testing.assert_frame_equal(allr,pd.read_csv(B/'tables/model_registry.tsv',sep='\t'),check_dtype=False)
print('Registry additions reproducible:',len(allr))
