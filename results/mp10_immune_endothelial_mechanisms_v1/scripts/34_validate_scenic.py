import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,pandas as pd,numpy as np
B=pathlib.Path('results/mp10_immune_endothelial_mechanisms_v1');T=B/'tables';O=B/'objects/scenic';checks=[]
def check(name,val,detail=''):
 checks.append(dict(check=name,passed=bool(val),detail=detail))
tfs=set(pathlib.Path('resources/mp10_regulatory_priors/human_tfs.txt').read_text().splitlines());m20=set(pathlib.Path('neighborhood_mp10_calibration_v5/01_score_fidelity/top20.txt').read_text().splitlines())
for lin in ['Epi','T_lineage','Endothelial']:
 genes=set(pd.read_csv(T/f'scenic_genes_{lin}.tsv',sep='\t').gene);m=pd.read_csv(T/f'scenic_membership_{lin}.tsv.gz',sep='\t');check(lin+'_membership_unique',m.barcode.is_unique)
 check(lin+'_training_cap_100',m[m.train].groupby('sample_uid').size().le(100).all())
 if lin=='T_lineage':check('T_training_excludes_NK',not m.analysis_celltype_c2l_final.eq('NK').any())
 adj=pd.read_csv(O/f'{lin}_seed123_adjacencies.tsv',sep='\t');check(lin+'_GRN_features',set(adj.TF).issubset(tfs&genes) and set(adj.target).issubset(genes));check(lin+'_GRN_nonnegative',np.isfinite(adj.importance).all() and adj.importance.ge(0).all())
 ed=pd.read_csv(T/f'scenic_{lin}_regulon_edges.tsv.gz',sep='\t');cov=pd.read_csv(T/f'scenic_{lin}_regulon_coverage.tsv',sep='\t');n=ed.groupby('regulon').target.nunique();check(lin+'_regulon_counts',np.array_equal(n.reindex(cov.regulon).values,cov.n_targets.values));check(lin+'_regulon_features',set(ed.target).issubset(genes));check(lin+'_AUCell_verified',(B/f'qc/aucell_{lin}.json').exists())
 st=pd.read_csv(T/f'scenic_{lin}_connection_seed_stability.tsv.gz',sep='\t');check(lin+'_seed_recovery_semantics',(st.recovered_seed456==(st.importance_seed456.notna()&st.importance_seed456.gt(0))).all())
 check(lin+'_self_targets_not_comparable',not st.loc[st.self_target,'seed_comparison_evaluable'].any());check(lin+'_comparison_has_original_edge',st.loc[st.seed_comparison_evaluable,'importance_seed123'].gt(0).all())
 if lin=='Epi':
  masked=pd.read_csv(T/'scenic_Epi_M20-masked_regulon_edges.tsv.gz',sep='\t');check('Epi_mask_excludes_Top20',not set(masked.target)&m20);check('Epi_mask_exact_edges',len(masked)==len(ed[~ed.target.isin(m20)]))
sc=pd.read_csv(T/'scenic_all_associations.tsv.gz',sep='\t');check('SCENIC_tests_finite',np.isfinite(sc.sample_model_p).all());check('SCENIC_P_range',sc.sample_model_p.between(0,1).all());check('SCENIC_sample_counts_min4',sc.n_recorded_patient_units.ge(4).all());check('SCENIC_unique_tests',not sc.duplicated(['model_id','mask','regulon']).any())
z=pd.DataFrame(checks);z.to_csv(B/'qc/scenic_validation.tsv',sep='\t',index=False);print(z[~z.passed].to_string(index=False));assert z.passed.all();print(len(z),'SCENIC checks passed')
