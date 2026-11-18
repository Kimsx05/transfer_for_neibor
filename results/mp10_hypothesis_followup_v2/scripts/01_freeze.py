import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import pathlib,json,hashlib,datetime,pandas as pd
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';O=P/'results/mp10_immune_endothelial_mechanisms_v1';V=P/'results/mp10_hypothesis_followup_v1';R=P/'resources/mp10_regulatory_priors'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(8*1024**2),b''):h.update(x)
 return h.hexdigest()
if (B/'config/M0_FROZEN.json').exists():raise SystemExit('M0 already frozen; do not overwrite')
sig=pd.read_csv(B/'tables/signature_registry.tsv',sep='\t');rows=[]
def add(kind,setting,family,x,y,state='H05',mask='Full',min_cells=50,**kw):
 rows.append(dict(kind=kind,setting=setting,FDR_family=family,x=x,y=y,state=state,mask=mask,min_cells=min_cells,patient='replacement' if setting=='replace_b' else 'main',primary=setting=='main' and mask=='Full',**kw))
sets=sig.loc[sig.role.eq('formal_process_primary'),'signature'].tolist()
for state in ['H05','H10','H50']:
 for mask in ['Full','M20']:
  for setting in ['main','min30','replace_b','common_H05','quality_RNA','composition_cluster0']+['dataset_'+s for s in ['CNP0000460','GSE135337','GSE326225','InhouseData','PRJNA662018']]:
   for s in sets:add('paired_score',setting,f'paired_{state}_{mask}_{setting}',state,s,state,mask,30 if setting=='min30' else 50)
  for s in sets:add('gsea','main',f'gsea_{state}_{mask}',state,s,state,mask)
for state in ['H05','H10','H50']:
 for mask in ['Full','M20']:
  for chol in ['CHOL_SYNTH','BLOCH']:
   for setting in ['main','min30','replace_b','common_H05','quality_RNA','composition_cluster0']:
    add('direct_D',setting,f'direct_D_{state}_{mask}_{setting}',chol,'FA_ACYL_COA',state,mask,30 if setting=='min30' else 50)
routes=[('cross_primary','CHOL_EFFLUX','CHOL_IMPORT','Epi_all','Treg'),('cross_primary','FA_ACYL_COA','LCFA_IMPORT','Epi_all','Treg'),('EC_exploratory','HALLMARK_MTORC1_SIGNALING','LCFA_IMPORT','Capillary EC','Capillary EC'),('EC_exploratory','HALLMARK_MTORC1_SIGNALING','LCFA_TRANSPORT','Capillary EC','Capillary EC'),('EC_exploratory','LCFA_TRANSPORT','LCFA_IMPORT','Capillary EC','Treg'),('EC_exploratory','FA_ACYL_COA','LCFA_IMPORT','Epi_all','Capillary EC'),('EC_exploratory','CHOL_EFFLUX','CHOL_IMPORT','Epi_all','Capillary EC')]
for x in ['CHOL_EFFLUX','FA_ACYL_COA']:
 for y in ['HALLMARK_OXIDATIVE_PHOSPHORYLATION','HALLMARK_G2M_CHECKPOINT','cycling_fraction']:routes.append(('Treg_secondary',x,y,'Epi_all','Treg'))
settings=['main','min30','replace_b','delta_H05_L50','M20','H10','H50','Epi_all','multi_specimen_datasets','quality_sender_RNA','quality_receiver_RNA','no_shared_genes']+['LODO_'+s for s in ['CNP0000460','GSE135337','GSE326225','InhouseData','PRJNA662018']]
for setting in settings:
 for fam,x,y,sender,receiver in routes:
  if setting in ['H10','H50','delta_H05_L50','Epi_all'] and sender!='Epi_all':continue
  if setting=='no_shared_genes' and sender!=receiver:continue
  add('cross',setting,f'{fam}_{setting}',x,y,'H05','M20' if setting=='M20' else 'Full',30 if setting=='min30' else 50,sender=sender,receiver=receiver)
reg=pd.DataFrame(rows);reg.insert(0,'model_id',[f'V2_{i:04d}' for i in range(1,len(reg)+1)]);reg['planned_n']=reg.groupby('FDR_family').model_id.transform('size');reg['statistical_unit']='recorded_patient_unit; specimen; a main / b replacement';reg.to_csv(B/'tables/model_registry.tsv',sep='\t',index=False)
plan=dict(task='mp10_hypothesis_followup_v2',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),design_status='hypothesis-directed internal analysis after historical discovery',primary_High='H05',Low='frozen L50',High_sensitivities=['H10','H50'],membership='reuse exact original pseudobulk membership; no rounded-score reassignment',min_cells=[50,30],min_units=4,min_residual_df=3,thread_total_max=4,internal_threads=1,patient='main T4463209a; replace with b sensitivity; recorded/provisional others',paired_counts='reuse edgeR robust QL ~ specimen + state; original global gene FDR preserved',paired_scores='OLS on specimen High-Low differences; HC3 t/CI, intercept test; no duplicate receiver',cross='Y ~ dataset + standardized X; OLS HC3 with residual-df t; quasi-binomial for cycling proportion; h>=1-1e-10 NE',score='reuse v1 full-transcriptome lineage TMM logCPM; ssGSEA GSVA 2.4.9 alpha .25 normalize FALSE; >=5 genes and >=50% expected members',scale='per-signature fixed mean/SD in main-a Epi_all>=50 (sender) or native receiver>=50; same fixed reference across all states/subsets; never center deltas',mask='M20 removed from downstream signature AND ranking universe only, RNA normalization unchanged; receiver M20 only sensitivity',gsea=dict(minSize=15,maxSize=500,seed=123,nPermSimple=10000,eps=1e-10,scoreType='std',numeric_budget='uniform 100000 for nonfinite estimates only'),FDR='BH with planned family n; NE remains NA',direct_D='standardized CHOL_SYNTH change minus FA_ACYL_COA change; BLOCH secondary; Full/M20; same paired samples',matched_gene_resampling=dict(seed=123,repeats=200,method='equal-size nearest-neighbor matching without replacement between measurable cholesterol and FA genes by main-a Epi_all mean logCPM and median detected fraction standardized across whole fixed Epi universe; draw k=min sizes, randomly permute assignment order; fixed 200 samples, do not select by effect; identical scoring universe, per-signature reference fixed before each delta; report robustness distribution not biological-replicate P'),scenic='exact frozen seed123 training/adjacency/context; min_genes20; original missing_fraction>=.20 excluded, including exactly80%; same annotation and NES>=3; first diagnose all SREBF1/2 modules; strict unique HGNC exact/previous-symbol correction only with retained denominator, no ambiguous mapping; no new GRNBoost2',scenic_alternative_reference='no alternate ranking preselected because no independent compatible alternative installed; original library limitations retained, do not select resource to rescue TF',EC_QC='original labels; endothelial PECAM1,VWF,CDH5,EMCN,KDR,RAMP2 versus EPCAM,KRT8,KRT18,KRT19; single-cell coexpression and contribution descriptive; no ambient correction claimed',mTOR_lipid_diagnostic='union of frozen 13 formal lipid processes; label MP10 overlap separately; remainder other Hallmark response',no_FFA_export_signature='no accepted generic free-fatty-acid secretion assay: explicit gap, lipolysis not substituted',clinical='no additional clinical scan',NE_conditions=['profile<cell threshold','<4 recorded units','rank deficient','df<3','zero exposure variance','signature coverage below fixed rule','HC3 h=1','quasi separation/nonconvergence','GSEA size outside15..500'])
(B/'config/analysis_plan.json').write_text(json.dumps(plan,indent=2,ensure_ascii=False));(B/'config/analysis_plan.yaml').write_text(json.dumps(plan,indent=2,ensure_ascii=False))
paths=[O/'objects/pseudobulk.qs',V/'objects/scoring_matrices.qs',O/'tables/pseudobulk_membership.tsv.gz',O/'tables/analysis_cell_metadata.tsv.gz',O/'tables/patient_relationship_audit.tsv',O/'tables/model_plan.tsv',O/'tables/model_membership.tsv.gz',O/'objects/scenic/Epi_training.tsv.gz',O/'objects/scenic/Epi_seed123_adjacencies.tsv',O/'objects/scenic/Epi_regulons.metadata.json',O/'objects/scenic/Epi_regulons.pkl',O/'objects/scenic/Epi_motif_enrichment.csv',R/'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather',R/'motifs-v9-nr.hgnc-m0.001-o0.0.tbl',R/'reactome.tsv.gz',R/'hallmark.tsv.gz',pathlib.Path('/home/data/t070721/.cache/R/msigdbr/msigdb.2026.1.Hs.C5.rds')]
a=[]
for p in paths:
 assert p.exists(),p
 st=p.stat();a.append(dict(path=str(p),bytes=st.st_size,mtime_utc=datetime.datetime.fromtimestamp(st.st_mtime,datetime.timezone.utc).isoformat(),sha256=sha(p)))
pd.DataFrame(a).to_csv(B/'tables/input_manifest.tsv',sep='\t',index=False)
(B/'config/M0_FROZEN.json').write_text(json.dumps({str(p.relative_to(B)):sha(p) for p in [B/'config/analysis_plan.json',B/'tables/model_registry.tsv',B/'tables/signature_registry.tsv',B/'tables/official_gene_sets.tsv.gz',B/'tables/input_manifest.tsv']},indent=2))
print('FROZEN',len(reg),'models',reg.groupby(['kind']).size().to_dict())
