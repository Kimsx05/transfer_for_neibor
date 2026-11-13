import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pathlib,hashlib,json,subprocess,datetime,importlib.metadata
import pandas as pd
P=pathlib.Path.cwd();B=P/'results/mp10_immune_endothelial_mechanisms_v1';C=P/'mp10_immune_endothelial_mechanisms_v1';R=P/'resources/mp10_regulatory_priors'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1048576),b''):h.update(x)
 return h.hexdigest()
def rd(n):return pd.read_csv(B/'tables'/n,sep='\t')
checks=rd('validation_checks.tsv');assert checks['pass'].all(),'Validation failures remain'
scenic_checks=pd.read_csv(B/'qc/scenic_validation.tsv',sep='\t');assert scenic_checks.passed.all(),'SCENIC validation failures remain'
render=pd.read_csv(B/'qc/render_audit.tsv',sep='\t');assert render.exists.all() and render.direct_data_exists.all() and render.render_ok.all()
assert not render.outside_text.fillna('').str.len().gt(0).any(),'Clipped PDF text remains'
assert render.loc[render.file.str.endswith('.pdf'),'raster_images'].eq(0).all(),'PDF contains raster plot elements'
assert render.loc[render.file.str.endswith('.png'),'dpi_x'].between(299,301).all()
for lin in ['Epi','T_lineage','Endothelial']:
 for a in json.loads((B/f'qc/aucell_{lin}.json').read_text()):assert a['exact_membership_order'] and a['finite_unit_interval']

assert json.loads((B/'qc/patient_relationships.json').read_text())['no_joint_a_b_in_any_model']
status=rd('model_status.tsv');plan=rd('model_plan.tsv');modules=[]
def add(module,state,scope,reason='',n=None,datasets=None):modules.append(dict(module=module,status=state,actual_scope=scope,reason=reason,n_recorded_patient_units=n,n_datasets=datasets))
add('input_and_pseudobulk','completed','40 original specimens; 39 recorded patient units; full counts and original normalized RNA; exact final-reference membership',n=39,datasets=5)
for r in status[status.primary].itertuples():add('DE_'+r.target,r.status,f'{r.model_id}; n_cells >=50; own eligible sample set',str(r.reason) if pd.notna(r.reason) else '',r.n_recorded_patient_units,r.n_datasets)
add('DE_prespecified_sensitivity','completed',f'{(status.status=="completed").sum()} completed and {(status.status=="not_evaluable").sum()} explicitly NE models of {len(status)} planned; all legal settings evaluated')
for name,f in [('Hallmark_Reactome','gsea_results.tsv.gz'),('PROGENy_CollecTRI_DoRothEA','activity_associations.tsv.gz'),('abundance','abundance_associations.tsv'),('Epi_composition_quality','epi_within_cluster_candidate_changes.tsv')]:assert (B/'tables'/f).exists();add(name,'completed','Main and applicable prespecified sensitivities; full positive/negative results retained')
ns=rd('nichenet_analysis_status.tsv');main=ns.merge(plan[['model_id','primary']],on='model_id');main=main[main.primary&(main.detection_threshold==.1)]
for r in main.itertuples():add(f'NicheNet_{r.receiver}_{r.model_id}_{r.mask}_{r.target_direction}',r.status,f'{r.direction}; {r.n_targets} targets / {r.n_background} testable background; native and focused ranks',str(r.reason) if pd.notna(r.reason) else '',r.n_specimens,int(status.loc[status.model_id.eq(r.model_id),'n_datasets'].iloc[0]))
add('NicheNet_sensitivity','completed','Full/M20-masked, 10%/5% detection, positive/negative/union and legal patient/dataset/threshold/exposure models; sparse targets flagged')
for lin in ['Epi','T_lineage','Endothelial']:
 for f in [f'scenic_{lin}_regulon_edges.tsv.gz',f'scenic_{lin}_associations.tsv.gz',f'scenic_{lin}_RSS.tsv.gz',f'scenic_{lin}_connection_seed_stability.tsv.gz']:assert (B/'tables'/f).exists(),f
 meta=json.loads((B/'objects/scenic'/f'{lin}_regulons.metadata.json').read_text());training=rd('scenic_training_audit.tsv');tr=training[training.lineage==lin].iloc[0]
 add('pySCENIC_'+lin,'completed',f'{tr.n_training} training cells; {tr.n_full} full scored cells; {meta["n_regulons"]} motif-pruned regulons; seed456 key-connection check','Second seed is not a second full network/AUCell experiment')
seedcov=rd('scenic_second_seed_coverage.tsv')
for lin,z in seedcov.groupby('lineage'):
 n=int(z.n_evaluable_GRN_edges.sum());add('pySCENIC_key_connection_stability_'+lin,'completed' if n else 'not_evaluable',f'{n} comparable non-self GRN edges in prespecified TFs; individual TF coverage retained','Preselected TFs lack motif-pruned regulons' if not n else 'Key connections only; not a full second network')
assert (B/'tables/scenic_Epi_M20-masked_associations.tsv.gz').exists();add('pySCENIC_Epi_M20_mask','completed','Frozen network with Top20 removed only from downstream regulon/scoring universe')
add('matched_scRNA_CNV_sensitivity','not_evaluable','No new CNV inferred','No verified matched existing per-cell CNV used in this task; optional sensitivity unavailable')
add('resource_runtime_audit','completed','Runtime thread configuration audited and corrected','Early global Rprofile override violated requested cap; see incident log')
add('evidence_integration_and_rendering','completed','Full unweighted candidate table, direct plotting data, all figure parts and render audit')
pd.DataFrame(modules).to_csv(B/'module_status.tsv',sep='\t',index=False)
cachefiles=[B/'objects/pseudobulk.qs',B/'tables/pseudobulk_membership.tsv.gz',B/'tables/analysis_cell_metadata.tsv.gz',B/'config/run_config.json']+[B/'objects/scenic'/f'{l}_training.tsv.gz' for l in ['Epi','T_lineage','Endothelial']]
(B/'config/delivered_cache_hashes.json').write_text(json.dumps({str(p.relative_to(B)):sha(p) for p in cachefiles},indent=2))
resources=[]
for p in sorted(R.iterdir()):
 if p.is_file():resources.append(dict(file=str(p),bytes=p.stat().st_size,mtime=datetime.datetime.fromtimestamp(p.stat().st_mtime).isoformat(),sha256=sha(p)))
pd.DataFrame(resources).to_csv(B/'config/resource_manifest.tsv',sep='\t',index=False)
code=[dict(file=str(p),sha256=sha(p)) for p in sorted((C/'scripts').glob('*')) if p.is_file()];pd.DataFrame(code).to_csv(B/'config/code_manifest.tsv',sep='\t',index=False)
(B/'config/python_freeze.txt').write_text(subprocess.check_output([str(C/'venv/bin/python'),'-m','pip','freeze'],text=True))
r='source("mp10_immune_endothelial_mechanisms_v1/scripts/common.R"); for(p in c("edgeR","limma","Seurat","qs","fgsea","nichenetr","decoupleR","RhpcBLASctl")){cat(p,as.character(packageVersion(p)),"\\n")}; print(sessionInfo()); cat("BLAS threads",RhpcBLASctl::blas_get_num_procs(),"OMP max",RhpcBLASctl::omp_get_max_threads(),"\\n")'
(B/'logs/sessionInfo_final.txt').write_text(subprocess.check_output(['Rscript','--vanilla','-e',r],text=True))
(B/'config/completion.json').write_text(json.dumps(dict(completed_at=datetime.datetime.now().isoformat(),narrative_finalized=('工作报告正在补齐' not in (B/'RESULTS.md').read_text()),validation_checks=len(checks),scenic_validation_checks=len(scenic_checks),figures=int(render.figure_id.nunique()),figure_files=len(render),all_passed=True,models_completed=int((status.status=='completed').sum()),models_NE=int((status.status=='not_evaluable').sum()),large_originals_copied=False,HD_evidence='pending',resource_constraint_deviation='See logs/resource_constraint_incident.json; early Rprofile override corrected, no concealment'),indent=2))
docs=[B/x for x in ['README.md','METHODS.md','RESULTS.md','FIGURE_INDEX.md','module_status.tsv']]
pd.DataFrame([dict(file=str(p),sha256=sha(p)) for p in docs]).to_csv(B/'config/document_manifest.tsv',sep='\t',index=False)
print('DELIVERY MANIFEST COMPLETE')
