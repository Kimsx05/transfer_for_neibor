"""Hash-validated resume. Frozen M0 and documented registry addenda are never rewritten.
A completed run is reused only if real inputs, parameters, code and every retained output match.
Partial-run checkpoints are conservative: downstream modifications can trigger safe recomputation.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_THREAD_LIMIT','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import pathlib,hashlib,json,subprocess,sys,datetime,pandas as pd
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v1';S=P/'mp10_hypothesis_followup_v1/scripts';OLD=P/'results/mp10_immune_endothelial_mechanisms_v1';STATE=B/'config/resume_state.json'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
 return h.hexdigest()
def sha_json(x):return hashlib.sha256(json.dumps(x,sort_keys=True).encode()).hexdigest()
def input_signature():
 paths=set();audit=pd.read_csv(B/'tables/input_reuse_audit.tsv',sep='\t');orig=pd.read_csv(B/'tables/original_source_input_manifest.tsv',sep='\t');hashes={}
 for row in pd.concat([audit,orig],ignore_index=True).drop_duplicates('path').itertuples():
  path=pathlib.Path(row.path);actual=sha(path)
  if actual!=row.sha256:raise RuntimeError('Frozen input identity conflict; dependent analysis stopped: '+str(path))
  hashes[str(path)]=actual
 paths.update(OLD.glob('tables/**/*.tsv'));paths.update(OLD.glob('tables/**/*.tsv.gz'))
 paths.update(S.glob('*.R'));paths.update(S.glob('*.py'))
 paths.update(B/'config'/n for n in ['analysis_plan.yaml','M0_FROZEN.sha256','clinical_addendum1.yaml','matched_control_addendum2.json','technical_revision_cluster.md','treg_gene_decomposition_execution.json','plot_config.json','environment_manifest.tsv'])
 paths.update(B/'tables'/n for n in ['model_registry.tsv','model_registry_v1_frozen.tsv','signature_registry.tsv','clinical_covariate.tsv','literature_evidence.tsv'])
 paths.update(B/'logs'/n for n in ['literature_source_records.json','branch_uniprot_annotations.tsv','respiratory_subunits_uniprot.tsv','respiratory_aliases_uniprot.tsv'])
 for path in sorted(paths):
  if path.is_file() and str(path) not in hashes:hashes[str(path)]=sha(path)
 return sha_json(hashes),hashes
def outputs():
 out={}
 for root in ['tables','objects','figures','plot_data','qc','scripts','config']:
  for p in (B/root).rglob('*'):
   if p.is_file() and p!=STATE and p.name not in ['delivery_inventory.json'] and '__pycache__' not in p.parts:out[str(p.relative_to(B))]=sha(p)
 for p in list(B.glob('*.md'))+[B/'module_status.tsv']:out[str(p.relative_to(B))]=sha(p)
 return out
def verified():
 for n in ['validation.tsv','counts_validation.tsv','render_audit.tsv']:
  d=pd.read_csv(B/'qc'/n,sep='\t')
  if not d.passed.all():raise RuntimeError('Unresolved validation: '+n)
STAGES=['01b_registry_addenda.py','02_cycle_pb.R','03_gsea_repair.R','06_branch_enrichment.R','04_scores.R','05_models.R','07_treg_gene_models.R','07b_treg_gsea.R','08_integrate.R','08b_provenance.R','08c_complex_aliases.R','08d_display_tables.R','15_evidence_schema.py','10_report.py','09_replot.R','11_validate_counts.R','11_validate.py','12_render_audit.py','10_report.py','13_delivery.py']
subprocess.run(['taskset','-c','220','Rscript','--vanilla',str(S/'14_environment.R'),'check'],check=True)
print('Verifying frozen input bytes, code and parameters...',flush=True);sig,inp=input_signature();state=json.loads(STATE.read_text()) if STATE.exists() else {}
if '--record-completed' in sys.argv:
 verified();state=dict(status='completed',input_signature=sig,input_hashes=inp,output_hashes=outputs(),recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),note='Recorded only after successful independent counts/numeric/render validation and final delivery; M0 remains immutable');STATE.write_text(json.dumps(state,indent=2));print('Recorded verified completed task',len(state['output_hashes']),'outputs');sys.exit()
if state.get('status')=='completed' and state.get('input_signature')==sig and state.get('output_hashes')==outputs():
 print('All real input, parameter, member/output and code hashes match. Completed run reused.');sys.exit()
print('Hash mismatch or partial task: executing the frozen pipeline, preserving internal input-hashed caches.',flush=True)
for i,name in enumerate(STAGES):
 logfile=B/'logs'/f'resume_{i:02d}_{name}.log';command=['/usr/bin/time','-v','taskset','-c','220']+(['Rscript','--vanilla',str(S/name)] if name.endswith('.R') else [sys.executable,str(S/name)])
 print('Run',name,flush=True)
 with logfile.open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
 STATE.write_text(json.dumps(dict(status='partial',input_signature=sig,last_completed_stage=i,last_completed_name=name),indent=2))
verified();sig,inp=input_signature();STATE.write_text(json.dumps(dict(status='completed',input_signature=sig,input_hashes=inp,output_hashes=outputs(),recorded_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2));subprocess.run([sys.executable,str(S/'13_delivery.py'),'--inventory-only'],check=True);print('Pipeline, validation and hashed delivery completed; no large raw objects copied.')
