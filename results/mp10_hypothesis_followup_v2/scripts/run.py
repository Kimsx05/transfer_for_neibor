#!/usr/bin/env python3
"""Sequential frozen-design runner. Refuse changed inputs/outputs; resume missing outputs."""
import os,sys,pathlib,json,hashlib,subprocess,argparse,datetime,shutil
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS','RCPP_PARALLEL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS']:os.environ[k]='1'
P=pathlib.Path.cwd();B=P/'results/mp10_hypothesis_followup_v2';S=P/'mp10_hypothesis_followup_v2/scripts';PY=P/'mp10_immune_endothelial_mechanisms_v1/venv/bin/python';CHECK=B/'config/completion_checkpoint.json'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(8*1024**2),b''):h.update(c)
 return h.hexdigest()
def verify(entries,allow_missing=False):
 missing=[]
 for path,digest in entries.items():
  p=pathlib.Path(path)
  if not p.exists():
   if allow_missing:missing.append(path);continue
   raise RuntimeError('Missing dependency: '+path)
  if sha(p)!=digest:raise RuntimeError('Hash mismatch; refusing incompatible overwrite: '+path)
 return missing
def execute(name):
 cmd=['Rscript','--vanilla',str((S/name).relative_to(P))] if name.endswith('.R') else [str(PY),str(S/name)]
 allowed=sorted(os.sched_getaffinity(0));cpu=220 if 220 in allowed else allowed[0]
 cmd=['taskset','-c',str(cpu)]+cmd
 with open(B/'logs'/('resume_'+name.replace('.','_')+'.log'),'w') as log:
  subprocess.run(cmd,cwd=P,env=os.environ.copy(),stdout=log,stderr=subprocess.STDOUT,check=True)
 print('completed',name,flush=True)
def seal():
 import pandas as pd
 inp=pd.read_csv(B/'tables/input_reuse_audit.tsv',sep='\t');sources=dict(zip(inp.path,inp.sha256))
 controls={str(p):sha(p) for p in list(S.glob('*'))+list((B/'config').glob('*')) if p.is_file() and p.name not in ['completion_checkpoint.json','cache_manifest.json']}
 outputs={str(p):sha(p) for root in ['tables','objects','figures','plot_data','qc','scripts'] for p in (B/root).rglob('*') if p.is_file()}
 outputs.update({str(p):sha(p) for p in B.glob('*.md')});outputs[str(B/'module_status.tsv')]=sha(B/'module_status.tsv')
 key=hashlib.sha256(json.dumps([sources,controls],sort_keys=True).encode()).hexdigest()
 CHECK.write_text(json.dumps(dict(task='mp10_hypothesis_followup_v2',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),dependency_hash=key,inputs=sources,controls=controls,outputs=outputs),indent=2))
 caches={p:h for p,h in outputs.items() if '/objects/' in p}
 (B/'config/cache_manifest.json').write_text(json.dumps(dict(dependency_hash=key,inputs=sources,parameters_and_code=controls,cache_outputs=caches,reuse_policy='all input/code/config and cached output hashes must match; filenames alone never authorize reuse'),indent=2))
 print('sealed',len(outputs),'artifacts;',len(sources),'input dependencies',flush=True)
def main():
 ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True)
 for x in ['resume','recompute','replot','seal']:g.add_argument('--'+x,action='store_true')
 a=ap.parse_args()
 if a.seal:seal();return
 if not CHECK.exists():raise RuntimeError('Missing completion checkpoint; initialize only after full provenance validation, do not infer from filenames.')
 cp=json.loads(CHECK.read_text())
 if a.replot:
  # This branch deliberately never opens original expression/Seurat/resources.
  verify({p:h for p,h in cp['outputs'].items() if '/tables/' in p})
  verify({p:h for p,h in cp['controls'].items() if pathlib.Path(p).name in ['11_replot.R','12_render_audit.py','theme.R','export.R','plot_config.json']})
  execute('11_replot.R');execute('12_render_audit.py');execute('16_final_audit.py');seal();return
 verify(cp['inputs']);verify(cp['controls']);missing=verify(cp['outputs'],allow_missing=True)
 if a.resume and not missing:print('Verified compatible completed result; all cached inputs, parameters, members and outputs match. No inference rerun.');return
 if a.resume and missing and all(any('/'+d+'/' in p for d in ['figures','plot_data','qc']) for p in missing):
  execute('11_replot.R');execute('12_render_audit.py');execute('13_validate.py');execute('16_final_audit.py');seal();return
 # Missing files only: deterministic owned result regeneration; modified existing files never overwritten.
 print('Frozen pipeline regeneration; missing files:',len(missing),flush=True)
 for name in ['02_scenic_audit.py','02b_mapping_repair.py','02c_retained_motif_details.py','03_scores_audit.R','04_EC_cell_QC.R','05_models.R','06_formal_gsea.R','07_matched_resampling.R','08_spatial_audit.py','08b_expression_readouts.R','08c_annotation_summary.py','08e_EC_masked_gsea.R','09_aucell_recovered.R','09c_scenic_effects.R','14_delivery.py','11_replot.R','13_validate.py','13b_counts_validate.R','12_render_audit.py','16_final_audit.py']:
  execute(name)
 for p in S.iterdir():
  if p.is_file():shutil.copy2(p,B/'scripts'/p.name)
 seal()
if __name__=='__main__':main()
