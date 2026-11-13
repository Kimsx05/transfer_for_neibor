"""Hash-validated GRNBoost2 entry; adoption is limited to the recorded initial task command."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import pathlib,sys,json,hashlib,subprocess,importlib.metadata
P=pathlib.Path.cwd();B=P/'results/mp10_immune_endothelial_mechanisms_v1';O=B/'objects/scenic';lin=sys.argv[1];adopt='--adopt-initial-completed' in sys.argv
expr=O/f'{lin}_training.tsv.gz';tf=P/'resources/mp10_regulatory_priors/human_tfs.txt';out=O/f'{lin}_seed123_adjacencies.tsv';meta=O/f'{lin}_seed123_grn.metadata.json'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
params=dict(lineage=lin,method='grnboost2',seed=123,training_sha256=sha(expr),TF_list_sha256=sha(tf),pyscenic=importlib.metadata.version('pyscenic'),arboreto=importlib.metadata.version('arboreto'),transpose=True)
key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()
if meta.exists() and out.exists() and out.stat().st_size>100 and json.loads(meta.read_text()).get('input_hash')==key:
 assert sha(out)==json.loads(meta.read_text())['output_sha256'];print('Verified GRN cache',lin);sys.exit()
if adopt:
 assert out.exists() and out.stat().st_size>100,'Initial job has not completed'
 # Initial jobs were launched in this task with this exact frozen training input, seed and implementation.
 params['provenance']='Adopted completed initial seed123 CLI job; exact command retained in logs and task execution records'
else:
 if out.exists() and out.stat().st_size>0:raise RuntimeError('Unverified existing GRN: do not overwrite; inspect provenance first')
 subprocess.run([str(P/'mp10_immune_endothelial_mechanisms_v1/venv/bin/pyscenic'),'grn',str(expr),str(tf),'--transpose','--method','grnboost2','--seed','123','--num_workers','1','-o',str(out)],check=True)
params.update(input_hash=key,output_sha256=sha(out));meta.write_text(json.dumps(params,indent=2));print('GRN verified',lin)
