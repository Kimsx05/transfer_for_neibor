import os,sys,time,subprocess,hashlib,json
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
R=Path(__file__).resolve().parents[1];py=R.parent/'hd_epi_niche_multimethod_v2/.venv/bin/python';force='--force'in sys.argv
if os.environ.get('HD_ALLOW_DAYTIME')!='1':
 while 7<=datetime.now(ZoneInfo('Asia/Shanghai')).hour<23:
  print('Waiting for project compute window; not completed',flush=True);time.sleep(30)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','OMP_THREAD_LIMIT']:os.environ[k]='1'
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:4]));(R/'logs').mkdir(exist_ok=True);cp=R/'.checkpoints.json';cache=json.loads(cp.read_text())if cp.exists()else{};dirty=force
steps=[('01_freeze_baseline.py','tables/baseline_reproduction.tsv'),('02_mixtures.py','tables/expanded_mixtures.parquet'),('03_paired_models.py','tables/paired_effects.tsv'),('04_spatial_support.py','tables/source_coverage_summary.tsv'),('05_calibrate_evidence.py','evidence_decision.yaml'),('06_diagnostics.py','tables/new_test_registry.tsv'),('07_roi_plots.py','tables/ROI_evidence_notes.tsv'),('08_common_counts.py','tables/common_support_paired_effects.tsv'),('08b_count_accounting.py','tables/common_support_count_accounting.tsv'),('09_summary_plots.py','figures/mixture_detectability.png'),('10_validate.py','tables/validation.tsv')]
def signature(name):
 # Only small implementation/config inputs, never full expression hashes.
 z=(R/'code'/name).read_bytes()+(R/'code/fu_common.py').read_bytes()
 for p in [R/'frozen_config.yaml',R.parent/'hd_treg_mp10_spatial_all_v1/frozen_config.yaml',R.parent/'hd_treg_mp10_spatial_all_v1/tables/input_manifest.tsv']:
  if p.exists():z+=p.read_bytes()
 return hashlib.sha256(z).hexdigest()
for name,out in steps:
 sig=signature(name)
 if not dirty and cache.get(name)==sig and(R/out).exists():print('RESUME',name,flush=True);continue
 dirty=True
 if name=='02_mixtures.py':os.environ['FORCE_MIX']='1'
 print('RUN',name,flush=True)
 with open(R/'logs'/('rerun_'+name+'.log'),'w')as log:subprocess.run(['/usr/bin/time','-v',str(py),str(R/'code'/name)],stdout=log,stderr=subprocess.STDOUT,check=True)
 cache[name]=signature(name);cp.write_text(json.dumps(cache,indent=2))
subprocess.run([str(py),str(R/'code/11_report.py')],check=True);subprocess.run([str(py),str(R/'code/12_package.py')],check=True)
