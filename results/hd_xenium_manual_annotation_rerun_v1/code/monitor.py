import os,sys,subprocess,time,pathlib,csv,datetime
O=pathlib.Path(__file__).resolve().parents[1];os.chdir(O)
cg=pathlib.Path('/sys/fs/cgroup')/pathlib.Path('/proc/self/cgroup').read_text().split('::')[-1].strip().lstrip('/')
assert int((cg/'pids.max').read_text())<=10 and int((cg/'memory.max').read_text())<=200000000000
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
for n in 'OPENBLAS_NUM_THREADS MKL_NUM_THREADS OMP_NUM_THREADS OMP_THREAD_LIMIT NUMBA_NUM_THREADS NUMEXPR_NUM_THREADS NUMEXPR_MAX_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS ARROW_NUM_THREADS RAYON_NUM_THREADS POLARS_MAX_THREADS RCPP_PARALLEL_NUM_THREADS LOKY_MAX_CPU_COUNT TF_NUM_INTEROP_THREADS TF_NUM_INTRAOP_THREADS'.split():os.environ[n]='1'
os.environ.update(CUDA_VISIBLE_DEVICES='',PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(O/'logs/mpl'),NUMBA_CACHE_DIR=str(O/'logs/numba'),PYTHONHASHSEED='1729',MPLBACKEND='Agg')
py=O.parents[1]/'highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python'
for arg in sys.argv[1:]:
 script,platform=arg.split(':');stage=script+'_'+platform;start=time.monotonic();threads=rss=memory=0
 with open(O/'logs'/(stage+'.log'),'a') as f:
  proc=subprocess.Popen([py,'-u','code/'+script+'.py',platform],stdout=f,stderr=subprocess.STDOUT)
  while proc.poll() is None:
   tasks=int((cg/'pids.current').read_text());threads=max(threads,tasks);memory=max(memory,int((cg/'memory.current').read_text()));total=0
   for pid in (cg/'cgroup.procs').read_text().split():
    try:total+=int(next(l.split()[1] for l in pathlib.Path('/proc',pid,'status').read_text().splitlines() if l.startswith('VmRSS:')))*1024
    except (OSError,StopIteration):pass
   rss=max(rss,total);time.sleep(.4)
 with open(O/'logs/resources.tsv','a') as f:print(stage,proc.returncode,round(time.monotonic()-start,2),threads,rss,memory,(cg/'pids.peak').read_text().strip(),(cg/'memory.peak').read_text().strip(),sep='\t',file=f)
 if proc.returncode:sys.exit(proc.returncode)

if 'finish:both' in sys.argv[1:]:
 from resources_report import summarize
 summarize(O)
