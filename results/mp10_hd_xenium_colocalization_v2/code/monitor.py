import sys,os,time,subprocess,fcntl,pathlib,datetime
base=pathlib.Path(__file__).resolve().parents[1]
root=pathlib.Path(os.environ.get('MP_V2_OUTPUT',base));os.chdir(root)
lock=open(base/'logs/task.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
scope=pathlib.Path('/sys/fs/cgroup')/pathlib.Path('/proc/self/cgroup').read_text().strip().split('::')[1].lstrip('/')
assert int((scope/'pids.max').read_text())<=7 and int((scope/'memory.max').read_text())<=165000000000
os.sched_setaffinity(0,{min(os.sched_getaffinity(0))})
for k in 'OPENBLAS_NUM_THREADS MKL_NUM_THREADS OMP_NUM_THREADS OMP_THREAD_LIMIT NUMBA_NUM_THREADS NUMEXPR_NUM_THREADS NUMEXPR_MAX_THREADS VECLIB_MAXIMUM_THREADS RCPP_PARALLEL_NUM_THREADS BLIS_NUM_THREADS ARROW_NUM_THREADS R_DATATABLE_NUM_THREADS LOKY_MAX_CPU_COUNT POLARS_MAX_THREADS RAYON_NUM_THREADS TBB_NUM_THREADS'.split():os.environ[k]='1'
os.environ.update(PYTHONDONTWRITEBYTECODE='1',MPLCONFIGDIR=str(root/'logs/mpl'),NUMBA_CACHE_DIR=str(root/'logs/numba'),CUDA_VISIBLE_DEVICES='',R_LIBS=str(root/'code/Rlib')+':'+str(base.parents[1]/'highres_VisiumHD/analysis/hd_bin8_comparison_v1/Rlib')+':/refdir/Rlib')
py=base.parents[1]/'highres_VisiumHD/analysis/hd_epi_niche_multimethod_v2/.venv/bin/python'
cmd=[str(py),'-u',*sys.argv[1:]] if sys.argv[1].endswith('.py') else sys.argv[1:]
stage='_'.join(pathlib.Path(a).name for a in sys.argv[1:])[:120]
start=time.time()
with open(root/'logs'/(stage+'.log'),'a',buffering=1) as f:
 p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT)
 while p.poll() is None:time.sleep(.3)
peak=(scope/'memory.peak').read_text().strip()
tasks=(scope/'pids.peak').read_text().strip() if (scope/'pids.peak').exists() else 'NA'
with open(root/'audit/resources.tsv','a') as f:print(datetime.datetime.now().isoformat(),stage,p.returncode,round(time.time()-start,2),peak,tasks,sep='\t',file=f)
if p.returncode:print((root/'logs'/(stage+'.log')).read_text()[-2500:])
sys.exit(p.returncode)
