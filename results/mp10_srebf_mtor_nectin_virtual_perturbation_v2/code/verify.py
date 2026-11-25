from common import *
import platform,re
from PIL import Image
from datetime import datetime,timezone,timedelta
T=Path('tables');F=rd('FIGURE_INDEX.tsv');checks=[]
def checked(description,truth):
 assert truth,description
 checks.append(description)
a=rd('qc/input.tsv').iloc[0];checked('exact 8000 cells, 38 patients, 39 specimens, full/clean 2233/2194',(a.n_cells,a.n_patients,a.n_specimens,a.full,a.clean)==(8000,38,39,2233,2194))
co=rd(T/'celloracle_job_status.tsv');checked('two seeds plus five omissions with ten original clusters each',len(co)==70 and (co.status=='complete').sum()==63 and (co[co.status=='NE'].state=='C9').all())
sj=rd(T/'knk_job_status.tsv');checked('38 fitted WT and one prespecified small cluster NE',len(sj)==39 and (sj.status=='complete').sum()==38)
mat=rd(T/'control_matching.tsv');checked('five frozen matches per source',(mat.groupby('ko').size()==5).all())
for marker in Path('cache/jobs').glob('*/*.json'):
 j=json.loads(marker.read_text());assert all(sha(f)==h for f,h in j['outputs'].items())
for marker in Path('cache/oracle').glob('*/*/complete.json'):
 j=json.loads(marker.read_text());assert all(sha(marker.parent/f)==h for f,h in j['outputs'].items())
checks.append('all scTenifold and CellOracle job checkpoint SHA-256 verified')
for ix in [T/'knk_all_gene_index.tsv',T/'celloracle_cluster_index.tsv']:
 for r in rd(ix).itertuples():assert sha(r.path)==r.sha256
checks.append('signed full-gene file indices verified')
cross=rd(T/'cross_method_summary.tsv');x=cross[(cross.analysis=='main')&(cross.gene=='NECTIN2')];checked('NECTIN2 present for all 5 and non-TF method NE',set(x.ko)==set(K) and x[x.ko.isin(K[2:])].signed_delta.isna().all())
z=rd(T/'celloracle_all_genes.tsv.gz');checked('CellOracle supplies no native P/q and excludes source from ranking',z['p.adj'].isna().all() and z['p.value'].isna().all() and z[z.is_KO]['rank'].isna().all())
checked('all three propagation steps and a randomized reference',set(z[(z.ko=='SREBF1')&(z.analysis=='main')].steps)=={2,3,4} and 'randomized' in z.graph.unique())
q=rd(T/'celloracle_simulation_qc.tsv');checked('no-op max abs under 1e-10',q.no_perturbation_max_abs.max()<1e-10)
for r in F.itertuples():
 assert Path(r.pdf).read_bytes()[:4]==b'%PDF' and '<svg' in Path(r.svg).read_text()[:3000] and Path(r.plot_data).is_file()
 with Image.open(r.png) as im:assert min(im.size)>700;im.verify()
checked('8 plot sets have validated PDF, SVG, PNG and source data',len(F)==8)
resources={n:(cg/n).read_text().strip() for n in ['memory.max','memory.high','memory.peak','memory.events','cpu.max','cpu.stat']}
resources['cpuset_controller_effective']=(cg/'cpuset.cpus.effective').read_text().strip() if (cg/'cpuset.cpus.effective').is_file() else None
resources['scoped_taskset_affinity']=sorted(os.sched_getaffinity(0));assert 0<len(resources['scoped_taskset_affinity'])<=C['threads_total']
sam=rd('logs/cgroup_usage.tsv');resources.update(observation_start_local=datetime.fromtimestamp(sam.epoch.iloc[0],timezone(timedelta(hours=8))).isoformat(),observation_last_local=datetime.fromtimestamp(sam.epoch.iloc[-1],timezone(timedelta(hours=8))).isoformat(),maximum_telemetry_gap_seconds=float(sam.epoch.diff().max()),sample_window_seconds=round(sam.epoch.iloc[-1]-sam.epoch.iloc[0],2),sampled_disk_peak_bytes=int(sam.disk_bytes_periodic.max()),result_tree_bytes=int(subprocess.check_output(['du','-sb','.']).split()[0]),budget_decimal_bytes=C['memory_bytes'])
checked('task slice maximum memory and peak within budget',int(resources['memory.max'])<=C['memory_bytes'] and int(resources['memory.peak'])<=C['memory_bytes'])
resources['early_setup_scope']='Initial dependency acquisition and report editing were outside the model measurement window.'
Path('qc/resources.json').write_text(json.dumps(resources,indent=2)+'\n')
times=[]
for f in sorted(Path('logs').glob('*.usage')):
 txt=f.read_text()
 def attr(k):
  for l in txt.splitlines():
   if l.strip().startswith(k):return l.strip()[len(k):].strip()
 times.append(dict(job=f.stem,elapsed=attr('Elapsed (wall clock) time (h:mm:ss or m:ss):'),peak_rss_kib=attr('Maximum resident set size (kbytes):'),user_seconds=attr('User time (seconds):'),system_seconds=attr('System time (seconds):'),exit_status=attr('Exit status:')))
save(pd.DataFrame(times),'qc/stage_resources.tsv')
checked('no source inputs or v1 artifacts overwritten',all(sha(f.path)==f.sha256 for f in rd('config/input_fingerprints.tsv').itertuples()))
Path('qc/validation.json').write_text(json.dumps(dict(checks=checks,status='passed'),ensure_ascii=False,indent=2)+'\n')
cache=[p for p in Path('cache').rglob('*') if p.is_file() and not any(n in p.parts for n in ('co310','oracle_env','__pycache__'))]
files=cache+[p for name in ['code','config','qc','tables','figures','plot_data'] for p in Path(name).rglob('*') if p.is_file() and '__pycache__' not in p.parts]+[Path(t) for t in ['METHODS.md','RESULTS.md','FIGURE_INDEX.tsv','run_all','replot','scope']]
# Task logs remain available in logs/. The live cgroup sampler is recorded as a measured window in qc/resources.json.
files+= [p for p in Path('logs').glob('*') if p.is_file() and p.name not in ['cgroup_usage.tsv','verify.log'] and p.stat().st_size>0]
out=dict(task='MP10–SREBF–mTOR–NECTIN virtual perturbation v2',status='complete_with_prespecified_NE',n_sc_WT=38,n_oracle_cluster_fits=63,methods='METHODS.md',results='RESULTS.md',tasks='tables/task_status.tsv',resources='qc/resources.json',unhashed=['tmp/','cache/co310/ (environment lock in config)','logs/cgroup_usage.tsv (live telemetry)','logs/verify.log (written by verifier)'],files=[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(files))]);new=Path('manifest.new');new.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n');os.replace(new,'manifest.json')
print('PASS',len(checks),'checks;',len(out['files']),'artifacts',flush=True)

