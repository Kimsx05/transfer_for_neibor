from lib import *
import re,sys,threadpoolctl
from PIL import Image
import pyarrow.parquet as par
p=O/'classification';lock=json.loads((p/'classification_lock.json').read_text())
assert sha(O/'provisional_classification.json')==lock['config_sha256'] and sha(p/'section_classification.tsv')==lock['class_sha256']
for k,f in [('immune_table_hash','immune_input.tsv.gz'),('window_table_hash','immune_windows.tsv.gz')]:assert sha(p/f)==lock[k]
features=json.loads((p/'natural_model.json').read_text())['features'];assert not any('mp10' in k.lower() or 'p95' in k.lower() or 'score' in k.lower() for k in features)
for f in ['discover.py','classify.py']:
    txt=(O/'scripts'/f).read_text()
    assert not any(k in txt for k in ["'spots/","'tables/section_effects","'tables/section_exposures","'tables/section_metrics_all"])
cs=read(p/'section_classification.tsv').set_index('section_id');eff=read(O/'tables/section_effects.tsv.gz');q=read(O/'audit/score_coverage_audit.tsv');bandtab=read(O/'tables/band_coverage.tsv.gz')
levels=read(O/'tables/all_section_levels.tsv');tests=read(O/'tables/main_six_tests.tsv');geom=read(O/'module_status.tsv')
assert len(q)==len(cs)==len(M)==63 and len(geom)==126 and len(tests)==6
assert q.E.sum()==q.eligible.sum()==133064 and q.n.sum()==265002 and q.max_error.max()<1e-7
assert bandtab.query("boundary=='all_E' and band in ['G1','G2'] and domain=='E_scoring'").n.sum()==0
total={'spots':0,'segments':0,'local_records':0}
for sec in M.index:
    tables=[par.ParquetFile(O/k/f'{sec}.parquet').metadata.num_rows for k in ['spots','boundary','local']]
    assert min(tables)>0
    for k,n in zip(total,tables):total[k]+=n
    d=load(f'spots/{sec}.parquet');assert (d.all_E & d.expression_QC & d.MP10.notna()).equals(d.new_epi_domain)
    guard=d.new_epi_domain;cut=d.MP10[guard].quantile(.95);assert np.array_equal((d.MP10[guard]>=cut).astype(int),d.P95_new[guard].astype(int))
    for b in B:
        r=eff.query("section_id==@sec and boundary==@b and width==200 and variant=='native'").iloc[0]
        dist=d['dist_'+b];near=d[guard&d[b]&(dist>=-200)&(dist<=0)];deep=d[guard&d[b]&(dist< -200)]
        assert (len(near),len(deep))==(r.n_B1,r.n_B2)
        assert np.isclose(r.score,near.MP10.mean()-deep.MP10.mean(),equal_nan=True)
    assert (O/'review'/f'{sec}.jpg').stat().st_size>4000
    assert levels[levels.section_id==sec].provisional_architecture_group.item()==cs.at[sec,'provisional_architecture_group']
f=read(O/'FIGURE_INDEX.tsv')
for name in f.figure:
    path=O/'figures'/(name+'.png')
    with Image.open(path) as pic:assert pic.size[0]>800 and min(pic.info['dpi'])>299
    assert (O/'figures'/(name+'.pdf')).read_bytes().startswith(b'%PDF')
for name in ['all_metric_associations.tsv.gz','major_lineage_contrasts.tsv.gz','paired_summary.tsv','level_composition_tests.tsv.gz','matched_architecture_comparisons.tsv','window_component_metrics_with_state.parquet']:
    assert (O/'tables'/name).stat().st_size>80
# Count full domain and second mask empty/equal cases independently in audit.py.
vals=dict(validated=True,n_sections=len(cs),regions_score_domain=int(q.eligible.sum()),raw_score_rechecked=int(q.n_check.sum()),max_abs_score_error=q.max_error.max(),formal_figures=len(f),reviews=len(cs),numeric_main_tests=int(tests.p.notna().sum()),planned_main_tests=len(tests),observed_empty_second_masks=int(((geom.boundary==B[1])&~geom.mask_exists).sum()),observed_equal_masks=int(geom.equal_masks.sum()/2),classification_lock_verified=True,classification_inputs='An immune-only matrix and windows; score column names excluded from the fixed feature list.',single_shared_classification_P95_exposure=True,all_E_outer_E_domain_count=0,**total)
(O/'audit/delivery_validation.json').write_text(json.dumps(vals,indent=2))
data=[]
for log in sorted((O/'logs').glob('time_*.log')):
    s=log.read_text();v=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',s)
    if v:data.append(dict(stage=log.stem[5:],rss_kib=int(v[1]),status=int(re.search(r'Exit status:\s*(\d+)',s)[1])))
assert all(a['status']==0 for a in data)
save(data,'audit/resources_by_stage.tsv')
rs=dict(max_recorded_stage_RSS_KiB=max(a['rss_kib'] for a in data),sum_of_stage_RSS_maxima_KiB_conservative_upper_bound=sum(a['rss_kib'] for a in data),cpu_affinity=list(os.sched_getaffinity(0)),native_threadpools=[(j['internal_api'],j['num_threads']) for j in threadpoolctl.threadpool_info()],lib_versions_file='software.tsv',per_process_address_space_limit_KiB=185000000,available_memory_at_verification_KiB=int(re.search(r'MemAvailable:\s*(\d+)',Path('/proc/meminfo').read_text())[1]),physical_ram_KiB=int(re.search(r'MemTotal:\s*(\d+)',Path('/proc/meminfo').read_text())[1]),resource_scope='One pinned compute core and one BLAS/OMP thread per invocation, sequential section processing. Per-stage maxima from GNU time; concurrent RSS was not continuously sampled. Conservative bound adds all recorded stage maxima. Excludes caller shell, desktop and unrelated users.')
(O/'audit/resources.json').write_text(json.dumps(rs,indent=2));print(json.dumps(vals))

