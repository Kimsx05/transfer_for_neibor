from common3 import *
import time,re
from PIL import Image
from threadpoolctl import threadpool_info
import pyarrow.parquet as pq
h=read(O/'tables/section_bins.tsv.gz');pr=read(O/'tables/section_predictions.tsv.gz')
files=read(O/'audit/section_figure_index.tsv');status=read(O/'audit/section_boundary_status.tsv')
cache=read(O/'cache/point_index.tsv').set_index('section_id')
assert M.tumor.sum()==60 and len(M)==63 and len(files)==504 and len(status)==126
assert len(files.dropna(subset='file'))==488
for s in M.index:
    d=spot(s);v=d.expression_QC&np.isfinite(d.MP10);q=v&d.all_E;mu=d.MP10[q].mean()
    assert d.new_epi_domain.equals(q) and np.isclose(mu,cache.at[s,'center']) and cache.at[s,'zero']==(v&d.MP10.eq(0)).sum()
    for b in B:
        ok=v&d['dist_'+b].notna();assert not ((d['dist_'+b]<0)&~d[b]|(d['dist_'+b]>0)&d[b]).any()
        assert not d.loc[ok,'region_main'].eq('Unknown').any()
        side=d['dist_'+b]>0
        assert not (q&side).any() if b==B[0] else not (q&side&~d.region_main.isin(['EIS','EIN'])).any()
        bins=h[(h.section_id==s)&(h.boundary==b)]
        assert bins.n.sum()==ok.sum() and bins[bins.domain=='nonE_outer'].p95_fraction.isna().all()
        assert np.allclose(bins.centered,bins['median']-mu,atol=1e-10)
        a=d[ok].copy();di=a['dist_'+b];a['dom']=np.where(a[b],'E_inner',np.where(a.all_E,'EIN_EIS_outer','nonE_outer'))
        a['bin']=np.where(a[b],-1,1)*(np.maximum(0,np.ceil(abs(di)/100-1e-10)-1)*100+50)
        if len(a):
            agg=a.groupby(['dom','bin']).MP10.agg(['size','median','mean']).sort_index();saved=bins.set_index(['domain','center']).sort_index();assert np.allclose(agg.to_numpy(),saved[['n','median','mean']].to_numpy())
        assert len(files[(files.section_id==s)&(files.boundary==b)].dropna(subset='file'))==(4 if ok.any() else 0)
for (s,b,dom,seg),g in pr.groupby(['section_id','boundary','domain','segment']):
    z=h[(h.section_id==s)&(h.boundary==b)&(h.domain==dom)&h.center.between(g.x.min(),g.x.max())]
    assert z.n.ge(5).all() and (np.diff(sorted(z.center))==100).all()
    if g.method.iloc[0]=='penalized_cubic':assert len(z)>=5
    assert (g.x<0).all() if dom=='E_inner' else (g.x>0).all()
for path in files.file.dropna():
    q=O/path;assert q.stat().st_size>1000 and q.read_bytes().startswith(b'%PDF')
cl=read(O/'classification/section_classification.tsv');lock=json.loads((O/'classification_lock.json').read_text())
for f,sign in lock['hashes'].items():assert sha(Path(f))==sign
assert lock['labels']==cl.set_index('section_id').exploratory_cluster.to_dict()
for f in ['cache/immune_metadata.tsv','cache/immune_features.tsv','cache/immune_windows.parquet','cache/immune_links.tsv.gz','classification/X_z.tsv','classification/X_pam.tsv']:
    cols=pq.read_schema(O/f).names if f.endswith('parquet') else pd.read_csv(O/f,sep='\t',nrows=0).columns
    assert not any(any(k in c.lower() for k in ['mp10','p95','highlow']) for c in cols)
assert (O/'classification_lock.json').stat().st_mtime <= (O/'tables/four_primary_tests.tsv').stat().st_mtime
main=read(O/'tables/four_primary_tests.tsv');assert len(main)==4 and np.allclose(main.q,bh(main.p,4),equal_nan=True)
assert len(read(O/'review/index.tsv'))==63 and set(files.section_id)==set(M.index)
boot=read(O/'classification/bootstrap_runs.tsv.gz');assert boot.groupby(['rule','method','K']).size().eq(100).all()
gf=read(O/'audit/all_group_figure_index.tsv')
for r in gf.itertuples():
    assert (O/r.pdf).read_bytes().startswith(b'%PDF')
    with Image.open(O/r.png) as im:assert min(im.info['dpi'])>299
def jd(p,data):p.write_text(json.dumps(data,indent=2,ensure_ascii=False))
jd(O/'audit/acceptance.json',dict(sections=len(M),tumor=int(M.tumor.sum()),normals=int((~M.tumor).sum()),source_points=int(cache.n.sum()),scored_E=int(cache.D_MP10.sum()),real_zero_points_tissue=int(cache.zero.sum()),raw_bins=len(h),valid_section_PDFs=int(files.file.notna().sum()),group_figures=len(gf),review_maps=len(cl),status_rows=len(status),no_interface=int(status.display_n.eq(0).sum()),identical_to_v2_score_and_sign_in_bins=True,no_cross_segment_or_side=True,shared_classifications=True,primary_tests=4,bootstrap_success=main.n_boot.tolist(),review_commit_status=read(O/'audit/actual_inputs_vs_commit.tsv').status.value_counts().to_dict(),notes='Full scatter reads exact v2 scores; every bin mean, median and denominator recomputed. Separate density/entry labels from inferred malignancy.'))
rss=[];elapsed=[]
def seconds(v):return sum(float(t)*60**j for j,t in enumerate(v.split(':')[::-1]))
for log in sorted((O/'logs').glob('time_*.log')):
    st=log.read_text();mem=re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)',st);tm=re.search(r'Elapsed .*?: ([\d:.]+)',st);ec=re.search(r'Exit status: (\d+)',st)
    if mem and tm and ec:rss.append(dict(stage=log.stem[5:],rss_KiB=int(mem[1]),s=seconds(tm[1]),exit=int(ec[1])))
save(rss,O/'audit/stage_resources.tsv')
assert all(r['exit']==0 for r in rss)
paths=sorted(p for p in O.rglob('*') if p.is_file() and not any(x in p.parts for x in ['__pycache__','versions']) and 'light_replay_' not in str(p))
jd(O/'audit/resources.json',dict(largest_timed_process_RSS_KiB=max(v['rss_KiB'] for v in rss),sum_recorded_stage_maxima_KiB=sum(v['rss_KiB'] for v in rss),stage_seconds_sum_including_replay_and_repeated_runs=sum(v['s'] for v in rss),wall_span_seconds_at_validation=int(time.time()-min(p.stat().st_ctime for p in paths)),disk_bytes_at_validation=sum(p.stat().st_size for p in paths),per_stage_affinity=list(os.sched_getaffinity(0)),native_threads=[(j['internal_api'],j['num_threads']) for j in threadpool_info()],peak_scope='GNU time per stage, including lightweight replay; aggregate concurrent RSS not continuously sampled. Run entry pins compute to one CPU, sequential sections and BLAS/OMP=1; all stages use same assigned CPU.'))
paths=[p for p in O.rglob('*') if p.is_file() and not any(x in p.parts for x in ['versions','__pycache__','logs']) and p.name!='file_manifest.tsv']
save([dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)],O/'file_manifest.tsv')
print((O/'audit/acceptance.json').read_text())

