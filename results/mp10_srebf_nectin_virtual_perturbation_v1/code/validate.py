from common import *
import datetime, sys,subprocess
from PIL import Image
import fitz
s=rd('tables/virtual_perturbation_summary.tsv');r=rd('tables/virtual_perturbation_robustness.tsv');c=rd('tables/virtual_perturbation_coverage.tsv')
need='analysis_id dataset_id method method_version state perturbed_gene endpoint_gene n_patients n_cells metric_type effect_scale effect_value direction_available predicted_direction uncertainty_type uncertainty_low uncertainty_high p_value q_value assessment_status ne_reason prior_id config_id'.split()
assert set(need)<=set(s) and len(r)==52 and set(r.analysis_type)=={'LODO','patient_subsample','alternative_seed'}
for t in [s,r]:
 assert not t.duplicated(['analysis_id','perturbed_gene','endpoint_gene']).any()
 assert set(t.assessment_status)<=set(['EVALUABLE','NE','EXPLORATORY'])
 ne=t[t.assessment_status=='NE'];assert ne.ne_reason.notna().all() and ne[['effect_value','p_value','q_value','uncertainty_low','uncertainty_high','direction_score']].isna().all().all()
 assert t[t.assessment_status!='NE'].effect_value.notna().all()
 assert t.loc[~t.direction_available,'predicted_direction'].isna().all()
 assert not (t.p_value.dropna()<0).any() and not (t.q_value.dropna()>1).any()
 assert (t.uncertainty_low.dropna()<=t.loc[t.uncertainty_low.notna(),'uncertainty_high']).all()
main=s[(s.analysis_id=='main')&(s.endpoint_gene=='NECTIN2')&s.perturbed_gene.isin(['SREBF1','SREBF2'])]
assert len(main)==2 and main.n_cells.eq(8000).all()
src=rd('qc/jobs/main_effects.tsv.gz');m=main.merge(src,left_on=['perturbed_gene','endpoint_gene'],right_on=['perturbed_gene','gene']);assert np.allclose(m.distance,m.effect_value,rtol=1e-12) and np.allclose(m['p.adj'],m.q_value)
jobs=rd('tables/resource_usage.tsv');assert len(jobs)==38 and (jobs.reason=='fitted').sum()==37
for j in jobs.analysis_id:assert 'Exit status: 0' in Path('logs/usage_'+j+'.txt').read_text()
meta=rd('cache/metadata.tsv.gz'); assert meta[meta.main].barcode.nunique()==8000
for j in r.analysis_id.unique():
 m=rd('qc/jobs/'+j+'_membership.tsv.gz');assert m.barcode.is_unique
 if j.startswith('sub'):assert m.patient_id.nunique()==31
 if j.startswith('drop_'):assert j[5:] not in set(m.dataset)
for k in ['SREBF1','SREBF2']:
 assert len(r[(r.perturbed_gene==k)&(r.analysis_type=='patient_subsample')])==20
 assert len(rd('tables/matched_controls.tsv').query('matched_to==@k'))==5
info=json.loads(Path('qc/input.json').read_text());assert info['ncount_mismatch']==info['nfeature_mismatch']==0
assert c.raw_feature[c.gene=='NECTIN2'].eq('NECTIN2').all()
fig=rd('FIGURE_INDEX.tsv');assert len(fig)>=9;tests=[]
for v in fig.itertuples():
 d=rd(v.plot_data);assert len(d)>0
 with Image.open(v.png) as img:
  dpi=img.info['dpi'];assert all(abs(x-300)<1 for x in dpi);width,height=img.size
 with fitz.open(v.pdf) as doc:
  assert len(doc)==1;pg=doc[0];txt=pg.get_text();assert len(txt.strip())>15
  fonts=sorted({f[3] for f in pg.get_fonts()});assert all('Arial' in f for f in fonts)
  assert abs(pg.rect.width/72*300-width)<4 and abs(pg.rect.height/72*300-height)<4
  assert not pg.get_images();assert len(pg.get_drawings())>0
  words=pg.get_text('words');assert all(w[0]>=-1 and w[1]>=-1 and w[2]<=pg.rect.width+1 and w[3]<=pg.rect.height+1 for w in words)
 tests.append(dict(figure=v.figure,n_direct_rows=len(d),width_px=width,height_px=height,dpi=dpi[0],fonts=';'.join(fonts),vector_pdf=True))
save(pd.DataFrame(tests),'qc/figure_checks.tsv')
rel=['tables/virtual_perturbation_summary.tsv','tables/virtual_perturbation_robustness.tsv','tables/virtual_perturbation_coverage.tsv']
for p in ['METHODS.md','RESULTS.md','run_all','replot','config/input_fingerprints.tsv','config/scTenifold.json','FIGURE_INDEX.tsv']:assert Path(p).stat().st_size>0
peak=rd('logs/simultaneous_memory.tsv').total_RSS_bytes.max()
usage=jobs.sort_values('peak_RSS_KiB',ascending=False);outsize=int(subprocess.check_output(['du','-sb','.']).split()[0])
pstart=Path('logs/fit_pilot.log').stat().st_mtime-jobs.set_index('analysis_id').loc['pilot','wall_s']
end=max(Path('logs/usage_'+j+'.txt').stat().st_mtime for j in jobs.analysis_id)
stats=dict(delivery_bytes_before_manifest=outsize,fit_window_wall_s=end-pstart,successful_job_cpu_s=jobs.cpu_s.sum(),max_single_process_RSS_KiB=int(usage.iloc[0].peak_RSS_KiB),sampled_simultaneous_RSS_bytes=int(peak),n_fitted_configurations=int((jobs.reason=='fitted').sum()),n_attempted_configurations=len(jobs))
Path('qc/delivery_metrics.json').write_text(json.dumps(stats,indent=2))
results=Path('RESULTS.md')
section='\n## 资源核验\n'
results.write_text(results.read_text().split(section)[0]+section+f'\n从小试开始到最后一次正式作业退出，实测墙钟 {(end-pstart)/60:.1f} 分钟。成功计算作业累计 CPU {jobs.cpu_s.sum()/60:.1f} 分钟。校验时交付目录逻辑大小 {outsize/2**20:.1f} MiB。单进程 RSS 与并行采样值的统计口径见 qc/delivery_metrics.json。\n')
files=[dict(path=p,sha256=sha(p),bytes=Path(p).stat().st_size,rows=len(rd(p))) for p in rel]
manifest=dict(status='complete',analysis_id='mp10_srebf_nectin_virtual_perturbation_v1',config_id='epi-v1-seed123/scTenifold-v2',method='scTenifoldKnk',method_version='2.0.0',root=str(D),access='read-only consumer interface',files=files,config_fingerprints={p:sha(p) for p in ['config/design.json','config/scTenifold.json','config/input_fingerprints.tsv','config/prior_audit.json']},finished_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),figures=int(len(fig)))
p=Path('config/virtual_perturbation_manifest.json');tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n');tmp.chmod(0o444)
for entry in files:Path(entry["path"]).chmod(0o444)
os.replace(tmp,p)
print(json.dumps(stats))

