from common import *
import hashlib,shutil,zipfile,argparse
from load_result import load_enrichment,attach_to_regions
p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
# Hash only small task code/config, not expression files or the server filesystem.
files=sorted(list((R/'code').glob('*.py'))+list((R/'code').glob('*.R'))+list(R.glob('*.sh'))+[R/'frozen_config.yaml'])
code={str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files}
audit=pd.read_csv(R/'definitions_audit.tsv',sep='\t');sources=[Path(x)for x in audit.path.dropna()]+[V/'frozen_config.yaml',F/'tables/expanded_mixtures.parquet']
sign={'code_config':code,'source_metadata':{str(f):{'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns}for f in sources}}
lock=R/'execution_identity.json'
if args.check:
 if lock.exists():
  assert sign==json.loads(lock.read_text()),'Task code/config or source identity changed. Review change and create explicit version; do not silently reuse cached results.'
 print('EXECUTION_IDENTITY_OK');sys.exit(0)
verify=[]
for s in SECTIONS:
 z=load_enrichment(s);t=z['T'];e=z['Epi'];n=z['neighbors'];a=attach_to_regions(t[['region_id']].copy(),s);assert np.array_equal(a.v3_FOXP3_supported,t.FOXP3_supported);assert len(z['effects'])==8
 assert len(n)==2*t.T_eligible.sum();assert not n.duplicated(['region_id','radius_um']).any()
 verify.append(dict(section=s,T_regions=len(t),T_eligible=int(t.T_eligible.sum()),Epi_eligible=len(e),neighbor_rows=len(n),effect_rows=8,loader_exact=True))
put(pd.DataFrame(verify),'loader_reread.tsv')
ix=pd.read_csv(R/'tables/figure_index.tsv',sep='\t');assert ix.figure.is_unique
for name in ix.figure:
 assert (R/'figures'/(name+'.png')).is_file()and(R/'figures'/(name+'.pdf')).is_file()
roi=pd.read_csv(R/'tables/fixed_ROI_v3_evidence.tsv',sep='\t');assert len(roi)==24 and not roi.duplicated(['section','roi']).any()
# All retained tabular deliveries readable. Parquet has already been loaded with IDs above.
for f in (R/'tables').glob('*.tsv'):pd.read_csv(f,sep='\t')
assert 'R_LOADER_OK' in (R/'logs/R_loader_validation.log').read_text()
# Check dimensions of actual PNGs and draw representative review statistics from reread tables.
from PIL import Image
for name in ['FE_H05_r200','Cervilla_ROI_candidate_1','D1_Johnson_program','distribution_Treg_main']:
 with Image.open(R/'figures'/(name+'.png'))as im:assert min(im.size)>900 and abs(im.info.get('dpi',(300,300))[0]-300)<1
# Native skill metadata included; exact source scoring is reused, no changed scoring package.
refs=pd.DataFrame([{'item':'UCell implementation','path':str(V/'code/04_native_UCell.R'),'note':'historical native2.14 check / average ties / maxRank1300; cached values reused'}, {'item':'historical rank scoring and panels','path':str(V/'code/common.py'),'note':'MP10 full18085 average ties K905'}, {'item':'frozen definitions','path':str(V/'frozen_config.yaml'),'note':'unchanged'}, {'item':'source pair semantics','path':str(V/'code/05_neighbors.py'),'note':'same-component, whole shared-bin pair removal'}, {'item':'plot skill','path':'/home/data/t070721/.agents/skills/jins-routine-plotting/SKILL.md','note':'Violin_distribution native template; custom spatial/forest layout uses same style'}]);put(refs,'implementation_sources.tsv')
# No big matrices or reference caches in review. Original raw atlas subset fixed by old ROI identity.
review=R/'review';review.mkdir(exist_ok=True);(review/'tables').mkdir(exist_ok=True);(review/'figures').mkdir(exist_ok=True);(review/'code').mkdir(exist_ok=True)
for f in list(R.glob('*.md'))+list(R.glob('*.yaml'))+list(R.glob('*.tsv'))+list(R.glob('*.sh')):shutil.copy2(f,review/f.name)
for f in (R/'code').iterdir():
 if f.suffix in ['.py','.R']:shutil.copy2(f,review/'code'/f.name)
for f in (R/'tables').glob('*.tsv'):
 if f.name not in ['null_histograms.tsv','continuous_and_ROI_robustness.tsv','spatial_robustness.tsv']:shutil.copy2(f,review/'tables'/f.name)
# Compact full per-eligible-region data allow independent label/denominator review, no expression matrices.
treview=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet');t=t[t.T_eligible];treview.append(t)
pd.concat(treview,ignore_index=True).to_parquet(review/'tables/T_eligible_regions.parquet',index=False)
selected=['definition_counts','label_overlap','distribution_Treg_main','distribution_log10_UMI','distribution_log1p_FOXP3','FE_H05_r50','FE_H05_r200','edge_Z_H05_r200','continuous_effect_r200','FE_block_sensitivity','hotspot_influence','calibration_full_rule','NMS_label_to_high_Epi_r200','NMS_high_Epi_to_label_r200']
for s in SECTIONS:selected +=[s+'_labels',s+'_evaluability',s+'_Epi_H05',s+'_program',s+'_distance_FOXP3_supported',s+'_ROI_candidate_1',s+'_ROI_control_1']
for name in selected:shutil.copy2(R/'figures'/(name+'.png'),review/'figures'/(name+'.png'))
for s,roi_id in [('Cervilla','candidate_1'),('D1_Johnson','candidate_1'),('UC_SCC','candidate_3'),('MPBC1','candidate_1'),('MPBC2','candidate_2')]:
 f=F/'figures'/f'{s}_{roi_id}_512atlas.png';shutil.copy2(f,review/'figures'/('v2_raw_'+f.name))
(review/'REVIEW_SCOPE.md').write_text('本包包含主要摘要、全部可评价T区域轻量信息、预定代表PNG和原固定ROI原始证据图。全101张图及PDF、各区成员和邻域明细在服务器目录；不是表达矩阵或全部原始数据。v2_raw_为原图原样复制用于对照，不是本轮重新分析。代码中的原依赖路径见manifest。\n')
# Avoid keeping a redundant duplicate of the distribution table in two formats.
f=R/'tables/distribution_plot_data.parquet'
if f.exists():f.unlink() # generated by our plotting step, gz TSV is the sole persistent plotting copy
lock.write_text(json.dumps(sign,indent=2))
# Lightweight review package; resource sizes measured after each final archive write.
for _ in range(2):
 total=sum(f.stat().st_size for f in R.rglob('*')if f.is_file()and f.name!='review.zip');review_bytes=sum(f.stat().st_size for f in review.rglob('*')if f.is_file());zip_bytes=(R/'review.zip').stat().st_size if (R/'review.zip').exists()else 0
 resources=pd.DataFrame([dict(metric='result_bytes_excluding_zip_including_review',value=total),dict(metric='review_directory_bytes',value=review_bytes),dict(metric='review_zip_bytes',value=zip_bytes),dict(metric='total_bytes_including_zip',value=total+zip_bytes),dict(metric='maximum_recorded_analysis_RSS_KiB',value=pd.read_csv(R/'tables/runtime_resources.tsv',sep='\t').peak_RSS_KiB.max()),dict(metric='total_threads_limit',value=4),dict(metric='figures',value=len(ix)),dict(metric='fixed_ROIs',value=24)])
 put(resources,'resource_and_delivery.tsv');shutil.copy2(R/'tables/resource_and_delivery.tsv',review/'tables/resource_and_delivery.tsv')
 with zipfile.ZipFile(R/'review.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6)as z:
  for f in sorted(review.rglob('*')):
   if f.is_file():z.write(f,str(f.relative_to(review)))
with zipfile.ZipFile(R/'review.zip')as z:assert z.testzip()is None;assert not any(x.endswith(('.h5ad','.rds','.h5','.npz'))for x in z.namelist())
print('FINAL_READ_AND_PACKAGE_OK',len(ix),'figures',len(roi),'ROIs',R/'review.zip')
