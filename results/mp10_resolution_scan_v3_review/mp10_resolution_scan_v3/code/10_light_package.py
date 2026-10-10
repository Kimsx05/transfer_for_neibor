from pathlib import Path
import csv,zipfile,hashlib,json
R=Path('/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_resolution_scan_v3')
with (R/'08_review/RECOMMENDED_CANDIDATES.tsv').open() as f:recs={x['candidate_id'] for x in csv.DictReader(f,delimiter='\t')}
with (R/'04_pseudobulk_markers/model_summary.tsv').open() as f:models=list(csv.DictReader(f,delimiter='\t'))
keys={x['model_key'] for x in models if x['candidate_id'] in recs}
readme='''轻量审阅包说明

打开REPORT.html。包含全部候选/排名/短名单/功能摘要、报告引用的主要ALL图、所有dataset的marker和来源覆盖图，以及两名全局F1冠军的完整DE/富集表。
未纳入其他候选的完整逐基因/逐通路表和大型缓存；这些完整表在服务器输出目录及mp10_resolution_scan_v3_review.zip中。NOT_TESTABLE仍有状态与资格记录。
所有69套原图、逐图数据及105个精确模型键的完整结果保存在完整输出目录；未新增分析、删细胞或重新选候选。
'''
(R/'LIGHT_PACKAGE_README.txt').write_text(readme)
files=[]
for p in sorted(R.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(R).as_posix();parts=p.relative_to(R).parts
 if p.suffix in ('.qs','.zip','.pyc') or '__pycache__' in rel:continue
 if rel in ['light_package_manifest.tsv','review_package_manifest.tsv','logs/09_package.log','logs/10_light_package.log']:continue
 include=False
 if len(parts)==1:include=True
 elif parts[0] in ['code','08_review','02_coverage','03_candidates','06_stability']:include=True
 elif parts[0]=='00_audit':include=p.name not in ['RNA_features.tsv.gz','gene_annotation.tsv.gz']
 elif parts[0]=='01_partitions':include=p.name.endswith('_sizes.tsv') or p.name=='all_cluster_sizes.tsv'
 elif parts[0]=='04_pseudobulk_markers':
  include=(len(parts)==2 and p.suffix=='.tsv') or (len(parts)>=4 and parts[2] in keys and p.name in ['DE_full.tsv.gz','shortlist.tsv','eligibility.tsv','status.tsv','model_samples.tsv','design.tsv','rank_mapping.tsv.gz'])
 elif parts[0]=='05_function':include=len(parts)==2 or (len(parts)>=3 and parts[1] in keys)
 elif parts[0]=='07_figures':
  n=p.name
  include=n in ['figure_index.tsv','color_limits.tsv','sample_display_key.tsv'] or n.startswith('ALL__') or any(t in n for t in ['marker_DotPlot','sample_precision_recall','candidate_source_composition','candidate_PR'])
  if n.endswith('_plotdata.tsv.gz') and any(t in n for t in ['candidate_highlights','shared_continuous']):include=False
 elif parts[0]=='logs':include=p.name.startswith('sessionInfo')
 if include:files.append(p)
manifest=[dict(relative_path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files]
with (R/'light_package_manifest.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=manifest[0],delimiter='\t');w.writeheader();w.writerows(manifest)
files.append(R/'light_package_manifest.tsv');zpath=R/'mp10_resolution_scan_v3_review_light.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,arcname='mp10_resolution_scan_v3/'+p.relative_to(R).as_posix())
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
print('Light ZIP',zpath,'bytes',zpath.stat().st_size,'files',len(files))
