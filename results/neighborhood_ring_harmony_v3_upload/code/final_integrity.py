from base import *
import zipfile,shutil
C=OUT/'04_comparison';R=OUT/'reports';s=pd.read_csv(C/'eighteen_scheme_comparison.tsv',sep='\t');assert len(s)==18 and s.config_id.nunique()==18
s['recommendation_status']=np.where((s.stability_ARI_median>=.8)&(s.silhouette_embedding>0)&(s.silhouette_raw_Hellinger>0)&s.convergence_status.ne('iteration_limit'),'niche_working_candidate','exploratory_with_limitations');tab(C/'eighteen_scheme_comparison.tsv',s)
text=(R/'RUN_SUMMARY.md').read_text()
for r in s.itertuples():
 p=OUT/'03_main'/r.config_id/'working_config.json';d=json.loads(p.read_text());d['recommendation_status']=r.recommendation_status;js(p,d)
 lines=text.splitlines()
 for i,line in enumerate(lines):
  if line.startswith('|'+r.config_id+'|'):lines[i]=f'|{r.config_id}|{r.n_clusters}|{r.stability_ARI_median:.3f}|{r.accepted_fraction:.1%}|{r.recommendation_status}|'
 text='\n'.join(lines)+'\n'
text+='\nSAMPLE R1 Leiden 的工作候选在原始 Hellinger 空间仍为负轮廓系数，属于探索性结果；高条件稳定性并不能消除这一限制。DATASET R1 最终取8群，SAMPLE R2取5群，均是既定有限网格中的组成分离受支持候选。\n'
(R/'RUN_SUMMARY.md').write_text(text)
methods=(R/'METHODS.md').read_text();methods+='\nFinite-repair deviation: SAMPLE_WEAK_R2 used an initial extra same-K K-means9 diagnostic, then a final extra K-means5 diagnostic after the composition-based Leiden selection revision. This exceeds the planned one extra same-K fit by one on this branch only; no extra representation fit or grid search was introduced. Final sameK_diagnostics.tsv refers to K5.\n';(R/'METHODS.md').write_text(methods)
# Ensure all requested fixed-section maps, ROI evidence and label checks exist before declaring complete.
b=pd.read_csv(OUT/'00_preflight/batch_key_mapping.tsv',sep='\t');assert len(b)==63
for row in b.itertuples():
 for ring in ['R1','R2','R12']:assert (C/'by_dataset'/row.dataset_id/(row.section_id+'__'+ring+'__six_modes.png')).exists()
roi=pd.read_csv(OUT/'05_user_review/ROI_review_to_fill.tsv',sep='\t');assert len(roi)==133 and set(['reviewed_label','review_comment','reviewer','review_date']).issubset(roi.columns);assert all(Path(p).exists() for p in roi.evidence_figure)
images=pd.read_csv(OUT/'05_user_review/image_availability.tsv',sep='\t');assert len(images)==63 and images.availability.eq('embedded image').all() and all(Path(p).exists() for p in images.HE)
assert (pd.read_csv(C/'validation.tsv',sep='\t').status=='PASS').all();scope=json.loads((OUT/'00_preflight/RUN_SCOPE.json').read_text());assert scope['stop_after']=='niche' and not scope['mp10_values_loaded'] and not scope['mp10_analysis_executed']
manifest=json.loads((OUT/'frozen_niche_manifest.json').read_text());manifest['configs']=s.to_dict('records');manifest['deviations'].append('Finite selection repair: SAMPLE_WEAK_R2 initially used Leiden K9 and one extra same-K K-means9; final composition-supported Leiden K5 required K-means5. Thus this branch used two extra same-K fits instead of the planned one; no further grid expansion or representation refit.');manifest['niche_labels_frozen']=True;manifest['n_target_sections']=63;manifest['n_section_ring_spatial_panels']=189;manifest['n_priority_ROI']=133;manifest['n_registered_embedded_HE_sections']=63;manifest['final_integrity']='PASS';manifest['working_config_hashes']={r.config_id:sha(OUT/'03_main'/r.config_id/'working_config.json') for r in s.itertuples()};manifest['code_hashes']={str(p):sha(p) for p in (OUT/'code').glob('*') if p.is_file()};manifest['task_list']=str(R/'TASKLIST.md');js(OUT/'frozen_niche_manifest.json',manifest)
js(R/'FINAL_INTEGRITY.json',dict(status='PASS',main_configurations=18,label_validation='18 PASS',section_ring_maps=189,priority_ROI=133,HE_sections=63,scope='niche only',MP10_analysis=False,source_validation='DEFERRED_BY_USER'))
# Refresh lightweight package without duplicate entries, including the final state and minimal restart code.
updates=[C/'eighteen_scheme_comparison.tsv',R/'RUN_SUMMARY.md',R/'METHODS.md',OUT/'frozen_niche_manifest.json',R/'TASKLIST.md',R/'FINAL_INTEGRITY.json']+list((OUT/'code').glob('*.py'))+list((OUT/'code').glob('*.R'))+list((OUT/'code').glob('*.sh'))+list((OUT/'code').glob('*.json'))
u={str(p.relative_to(OUT)):p for p in updates};temp=R/'review.finalizing.zip'
with zipfile.ZipFile(R/'review.zip') as old,zipfile.ZipFile(temp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as new:
 for name in old.namelist():
  if name not in u:new.writestr(name,old.read(name))
 for name,p in u.items():new.write(p,name)
temp.replace(R/'review.zip')
js(R/'package_manifest.json',dict(review_zip_bytes=(R/'review.zip').stat().st_size,review_zip_sha256=sha(R/'review.zip'),dataset_packs=[p.name for p in R.glob('*__spatial_figures.zip')]))
log('FINAL INTEGRITY PASS')
