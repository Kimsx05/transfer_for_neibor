from run_v3_2 import *
import html
summ=[];blocks=[];small=[];evi=[];bd=[];labels=[];trans=[];panels=[];domains=[];identity=[]
cross=pd.read_csv(R/'spatial_ecology_top20_v2/patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv',sep='\t')
unit=pd.read_csv(OLD/'tables/v3_original_regions_unit_reaggregation.tsv',sep='\t',usecols=['dataset_id','analysis_unit_uid','analysis_unit_type_metadata','source_sections'])
for sec in M.index:
 dataset=M.loc[sec,'dataset_id'];p=O/'by_dataset'/dataset/sec
 if not (p/'complete.json').exists():continue
 info=json.loads((p/'complete.json').read_text());summ.append(info)
 for name,arr in [('tissue_component_summary.tsv',blocks),('small_roi_decisions.tsv',small),('evidence_distributions.tsv',evi),('boundary_evidence.tsv',bd),('domain_annotations.tsv',domains)]:arr.append(pd.read_csv(p/name,sep='\t'))
 z=pd.read_csv(p/'spot_labels.tsv.gz',sep='\t',keep_default_na=False);labels.append(z);a=z.groupby(['v3_1_lesion_location','lesion_location'],dropna=False).size().reset_index(name='n_spots');a['section_id']=sec;a['dataset_id']=dataset;trans.append(a)
 cp=cross[cross.section_id.eq(sec)];units=unit[(unit.dataset_id==dataset)&unit.source_sections.fillna('').apply(lambda t:sec in [x.strip() for x in t.split(';')])]
 uid=units.analysis_unit_uid.dropna().unique()
 if len(uid)>1:raise ValueError((sec,uid))
 status='; '.join(sorted({v.strip() for t in units.analysis_unit_type_metadata.dropna() for v in t.split(';')}))
 if len(cp):
  orig=cp.iloc[0].to_dict();orig['dataset_id']=dataset
 else:orig={'section_id':sec,'dataset_id':dataset}
 # Preserve v3.1 fallback, explicitly section-level and never a newly inferred patient.
 orig['analysis_unit_uid_v3_1']=uid[0] if len(uid) else (orig.get('analysis_unit_uid') if pd.notna(orig.get('analysis_unit_uid')) else '')
 orig['analysis_unit_type_v3_1']=status;orig['mapping_policy']='preserved_v3.1; no new patient or T/N pairing';identity.append(orig)
 # Nearby margin changes for each old small ROI, calculated from new bands, not inferred from counts.
 d=pd.read_parquet(p/'spot_regions_v3_2.parquet');adj,_=graph(d);sm=small[-1]
 for i,q in sm.iterrows():
  ids=np.flatnonzero(d.v3_1_lesion_roi_id.eq(q.old_roi_id));nb=sorted(set(ids)|{j for k in ids for j in adj[k]});oldb=d.v3_1_distance_band_200um.iloc[nb];newb=d.distance_band_200um.iloc[nb]
  sm.loc[i,'one_ring_neighborhood_n']=len(nb);sm.loc[i,'nearby_old_margin_n']=oldb.isin(['Inner margin','Outer margin']).sum();sm.loc[i,'nearby_new_margin_n']=newb.isin(['Inner margin','Outer margin']).sum();sm.loc[i,'nearby_band_label_changed_n']=oldb.ne(newb).sum()
 # Four core figure pairs must really exist before indexing a section as complete.
 figs=['01_location_comparison_HE','02_bands_distance_reliability','03_epithelial_molecular_evidence','04_domain_CNV_small_regions']
 if (p/'05_author_S1_comparison.png').exists():figs+=['05_author_S1_comparison']
 for f in figs:
  for ext in ['png','pdf']:assert (p/(f+'.'+ext)).is_file()
 title=html.escape(sec);page='<html><meta charset="utf-8"><title>'+title+'</title><style>body{font:16px sans-serif;max-width:1500px;margin:24px auto}img{width:100%}a{margin:6px}</style><h1>'+title+'</h1><p>v3.2 computational macro location. No new pathology ground truth. Gray/purple means extent unresolved, not non-lesion. Inside does not mean deep core.</p><p><a href="../../../index.html">Cohort index</a></p>'
 for f in figs:page+=f'<p><a href="{f}.png">PNG</a><a href="{f}.pdf">PDF</a></p><img src="{f}.png" loading="lazy">'
 page+='<p><a href="spot_labels.tsv.gz">Spot labels TSV.gz</a> <a href="spot_regions_v3_2.parquet">Full preserved molecular + region parquet</a></p></html>';(p/'index.html').write_text(page)
S=pd.DataFrame(summ);S['macro_outer_analysis_status']=np.where(S.supported_outer_segments.gt(0),'local_candidate_requires_extent_check','no_supported_macro_outer; do_not_interpret_legacy_margin_as_outer');S.to_csv(O/'tables/section_summary.tsv',sep='\t',index=False)
for name,arr in [('tissue_component_summary',blocks),('small_roi_decisions',small),('evidence_distributions',evi),('boundary_evidence',bd),('domain_annotations',domains),('label_transitions',trans)]:pd.concat(arr,ignore_index=True).to_csv(O/'tables'/f'{name}.tsv',sep='\t',index=False)
L=pd.concat(labels,ignore_index=True);assert not L[['section_id','barcode']].duplicated().any();L.to_csv(O/'tables/all_spot_labels_v3_2.tsv.gz',sep='\t',index=False)
# Re-read typed labels so missing geometry remains numeric NA, not the display empty string.
fields=list(L.columns);typed=[]
for sec in S.section_id:typed.append(pd.read_parquet(O/'by_dataset'/M.loc[sec,'dataset_id']/sec/'spot_regions_v3_2.parquet',columns=fields))
pd.concat(typed,ignore_index=True).to_parquet(O/'tables/all_spot_labels_v3_2.parquet',index=False)
pd.DataFrame(identity).to_csv(O/'tables/preserved_identity_crosswalk.tsv',sep='\t',index=False)
# Explicitly keep the v3.1 correction metadata, not its biological result values.
shutil.copyfile(OLD/'tables/unit_key_correction.tsv',O/'tables/v3_1_unit_key_correction_metadata.tsv');shutil.copyfile(OLD/'tables/missing_unit_identity.tsv',O/'tables/v3_1_missing_unit_identity.tsv')
B=pd.concat(blocks,ignore_index=True);B[B.status.str.startswith('all_')].to_csv(O/'tables/single_region_components.tsv',sep='\t',index=False)
index='<html><meta charset="utf-8"><title>BLCA pathology regions v3.2</title><style>body{font:15px sans-serif;margin:30px}table{border-collapse:collapse}td,th{padding:6px;border:1px solid #ddd}img{max-width:1400px;width:100%}</style><h1>BLCA v3.2 — macro location and interface revision</h1><p>Red: computational inside. Blue: positive non-lesion context. Purple: extent unresolved. No edge means no margin/core. No new MP10/TME analysis.</p><p><a href="METHODS.md">Methods</a> · <a href="RESULTS.md">Results</a> · <a href="tables/section_summary.tsv">Section summary</a> · <a href="tables/reference_comparison.tsv">Published annotation comparison</a></p>'
if not S.supported_outer_segments.gt(0).any():index+='<p style="padding:14px;background:#fff0d2"><b>本版没有获得足够支持的宏观病灶外缘；63张均不宜把旧margin直接解释为宏观肿瘤外缘。</b> 已保留内部候选、正向非病灶背景及具体范围未定区域。详见结果说明。</p>'
for dataset,z in S.groupby('dataset_id',sort=False):
 fig,axs=plt.subplots(int(np.ceil(len(z)/4)),4,figsize=(16,3.8*int(np.ceil(len(z)/4))),squeeze=False)
 for ax in axs.ravel():ax.axis('off')
 for ax,sec in zip(axs.ravel(),z.section_id):
  d=pd.read_parquet(O/'by_dataset'/dataset/sec/'spot_regions_v3_2.parquet',columns=['x','y','lesion_location']);st=S.set_index('section_id').loc[sec];ax.scatter(d.x,d.y,c=[LP[v] for v in d.lesion_location],s=2,linewidths=0);ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(sec+f'\nouter segments={st.supported_outer_segments}; unresolved={st.extent_unresolved}',fontsize=8)
 fig.suptitle(dataset+' | red inside; blue non-lesion-supported; purple extent unresolved');fig.tight_layout();path=O/'figures'/f'{dataset}_overview';fig.savefig(str(path)+'.png',dpi=300);fig.savefig(str(path)+'.pdf');plt.close(fig)
 index+=f'<h2>{html.escape(dataset)}</h2><a href="figures/{dataset}_overview.pdf">Overview PDF</a><img loading="lazy" src="figures/{dataset}_overview.png"><table><tr><th>Section</th><th>Spots</th><th>Inside</th><th>Outside</th><th>Unresolved</th><th>Supported segments</th><th>Exact distance</th><th>200um band determined</th></tr>'
 for q in z.itertuples():index+=f'<tr><td><a href="by_dataset/{dataset}/{q.section_id}/index.html">{q.section_id}</a></td><td>{q.n_spots}</td><td>{q.inside}</td><td>{q.outside}</td><td>{q.extent_unresolved}</td><td>{q.supported_outer_segments}</td><td>{q.exact_estimable}</td><td>{q.band_200_determined}</td></tr>'
 index+='</table>'
index+='</html>';(O/'index.html').write_text(index)
print(S.to_string(index=False));print('TOTAL',S.select_dtypes('number').sum().to_dict(),flush=True)
