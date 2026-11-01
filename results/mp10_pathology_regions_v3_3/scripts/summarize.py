from run_v3_3 import *
import html,importlib.metadata as md
S=[];blocks=[];small=[];changes=[];labels=[];ev=[];refs=[]
figs=['01_location_versions','02_location_HE_support','03_bands_distance_reliability','04_epithelial_molecular_evidence','05_reference_CNV_support']
for sec in M.index:
 dataset=M.loc[sec,'dataset_id'];p=O/'by_dataset'/dataset/sec
 if not (p/'complete.json').exists():raise RuntimeError('unfinished '+sec)
 info=json.loads((p/'complete.json').read_text());S.append(info);d=pd.read_parquet(p/'spot_regions_v3_3.parquet');labelcols=list(pd.read_csv(p/'spot_labels.tsv.gz',sep='\t',nrows=0).columns);labels.append(d[labelcols]);blocks.append(pd.read_csv(p/'tissue_component_summary.tsv',sep='\t'));ev.append(pd.read_csv(p/'evidence_distributions.tsv',sep='\t'))
 q=d.groupby(['v3_1_lesion_location','v3_2_lesion_location','lesion_location','change_class'],dropna=False).size().reset_index(name='n_spots');q['section_id']=sec;q['dataset_id']=dataset;changes.append(q)
 # Complete the small-ROI final disposition after author constraints, including pilot files.
 
 try:sm=pd.read_csv(p/'small_roi_decisions.tsv',sep='\t')
 except pd.errors.EmptyDataError:sm=pd.DataFrame(columns=['section_id','old_roi_id','n_spots','v3_2_weak_roi','decision','final_locations','final_changed_from_v3_1'])
 adj,_=graph(d)
 for i,r in sm.iterrows():
  ix=np.flatnonzero(d.v3_1_lesion_roi_id.eq(r.old_roi_id));near=sorted(set(ix)|{j for k in ix for j in adj[k]});sm.loc[i,'final_locations']=';'.join(sorted(set(d.lesion_location.iloc[ix])));sm.loc[i,'final_changed_from_v3_1']=int(d.lesion_location.iloc[ix].ne(d.v3_1_lesion_location.iloc[ix]).sum());sm.loc[i,'author_coarse_constrained_spots']=int(d.author_coarse_reference.iloc[ix].isin(['author_coarse_inside','author_coarse_outside']).sum());sm.loc[i,'old_margin_neighborhood_n']=int(d.v3_1_distance_band_200um.iloc[near].isin(['Inner margin','Outer margin']).sum());sm.loc[i,'new_margin_neighborhood_n']=int(d.distance_band_200um.iloc[near].isin(['Inner margin','Outer margin']).sum())
 sm.to_csv(p/'small_roi_decisions.tsv',sep='\t',index=False);small.append(sm)
 ref={'section_id':sec,'dataset_id':dataset,'form':'no_usable_same_section_coarse_reference','source':'existing project sources; no new patient/source mapping','registration_status':'unavailable','coarse_inside_n':0,'coarse_outside_n':0,'inside_to_outside':info['author_inside_to_outside'],'outside_to_inside':info['author_outside_to_inside'],'boundary_support':'candidate only; no pathology confirmation'}
 rp=O/'registration'/sec/'registration.json'
 if rp.exists():
  meta=json.loads(rp.read_text());ref.update(form='published S1 coarse colour regions registered by H&E; not barcode truth',source='https://doi.org/10.1158/2159-8290.CD-26-0099',registration_status=meta['status'],coarse_inside_n=meta['n_coarse_inside'],coarse_outside_n=meta['n_coarse_outside'],boundary_support='author regional calibration; buffer/others not transferred; not independent accuracy')
 elif sec.startswith('GSE285715'):
  ref.update(source='https://doi.org/10.1073/pnas.2505125123',boundary_support='marker/composition panels; Tumor-419a to T3 unverified; candidate interface may be intralesional')
 elif sec.startswith('NPJ'):
  ref.update(source='https://doi.org/10.1038/s41698-026-01417-x',boundary_support='local alias bridge unverified; shared SpaCET evidence; no transferred author invasive front')
 refs.append(ref)
 # Index links only to figures that actually exist.
 sf=figs+(['06_author_registration','07_author_comparison'] if rp.exists() else [])
 for f in sf:
  for ext in ['png','pdf']:assert (p/(f+'.'+ext)).exists()
 page='<html><meta charset="utf-8"><title>'+sec+'</title><style>body{font:16px sans-serif;max-width:1500px;margin:24px auto}img{width:100%}</style><h1>'+sec+'</h1><p>v3.3 候选位置。红蓝具有相同计算估计地位；距离以当前候选ROI为条件，不等于病理前沿已确认。</p><p><a href="../../../index.html">全队列</a></p>'
 for f in sf:page+=f'<p><a href="{f}.pdf">PDF</a> · <a href="{f}.png">PNG</a></p><img src="{f}.png" loading="lazy">'
 page+='<p><a href="spot_labels.tsv.gz">标签TSV.gz</a> · <a href="spot_regions_v3_3.parquet">完整点表Parquet</a></p></html>';(p/'index.html').write_text(page)
S=pd.DataFrame(S);S['boundary_use']=np.where(S.candidate_edges.gt(0),'candidate_geometry_available_not_histological_front','no_observed_candidate_interface');S.to_csv(O/'tables/section_summary.tsv',sep='\t',index=False)
for name,arr in [('tissue_component_summary',blocks),('small_roi_decisions',small),('position_changes',changes),('evidence_distributions',ev)]:pd.concat(arr,ignore_index=True).to_csv(O/'tables'/f'{name}.tsv',sep='\t',index=False)
B=pd.concat(blocks);B[B.status.str.startswith('all_')].to_csv(O/'tables/single_region_components.tsv',sep='\t',index=False)
L=pd.concat(labels,ignore_index=True);assert not L[['section_id','barcode']].duplicated().any();L.to_parquet(O/'tables/all_spot_labels_v3_3.parquet',index=False);L.to_csv(O/'tables/all_spot_labels_v3_3.tsv.gz',sep='\t',index=False);pd.DataFrame(refs).to_csv(O/'tables/reference_comparison.tsv',sep='\t',index=False)
# Readable dataset-level availability, without a fake accuracy score.
index='<html><meta charset="utf-8"><title>BLCA v3.3</title><style>body{font:15px sans-serif;margin:30px}table{border-collapse:collapse}td,th{padding:6px;border:1px solid #ddd}img{width:100%;max-width:1400px}</style><h1>病灶候选区域 v3.3</h1><p>位置、支持程度、几何可计算性分开。红/蓝都是候选位置；可计算的margin不等于已确认宏观病理外缘。</p><p><a href="RESULTS.md">结果</a> · <a href="METHODS.md">方法/复现</a> · <a href="tables/section_summary.tsv">切片汇总</a> · <a href="tables/all_spot_labels_v3_3.parquet">全部标签</a> · <a href="tables/reference_comparison.tsv">作者对照</a></p>'
for dataset,z in S.groupby('dataset_id',sort=False):
 fig,axs=plt.subplots(int(np.ceil(len(z)/4)),4,figsize=(16,3.8*int(np.ceil(len(z)/4))),squeeze=False)
 for ax in axs.ravel():ax.axis('off')
 for ax,r in zip(axs.ravel(),z.itertuples()):
  d=pd.read_parquet(O/'by_dataset'/dataset/r.section_id/'spot_regions_v3_3.parquet',columns=['x','y','lesion_location']);ax.scatter(d.x,d.y,s=2,c=[LP[x] for x in d.lesion_location],linewidths=0);ax.invert_yaxis();ax.set_aspect('equal');ax.set_title(r.section_id+f'\ncandidate edges={r.candidate_edges}; exact n={r.exact_estimable}',fontsize=8)
 fig.suptitle(dataset+' | red/blue candidate positions, not confirmed malignancy/normality');fig.tight_layout();fp=O/'figures'/f'{dataset}_overview';fig.savefig(str(fp)+'.png',dpi=300);fig.savefig(str(fp)+'.pdf');plt.close(fig)
 index+=f'<h2>{dataset}</h2><a href="figures/{dataset}_overview.pdf">总览PDF</a><img loading="lazy" src="figures/{dataset}_overview.png"><table><tr><th>Section</th><th>Inside</th><th>Outside</th><th>Unlocalized</th><th>Candidate edges</th><th>Exact distances</th><th>200um bands</th><th>Author position changes</th></tr>'
 for r in z.itertuples():index+=f'<tr><td><a href="by_dataset/{dataset}/{r.section_id}/index.html">{r.section_id}</a></td><td>{r.inside}</td><td>{r.outside}</td><td>{r.unlocalized}</td><td>{r.candidate_edges}</td><td>{r.exact_estimable}</td><td>{r.band_200_determined}</td><td>{r.author_changed}</td></tr>'
 index+='</table>'
(O/'index.html').write_text(index+'</html>')
versions={p:md.version(p) for p in ['numpy','pandas','scipy','shapely','matplotlib','Pillow','pyarrow']};import cv2;versions['opencv']=cv2.__version__;versions['python']=sys.version;versions['seed']=20261030;(O/'software_versions.json').write_text(json.dumps(versions,indent=2))
print(S.to_string(index=False));print('TOTAL',S.select_dtypes('number').sum().to_dict(),flush=True)
