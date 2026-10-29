from revise import *
S=[];D=[];small=[];gaps=[];checks=[]
for sec,r in M.iterrows():
 p=OUT/'by_dataset'/r.dataset_id/sec;s=json.loads((p/'complete.json').read_text());S.append(s);d=pd.read_parquet(p/'spot_regions_v3_1.parquet');old=pd.read_parquet(OLD/'by_dataset'/r.dataset_id/sec/'spot_regions_v3.parquet');assert len(d)==len(old) and d.barcode.equals(old.barcode)
 for c in ['lesion_location','lesion_roi_id','malignant_fraction','MP10_original','CNV_support','spatial_domain']:assert d[c].equals(old[c]),(sec,c)
 assert d.signed_distance_um.notna().equals(d.exact_distance_status.eq('estimable'))
 assert (d.signed_distance_um[d.lesion_location=='inside'].dropna()<0).all() and (d.signed_distance_um[d.lesion_location=='outside'].dropna()>0).all()
 bounded=d.abs_distance_upper_bound_um.notna();assert (d.abs_distance_lower_bound_um[bounded]<=d.abs_distance_upper_bound_um[bounded]+1e-6).all()
 for w in [100,200,300]:
  near=d[f'distance_band_{w}um'].isin(['Inner margin','Outer margin']);far=d[f'distance_band_{w}um'].isin(['Interior (>w confirmed)','Outside (>w confirmed)']);assert (d.abs_distance_upper_bound_um[near]<=w+1e-6).all();assert (d.abs_distance_lower_bound_um[far]>w).all()
 for f in ['lesion_rois_um.geojson','lesion_outer_boundaries_um.geojson']:assert (p/f).read_bytes()==(OLD/'by_dataset'/r.dataset_id/sec/f).read_bytes()
 for f in ['01_location_HE','02_bands','03_distance_reliability']:
  assert (p/(f+'.pdf')).exists()
  with Image.open(p/(f+'.png')) as img:assert min(img.info['dpi'])>299
 core=pd.read_csv(p/'spot_labels.tsv.gz',sep='\t');core['band_status']=core.band_status_200um;D.append(core);sm=pd.read_csv(p/'small_ROI_review.tsv',sep='\t');sm['HE_foreground_median']=sm.lesion_roi_id.map(d.groupby('lesion_roi_id').HE_foreground_fraction.median());sm['image_review_status']='computational_image_context; no pathologist confirmation';sm.to_csv(p/'small_ROI_review.tsv',sep='\t',index=False);small.append(sm);gaps.append(pd.read_csv(p/'internal_gap_assessment.tsv',sep='\t'));checks.append({'section_id':sec,'barcodes_original_values_ROI_edges_bands':'pass'})
s=pd.DataFrame(S);s.to_csv(OUT/'tables/section_summary.tsv',sep='\t',index=False);pd.concat(D,ignore_index=True).to_parquet(OUT/'tables/all_spot_labels_v3_1.parquet',index=False);pd.concat(D,ignore_index=True).to_csv(OUT/'tables/all_spot_labels_v3_1.tsv.gz',sep='\t',index=False);pd.concat(small,ignore_index=True).to_csv(OUT/'tables/small_ROI_review.tsv',sep='\t',index=False);pd.concat(gaps,ignore_index=True).to_csv(OUT/'tables/internal_gap_assessment.tsv.gz',sep='\t',index=False);pd.DataFrame(checks).to_csv(OUT/'tables/checks.tsv',sep='\t',index=False)
# Only reaggregate existing v3 section contrasts. No v3.1 scores or comparisons are computed.
k=pd.read_csv(OLD/'tables/section_paired_contrasts.tsv',sep='\t');keys=['dataset_id','analysis_unit_uid','contrast'];missing=k.analysis_unit_uid.isna();k.loc[missing,['dataset_id','section_id','analysis_unit_uid','analysis_unit_type','contrast']].to_csv(OUT/'tables/missing_unit_identity.tsv',sep='\t',index=False);k.loc[missing,'analysis_unit_uid']=k.loc[missing,'dataset_id']+'::unresolved_section::'+k.loc[missing,'section_id'];k.loc[missing,'analysis_unit_type']='UNRESOLVED_SECTION_NOT_PATIENT';vals=[c for c in k if c.endswith('_median_difference')];u=k.groupby(keys,dropna=False)[vals].mean();g=k.groupby(keys,dropna=False)
u['analysis_unit_type_metadata']=g.analysis_unit_type.apply(lambda a:'; '.join(sorted(set(a.dropna()))));u['identity_status_conflict']=g.analysis_unit_type.nunique()>1;u['source_sections']=g.section_id.apply(lambda a:'; '.join(sorted(set(a))));u['n_sections']=g.size();u['source_region_version']='v3_original_regions';u['aggregation']='unweighted_mean_of_existing_v3_section_differences';u=u.reset_index();assert not u.duplicated(keys).any();u.to_csv(OUT/'tables/v3_original_regions_unit_reaggregation.tsv',sep='\t',index=False)
old_u=pd.read_csv(OLD/'tables/analysis_unit_paired_contrasts.tsv',sep='\t');count=old_u.groupby(keys).size().rename('old_rows').reset_index();audit=count.merge(u[keys+['n_sections','identity_status_conflict','analysis_unit_type_metadata','source_sections']],on=keys);audit['new_rows']=1;audit['source_region_version']='v3_original_regions';audit.to_csv(OUT/'tables/unit_key_correction.tsv',sep='\t',index=False)
# Analytic geometry checks: contour zero, signs, absent-unknown exact bounds, and outside unknown extension.
square=np.array([[[0,0],[1000,0]],[[1000,0],[1000,1000]],[[1000,1000],[0,1000]],[[0,1000],[0,0]]]);dist=nearest(np.array([[500,500],[1200,500],[0,500]]),square)[0];assert np.allclose(dist,[500,200,0]);assert np.allclose(nearest(np.array([[1200,500]]),np.array([[[1100,0],[1100,1000]]]))[0],100)
rows=[]
for z in s.itertuples():
 base=f'by_dataset/{z.dataset_id}/{z.section_id}';links=' · '.join(f'<a href="{base}/{f}.png">{label}</a>' for f,label in [('01_location_HE','位置/H&E'),('02_bands','100/200/300μm分带'),('03_distance_reliability','距离/可靠性')]);
 if z.section_id in CONFIG['pilot']:links+=f' · <a href="{base}/04_v3_v3_1_comparison.png">v3对照</a> · <a href="{base}/05_local_review.png">局部复核</a>'
 rows.append(f'<tr><td>{z.section_id}</td><td>{z.exact_restored}</td><td>{z.band_only_restored}</td><td>{z.band_200_unresolved}</td><td>{links}</td></tr>')
(OUT/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>BLCA v3.1</title><style>body{font:15px sans-serif;margin:30px}td{padding:8px;border-bottom:1px solid #ddd}</style><h1>病灶区域 v3.1 — 63张</h1><p>位置、固定宽度分带、精确距离可靠性分别展示。ROI和有效外缘保持v3；斜线为分带未定。</p><p><a href="RESULTS.md">结果说明</a> · <a href="METHODS.md">方法</a> · <a href="tables/all_spot_labels_v3_1.tsv.gz">全队列标签</a></p><table><tr><th>切片</th><th>精确距离恢复</th><th>仅分带恢复</th><th>200μm分带未定</th><th>图</th></tr>'+''.join(rows)+'</table>')
b=s.set_index('section_id').loc['GSE238145_B16'];sm=pd.concat(small);corr=audit[audit.old_rows>1]
report=f'''# v3.1 实际结果

完成{len(s)}/63张、{s.n_spots.sum():,}个原spot，全部保留原barcode、上游数值、主ROI和有效外缘。新结果入口：[全队列浏览](index.html)。没有开展新的MP10/TME比较或方法对照。

位置图统一内红外蓝，不再把距离不可估计显示成另一种组织区域。分带图独立区分内侧/外侧Margin，斜线表示分带未定；连续距离使用RdBu，负值红、正值蓝。可靠性图另示截断、内部缺损、未知缺口、路径遮挡及无有效外缘。

| 量 | 全队列 | B16 |
|---|---:|---:|
| v3精确距离点 | {s.exact_v3.sum():,} | {b.exact_v3:,} |
| v3.1精确距离点 | {s.exact_v3_1.sum():,} | {b.exact_v3_1:,} |
| 从v3删失恢复精确距离 | {s.exact_restored.sum():,} | {b.exact_restored:,} |
| 因全未知区域模型重新删失 | {s.exact_newly_censored.sum():,} | {b.exact_newly_censored:,} |
| 原未分带、现仅恢复200μm分带 | {s.band_only_restored.sum():,} | {b.band_only_restored:,} |
| 200μm分带可判断 | {s.band_200_determined.sum():,} | {b.band_200_determined:,} |
| 200μm分带未定 | {s.band_200_unresolved.sum():,} | {b.band_200_unresolved:,} |

“仅恢复分带”要求v3原分带未定、v3.1精确距离仍NA，但上下界足够判断阈值。无支持病灶和位置未定另存，不归入分带可判断数量。

B16原浅色区域主要是精确距离删失的显示结果。145个组织内部封闭缺位（含ROI外组织孔洞）满足本版技术缺位假设，允许几何路径穿过，不补表达/spot；原内部浅粉岛中的部分位置可恢复分带。外侧使用整个未知空间约束，因此部分原精确值转为删失；不以恢复指定数量为目标。边缘附近无法排除隐藏前沿的位置仍有斜线。

B16的L4（5点）位于左侧截断附近且无有效外缘；L6（1点）位于下方边缘，证据不足以证明生成伪影。两者保留待复核，不声称是独立肿瘤灶。其余非主体ROI也列入表。BC01碎片小ROI同样保留候选。全队列未实施ROI删除/合并：改变spot及原有效边段均为0。小ROI清单包含barcode、域来源、SpaCET/CNV、有效边数及复核图位置。

修复按状态文字拆分分析单位的问题：键改为dataset_id + analysis_unit_uid + contrast，同键切片差值仍取原有等权均值。发现{len(corr)}组原键重复，原单位表{len(old_u)}行合并为{len(u)}行。B5/B18的GSE238145::P1合并为一行，身份状态冲突作为元数据保留。还发现GSE285715的T1/T2/T3缺失analysis_unit_uid，v3曾因NA分组被合并；本版保留三张独立的未核验切片单位，不声称它们属于同一患者或三个已核验患者，缺失身份单独列出。结果文件明确标注 **v3_original_regions**，不是v3.1生物学结果。

交付：`tables/section_summary.tsv`、`all_spot_labels_v3_1.parquet/tsv.gz`、`small_ROI_review.tsv`、`internal_gap_assessment.tsv.gz`、`v3_original_regions_unit_reaggregation.tsv`和`unit_key_correction.tsv`。每片完整保留原数值的parquet在各dataset目录，核心图均有300dpi PNG/PDF。几何/单位/符号、原值和barcode、ROI边段不变、上下界与分带一致、单位键唯一性检查通过。

[看B16对照](by_dataset/GSE238145/GSE238145_B16/04_v3_v3_1_comparison.png) · [看B16局部](by_dataset/GSE238145/GSE238145_B16/05_local_review.png)。
'''
(OUT/'RESULTS.md').write_text(report);print(s.sum(numeric_only=True).to_string());print(corr.to_string(index=False))
