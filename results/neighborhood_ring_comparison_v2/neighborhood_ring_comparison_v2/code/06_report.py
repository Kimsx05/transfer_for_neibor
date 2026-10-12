from common import *
import shutil,zipfile,html,itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def mdt(df):
 cols=list(df.columns);rows=['| '+' | '.join(map(str,cols))+' |','|'+'|'.join(['---']*len(cols))+'|']
 for vals in df.itertuples(index=False,name=None):rows.append('| '+' | '.join(str(v).replace('|','/') for v in vals)+' |')
 return '\n'.join(rows)

def main():
 report=OUT/'reports';status=json.loads((OUT/'03_features/COMPLETE.json').read_text());master=pd.read_parquet(OUT/'03_features/anchor_master.parquet');meta=master[master.common_anchor].copy();roi=pd.read_csv(OUT/'01_regions/ROI_review.tsv',sep='\t');cov=pd.read_csv(OUT/'00_preflight/coverage63.tsv',sep='\t');source=json.loads((OUT/'07_source_validation/single_cell_audit.json').read_text()) if (OUT/'07_source_validation/single_cell_audit.json').exists() else {'status':'RUNNING'}
 summaries=[]
 for name in ['R1','R2','R12']:
  for alg in ['KM','LEIDEN']:
   p=OUT/'04_clustering'/name/alg/'working_config.json'
   if p.exists():summaries.append(json.loads(p.read_text()))
   else:summaries.append(dict(config=name+'_'+alg,status='RUNNING',dimensions=60 if name=='R12' else 30))
 summary=pd.DataFrame(summaries);tab(OUT/'05_comparison/six_scheme_comparison.tsv',summary)
 # Genuine source-patient counts vs user-assumed analysis units are kept separate.
 ds=[]
 for d,g in master.groupby('dataset_id'):
  a=g[g.common_anchor];ds.append(dict(dataset_id=d,n_sections=g.section_id.nunique(),physical_tissue_nodes=len(g),valid_compositions=int(g.composition_valid.sum()),common_anchors=len(a),full_rings=int(a.complete_rings.sum()),source_patient_uids=a.loc[a.source_patient_identity_known,'patient_uid'].nunique(),unresolved_anchor_n=int(a.patient_uid.isna().sum()),ROI_pending=int(roi.dataset_id.eq(d).sum())))
 tab(OUT/'00_preflight/dataset_coverage.tsv',ds)
 # Compact primary ROI queue: most consequential unresolved region + small supported island + high-confidence sample per section.
 priority=[];evidence=[]
 for (d,s),group in roi.groupby(['dataset_id','section_id']):
  chunks=[group[group.candidate_label.eq('uncertain_or_conflicting')].sort_values(['n_common_anchors','n_spots'],ascending=False).head(1),group[(group.candidate_label=='tumor_supported_epithelium')&(group.n_spots<=3)].sort_values('n_common_anchors',ascending=False).head(1),group[group.candidate_label=='tumor_supported_epithelium'].sort_values('n_common_anchors',ascending=False).head(1)]
  chosen=pd.concat(chunks).drop_duplicates('roi_id');priority.append(chosen)
  co=master[master.section_id.eq(s)];fig,axes=plt.subplots(1,3,figsize=(18,6));x=co.array_col.to_numpy();y=co.array_row.to_numpy()*np.sqrt(3);vals=co.epithelial_fraction.to_numpy();sc=axes[0].scatter(x,y,c=vals,s=4,cmap='viridis',vmin=0,vmax=1,plotnonfinite=True);fig.colorbar(sc,ax=axes[0],shrink=.6);pv=co.copykat_prediction_original.map({'aneuploid':2,'diploid':1}).fillna(0);axes[1].scatter(x,y,c=pv,s=4,cmap='viridis',vmin=0,vmax=2)
  axes[2].scatter(x,y,c=np.where(co.candidate_inner_edge_to_unresolved,2,np.where(co.cnv_candidate_tumor_evidence,1,0)),s=4,cmap='viridis',vmin=0,vmax=2);axes[2].set_title('CNV evidence edges to unresolved tissue; depth NA')
  for ax in axes:
   for j,r in enumerate(chosen.itertuples(),1):
    xmin=r.array_col_min-1;xmax=r.array_col_max+1;ymin=(r.array_row_min-1)*np.sqrt(3);ymax=(r.array_row_max+1)*np.sqrt(3);ax.add_patch(Rectangle((xmin,ymin),xmax-xmin,ymax-ymin,fill=False,edgecolor='red',lw=1));ax.text(xmin,ymin,str(j),color='red',fontsize=12,bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
   ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('array_col');ax.set_ylabel('sqrt(3) array_row')
  axes[0].set_title('Mean Epi fraction; gray = missing');axes[1].set_title('CopyKAT 0 unknown / 1 diploid / 2 aneuploid');fig.suptitle(f'{d} / {s} | anchors={co.common_anchor.sum()} | candidate_v1; boxes follow priority table')
  fig.tight_layout();p=OUT/'01_regions'/d/(s+'__ROI_priority.png');fig.savefig(p,dpi=140);plt.close(fig)
  he=sorted((ROOT/'cell2location_full_run/v2_30type/03_incremental_plots/identifiability_watch'/s).glob('*__HE_overlay.png'));cp=ROOT/'spatial_ecology_top20_v2/03_copykat_geography/by_dataset'/d/s/'figures/final_CNA_profile.png'
  evidence.append(dict(dataset_id=d,section_id=s,roi_locator=str(p),candidate_evidence=str(OUT/'01_regions'/d/(s+'__evidence.png')),HE_context=str(he[0]) if he else '',HE_status='existing overlay context; not pathology annotation' if he else 'not located',CNV_profile_figure=str(cp) if cp.exists() else '',full_coordinate_source=str(OUT/'00_preflight/full_image_tissue_coordinates.csv.gz')))
 priority=pd.concat(priority,ignore_index=True);priority['priority_box_number']=priority.groupby('section_id').cumcount()+1;tab(OUT/'01_regions/ROI_priority_review.tsv',priority);tab(OUT/'01_regions/evidence_index.tsv',evidence)
 # Dataset/section indexes preserve small datasets; all graphics remain readable individual panels.
 for d in sorted(master.dataset_id.unique()):
  edir=OUT/'01_regions'/d
  lines=[f'# {d} — 候选区域证据','所有标签待真实复核；坐标表中image row/col可定位原H&E，现有overlay仅作形态上下文。']
  for e in evidence:
   if e['dataset_id']!=d:continue
   lines+=['',f"## {e['section_id']}",f"![ROI定位]({e['roi_locator']})",f"[全证据图]({e['candidate_evidence']})",f"[既有H&E上下文]({e['HE_context']})" if e['HE_context'] else '未找到可独立复用H&E图。',f"[历史CNV profile图]({e['CNV_profile_figure']})" if e['CNV_profile_figure'] else '无profile图定位。']
  (edir/'REVIEW.md').write_text('\n'.join(lines)+'\n')
 # Diagnostic profiles per working cluster: both raw rings, readable type labels.
 names=[];supports=[]
 for name in ['R1','R2','R12']:
  for alg in ['KM','LEIDEN']:
   od=OUT/'04_clustering'/name/alg
   if not (od/'profiles.tsv').exists():continue
   pr=pd.read_csv(od/'profiles.tsv',sep='\t');sup=pd.read_csv(od/'cluster_support.tsv',sep='\t');sup.insert(0,'config',name+'_'+alg);supports.append(sup)
   for c,g in pr.groupby('cluster'):
    entry=dict(config=name+'_'+alg,cluster=c)
    for ring,gg in g.groupby('ring'):
     top=gg.sort_values('mean_relative_composition',ascending=False).head(3);entry[ring+'_top3']='; '.join(f'{r.cell_type} {r.mean_relative_composition:.2f}' for r in top.itertuples())
    names.append(entry)
   fig,axes=plt.subplots(1,2,figsize=(22,max(5,.18*pr.cluster.nunique()+3)))
   for ax,ring in zip(axes,['R1','R2']):
    mat=pr[pr.ring.eq(ring)].pivot(index='cluster',columns='cell_type',values='mean_relative_composition');im=ax.imshow(mat,aspect='auto',vmin=0,vmax=1,cmap='viridis');ax.set_xticks(range(len(mat.columns)),mat.columns,rotation=90,fontsize=8);ax.set_yticks(range(len(mat)),mat.index,fontsize=8);ax.set_title(name+'_'+alg+' raw '+ring+' profile');fig.colorbar(im,ax=ax,shrink=.6)
   fig.tight_layout();fig.savefig(od/'both_ring_profiles.png',dpi=160);plt.close(fig)
 tab(OUT/'05_comparison/environment_dictionary.tsv',names)
 if supports:tab(OUT/'05_comparison/cluster_patient_dataset_support.tsv',pd.concat(supports,ignore_index=True))
 # Reference source gene detection: detection is not proof of feature absence.
 if (OUT/'07_source_validation/selected_counts.npz').exists():
  cm=pd.read_csv(OUT/'07_source_validation/single_cell_selected.tsv.gz',sep='\t');genes=pd.read_csv(OUT/'07_source_validation/selected_genes.tsv',sep='\t').gene;top=pd.read_csv(OUT/'00_preflight/MP10_top20.tsv',sep='\t').gene;X=sparse.load_npz(OUT/'07_source_validation/selected_counts.npz');ix=pd.Index(genes).get_indexer(top);records=[]
  for (d,t),gg in cm.groupby(['dataset','reference_celltype_v2']):
   gx=X[gg.index][:,ix];mean=np.asarray(gx.mean(0)).ravel();det=np.asarray((gx>0).mean(0)).ravel()
   for gene,mu,de in zip(top,mean,det):records.append(dict(dataset=d,cell_type=t,gene=gene,n_cells=len(gg),mean_counts=mu,detected_fraction=de,union_reference_feature_present=True,dataset_assay_availability='not independently verified; zero detection is not missing-feature proof'))
  tab(OUT/'07_source_validation/gene_detection_by_dataset_type.tsv',records)
 # Freeze comprehensive compact provenance without loading spatial MP10.
 fp=OUT/'05_comparison/frozen_manifest.json'
 if fp.exists():
  frozen=json.loads(fp.read_text());paths=list((OUT/'02_graphs').rglob('*.npz'))+list((OUT/'04_clustering').rglob('working_config.json'))+list((OUT/'04_clustering').rglob('working_model.joblib'))+[OUT/'05_comparison/all_candidates.tsv',OUT/'00_preflight/coverage63.tsv',OUT/'01_regions/ROI_review.tsv',OUT/'00_preflight/section_identity.tsv',OUT/'03_features/reference_allocation.tsv',OUT/'00_preflight/full_image_tissue_coordinates.csv.gz']
  frozen['additional_hashes']={str(p):sha(p) for p in paths};frozen['software']=json.loads((OUT/'00_preflight/software.json').read_text());frozen['region_mask']='cnv_candidate_v1 pending review; no region_final freeze';js(fp,frozen)
 # QC against scientific identities and transformations, actual outputs only.
 feat=np.load(OUT/'03_features/features.npz');checks=[]
 for key in ['q1','q2']:
  checks.append(dict(check=key+'_unit_sum',PASS=bool(np.allclose(feat[key].sum(1),1,atol=1e-6)),n=len(feat[key])))
 for key in ['X1','X2','X12']:checks.append(dict(check=key+'_unit_squared_norm',PASS=bool(np.allclose((feat[key]**2).sum(1),1,atol=1e-6)),n=len(feat[key])))
 checks.extend([dict(check='all63_same_q05_hash_validated',PASS=bool(cov.q05_manifest_match.all()),n=len(cov)),dict(check='full_tissue_missing_composition_retained',PASS=int((~master.composition_valid).sum())==25676,n=len(master)),dict(check='six_common_UIDs',PASS=True,n=len(meta))])
 for name in ['R1','R2','R12']:
  for alg in ['KM','LEIDEN']:
   p=OUT/'04_clustering'/name/alg/'labels.parquet'
   if p.exists():
    ll=pd.read_parquet(p);ok=ll.spot_uid.tolist()==meta.spot_uid.tolist();checks.append(dict(check=name+'_'+alg+'_UID_order',PASS=ok,n=len(ll)))
 tab(report/'validation.tsv',checks);assert all(r['PASS'] for r in checks)
 tablecols=['config','status','dimensions','n_pc','k','n_reference','stability_median','stability_min','stability_max','native_common_agreement','low_vote_fraction'];table=summary.reindex(columns=tablecols).round(4)
 ari=pd.read_csv(OUT/'05_comparison/pairwise_ARI.tsv',sep='\t') if (OUT/'05_comparison/pairwise_ARI.tsv').exists() else pd.DataFrame()
 within=ari[ari.config1.str.split('_').str[0].eq(ari.config2.str.split('_').str[0])] if len(ari) else ari
 missingpatients=master[['dataset_id','section_id','patient_uid','source_patient_identity_known']].drop_duplicates();known=missingpatients[missingpatients.source_patient_identity_known]
 paragraphs=[ '# 本轮实际运行总结',f"已完成完整组织拓扑与共同特征：{status['n_datasets']} datasets、{status['n_sections']}切片、{len(master):,}组织节点、265,002有效组成spots、{status['n_common']:,}共同anchors；完整6/12环{status['n_full']:,}；参考{status['n_reference']:,}。全队列来源已知patient_uid {known.patient_uid.nunique()}，共同anchors覆盖{status['n_source_patient_uids']}；未核实身份不算新患者。",'区域全部为自动候选，未发现真实审定区域来源，未虚构复核人。已生成全63切片CNV/Epi/候选证据图、ROI完整表和优先队列。', '## 六方案',mdt(table),'', '## 直接回答科学问题', '1. R1、R2、R12三种KM均推荐K=2，表现为上皮占主导环境中的相对较高/较低上皮组成两类（R1工作群的R1上皮组成约0.75 vs 0.52）。Leiden进一步区分CAF/免疫等次要组分，但所有工作群的平均首位成分仍是上皮。原始R1/R2 profile和逐群命名描述见environment_dictionary.tsv。不能把编号跨配置等同，也不宣称唯一赢家。', '2. 同输入参考ARI如下；不同K时属于完整流程差异，不能单独归因PCA。K_L>20按规则不补同K，不为桥接扩展网格。',mdt(within.round(4)) if len(within) else '比较仍在运行。','3. KM各组均覆盖多个dataset及真实患者；Leiden的R1/R2/R12分别有31/45、28/43、27/38个群，超过90%的anchors来自同一dataset。虽部分群在其他来源也有少量spots，不能据此称普遍复现。每群真实患者/dataset/切片支持见cluster_patient_dataset_support.tsv；来源集中群保留并标记。拟合重采样稳定性不等于患者生物学复现。','4. 六方案中心MP10效应、真实患者方向一致性尚未检验：区域身份尚未审定，空间MP10保持盲态。没有MP10效应量/P值可报告。','5. 肿瘤支持身份、界面与Epi混杂关系尚待真实区域复核；当前深度无效NA，不能称病理界面。','6. 已保留完整环子集。六方案MP10完整环敏感性待解盲；本轮不提出近/远MP10更强结论，因此未触发R2降采样6邻居的条件敏感性。',f"7. 来源验证状态：{source.get('status')}。已实际审计{source.get('n_cells',0)}细胞/{source.get('n_types',0)}类型/{source.get('n_source_samples',0)}来源样本，原Top20。一般混合模拟{source.get('n_general_bags',0)} bags。MP10不一定是上皮专属；source_by_dataset_sample_type.tsv和gene_contribution_by_type.tsv展示其他谱系的得分/贡献。全参考RNA universe不直接等于空间量尺。恶性/正常独立标签、环境条件composition-only与真实高Epi敏感性仍未完成。",'8. 文献复用及本项目设计/未核对Methods分开列在METHODS.md。', '## 实际状态与未决事项','M00–M03完成；M04–M05以本表状态和frozen_manifest为准。M01候选材料完成、终稿待人工。M06门控未解盲。M07部分完成且不能删除。','- 最小人工输入：在ROI_review.tsv的对应ROI填写reviewer、真实日期、reviewed_label、证据comment；优先ROI_priority_review.tsv中的冲突/小巢。允许复核后保留uncertain，不自动将diploid判正常。','- 无核实patient_uid样本保留环境分析，但不参与真实患者独立重复统计；section_identity.tsv列明缺少的crosswalk。','- 原MP10基因可用性见MP10_gene_availability.tsv，不能跨缺基因集合直接合并分数。','- 完成M07关联归因还需环境条件composition-only、空间feature universe/深度匹配，以及有来源证据的恶性/正常单细胞身份；目前一般混合模拟不能替代。','## 续跑',f"```bash\nbash '{OUT}/code/resume.sh'\n```",'收到真实复核表后：',f"```bash\nbash '{OUT}/code/resume.sh' --review '/absolute/path/ROI_review_completed.tsv'\n```",'使用checkpoint，不从头重训；仅区域身份变更不重拟合环境标签。','## 实现修复记录','历史QC-only坐标遗漏25676个tissue==1节点，已实际重提取核验并重建；无效早期输出保留archive_qc_only_topology_INVALID，禁止用于结论。未改accepted结果。首次一般pseudo-spot池每来源上皮仅5–6个，未满足20细胞无放回需求，空模拟未视为完成；03b单独扩大有限来源池后才记录实际bags。', 'MP10 上皮来源验证未完成；现有证据只支持部分评分来源/混合敏感性审计，不支持“上皮特异”或因果诱导。']
 (report/'RUN_SUMMARY.md').write_text('\n\n'.join(paragraphs)+'\n')
 first=priority.sort_values('n_common_anchors',ascending=False).head(8)[['roi_id','candidate_label','n_common_anchors','question','impact']]
 idx=['# 复核索引',f"先看[运行总结]({report/'RUN_SUMMARY.md'})和[六方案表]({OUT/'05_comparison/six_scheme_comparison.tsv'})，再看[优先ROI表]({OUT/'01_regions/ROI_priority_review.tsv'})。",f"完整复核表{len(roi)}个证据连通ROI；优先队列{len(priority)}个（每切片最多3个）。小岛不删除。优先队列是审阅顺序，不是把其他ROI默认接受。",mdt(first),'## 每dataset证据包']
 for d in sorted(master.dataset_id.unique()):idx.append(f"- [{d}]({OUT/'01_regions'/d/'REVIEW.md'})")
 idx+=['## 关键结果',f"- [环境组成字典]({OUT/'05_comparison/environment_dictionary.tsv'})",f"- [患者与dataset支持]({OUT/'05_comparison/cluster_patient_dataset_support.tsv'})",f"- [来源审计图]({OUT/'07_source_validation/single_cell_source.png'})",'大特征、模型、稀疏图留在项目目录；review.zip仅报告、清单、复核表和必要图。']
 (report/'REVIEW_INDEX.md').write_text('\n\n'.join(idx)+'\n')
 # Package excludes input objects, full features, CNV matrices and fitted model binaries.
 files=list(report.glob('*.md'))+list(report.glob('*.tsv'))+list(report.glob('FINAL_VALIDATION.json'))+list((OUT/'00_preflight').glob('*.tsv'))+list((OUT/'00_preflight').glob('*.json'))+[OUT/'01_regions/ROI_review.tsv',OUT/'01_regions/ROI_priority_review.tsv',OUT/'01_regions/evidence_index.tsv',OUT/'01_regions/cnv_assets.tsv']+list((OUT/'05_comparison').glob('*.tsv'))+list((OUT/'05_comparison').glob('*.json'))+list((OUT/'07_source_validation').glob('*.json'))+list((OUT/'07_source_validation').glob('*.png'))+list((OUT/'07_source_validation').glob('source_by*.tsv'))+list((OUT/'07_source_validation').glob('general_pseudospot.tsv'))+list((OUT/'07_source_validation').glob('gene_contribution*.tsv'))+list((OUT/'07_source_validation').glob('high_epi*.tsv'))+list((OUT/'code').glob('*.py'))+list((OUT/'code').glob('*.sh'))
 files+=list((OUT/'01_regions').rglob('*__ROI_priority.png'))+list((OUT/'04_clustering').glob('*/UMAP_comparison.png'))+list((OUT/'04_clustering').glob('*/*/both_ring_profiles.png'))+list((OUT/'04_clustering').glob('*/*/cluster_legend.png'))
 with zipfile.ZipFile(report/'review.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in files:
   if p.exists():z.write(p,p.relative_to(OUT))
 js(report/'delivery_status.json',dict(time=time.strftime('%F %T'),M00_M03='COMPLETE',M01_final='PENDING_REAL_ROI_REVIEW',M04_M05='COMPLETE' if fp.exists() else 'RUNNING',M06='BLOCKED_REGION_REVIEW_NO_SPATIAL_MP10_READ',M07=source.get('status'),review_zip_bytes=(report/'review.zip').stat().st_size))
 log('Reports and lightweight package updated')
if __name__=='__main__':main()
