import os,sys,json,gzip,zipfile,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'mp10_cnv_regions_interface_v1';E=ROOT/'spatial_ecology_top20_v2';sys.path.insert(0,str(E/'runtime'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import pandas as pd,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
plt.rcParams.update({'font.size':9,'axes.titlecolor':'black','pdf.fonttype':42})
def read(p):return pd.read_csv(O/p,sep='\t',low_memory=False)
def tab(d,p):pd.DataFrame(d).to_csv(O/p,sep='\t',index=False,na_rep='NA')
cov=read('input_review/section_coverage.tsv');mp=read('target_profiles/region_mp10_summary.tsv');geo=read('input_review/coordinate_alignment.tsv');effects=read('target_profiles/patient_section_effects.tsv');eff=effects[effects.level=='section'];pat=read('target_profiles/verified_patient_effects.tsv');profiles=read('target_profiles/cross_dataset_summary.tsv');extra=read('input_review/additional_platform_coverage.tsv');tab(extra.drop_duplicates(['dataset_id','section_id']),'input_review/additional_platform_coverage.tsv')
for i,r in cov.iterrows():
 z=mp[(mp.section_id==r.section_id)&(mp.scope=='distance_band')];unavailable=z.loc[z.n_positions.eq(0),'group'].tolist();reason=r.missing_reason
 if r.n_epi_anchors==0:reason+='; no positions reach frozen Epi>=0.25, no region/interface can be defined'
 if unavailable:reason+='; not evaluable bands: '+','.join(unavailable)
 if geo.set_index('section_id').loc[r.section_id,'max_fit_error_dNN']>.1:reason+='; affine grid residual >0.1 d_NN; boundary localization less precise'
 cov.loc[i,'missing_reason']=reason
 cov.loc[i,'status']='completed_no_epi_region' if r.n_regions==0 else 'completed_with_limitations'
tab(cov,'input_review/section_coverage.tsv')
idp=E/'patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv';tab(pd.read_csv(idp,sep='\t'),'input_review/section_identity_crosswalk.tsv');genes=pd.read_csv(E/'00_inputs/gene_coverage.tsv',sep='\t');tab(genes[genes.MP=='EpiMP10'],'input_review/MP10_gene_coverage.tsv')
# Two label systems describe different spatial objects; overlap is descriptive, on their common accepted positions only.
ol=read('niche_review/local_R12_overlap.tsv');ol=ol[~ol.R12_niche.isin(['Not_in_anchor_universe','Unassigned'])];ns=[]
for ds,g in ol.groupby('dataset_id'):
 ct=g.groupby(['local_niche_id','R12_niche']).n_positions.sum().unstack(fill_value=0).to_numpy(float);prob=ct/ct.sum();row=prob.sum(1);col=prob.sum(0);expected=row[:,None]*col[None,:];nz=prob>0;mi=float((prob[nz]*np.log(prob[nz]/expected[nz])).sum());h1=-(row[row>0]*np.log(row[row>0])).sum();h2=-(col[col>0]*np.log(col[col>0])).sum();ns.append(dict(dataset_id=ds,n_shared_positions=int(ct.sum()),n_sections=g.section_id.nunique(),n_local_niches=ct.shape[0],n_R12_niches=ct.shape[1],normalized_mutual_information=mi/np.sqrt(h1*h2) if h1*h2>0 else np.nan,local_to_R12_modal_fraction=ct.max(1).sum()/ct.sum(),R12_to_local_modal_fraction=ct.max(0).sum()/ct.sum(),interpretation='descriptive shared accepted anchors; algorithms/resolutions differ; no forced bijection'))
ns=pd.DataFrame(ns);tab(ns,'niche_review/niche_system_overlap_summary.tsv')
sens=read('interfaces/representative_2w_sensitivity.tsv');se=[]
for (ds,sec,scale),g in sens.groupby(['dataset_id','section_id','bandwidth_multiplier']):
 q=g.set_index('distance_band')
 for l,r in [('outer_neighborhood','inner_edge'),('inner_edge','inner_core')]:
  for metric in ['MP10_mean','Treg_all_modeled_fraction','Treg_T_NK_fraction','Epi_fraction']:
   se.append(dict(dataset_id=ds,section_id=sec,bandwidth_multiplier=scale,contrast=l+' minus '+r,metric=metric,n_left=int(q.loc[l,'n_positions']),n_right=int(q.loc[r,'n_positions']),effect=q.loc[l,metric]-q.loc[r,metric] if q.loc[l,'n_positions']>0 and q.loc[r,'n_positions']>0 else np.nan))
se=pd.DataFrame(se);tab(se,'interfaces/representative_bandwidth_effects.tsv');wide=se.pivot(index=['dataset_id','section_id','contrast','metric'],columns='bandwidth_multiplier',values='effect').reset_index();wide.columns=[str(x) for x in wide.columns];wide['both_evaluable']=wide['1'].notna()&wide['2'].notna();wide['direction_changed']=(wide['1']*wide['2']<0)&wide.both_evaluable;tab(wide,'interfaces/representative_bandwidth_direction.tsv')
# Re-render paired summaries with explicit missingness, readable labels and per-panel coverage.
metrics=['MP10_mixed_AUCell','MP10_epi_component_existing','MP10_linear_epi_residual','major:T_NK:all_modeled_fraction','fine:Treg:all_modeled_fraction','fine:Treg:T_NK_fraction','major:Endothelial:all_modeled_fraction','major:Fibroblast:all_modeled_fraction','major:Mural:all_modeled_fraction'];labels=['Mixed-spot MP10 (AUCell)','Existing Epi-component MP10','MP10 residual vs Epi proportion','T/NK / all types','Treg / all types','Treg / T/NK','Endothelial / all types','Fibroblast / all types','Pericyte + SMC / all types'];datasets=sorted(cov.dataset_id.unique());colors=plt.cm.tab10(np.arange(7))
for contrast in eff.contrast.unique():
 f,axs=plt.subplots(3,3,figsize=(16,12))
 for ax,metric,label in zip(axs.ravel(),metrics,labels):
  xt=[]
  for i,ds in enumerate(datasets):
   z=eff[(eff.dataset_id==ds)&(eff.contrast==contrast)&(eff.metric==metric)].dropna(subset=['effect']);v=z.effect.to_numpy();ax.scatter(i+np.linspace(-.18,.18,len(v)),v,s=13,color=colors[i],alpha=.7);pv=pat[(pat.dataset_id==ds)&(pat.contrast==contrast)&(pat.metric==metric)].effect.dropna();ax.scatter(np.repeat(i+.27,len(pv)),pv,s=18,marker='D',facecolors='none',edgecolors='black')
   if len(v):ax.plot([i-.2,i+.2],[np.median(v)]*2,color='black',lw=2)
   xt.append(ds.replace('__UNVERIFIED','')+f'\nn={len(v)}')
  if not eff[(eff.contrast==contrast)&(eff.metric==metric)].effect.notna().any():ax.text(.5,.5,'Not evaluable\nExisting component scores are restricted\nto Epi anchors; no outer-neighborhood scores',transform=ax.transAxes,ha='center',va='center',fontsize=10)
  ax.axhline(0,c='grey',lw=.7);ax.set_xticks(range(7),xt,rotation=35,ha='right',fontsize=7);ax.set_title(label,loc='center');ax.set_ylabel('within-section paired difference')
 f.suptitle(contrast.replace('_',' ')+' | circles: sections; diamonds: verified patients; bar: section median',color='black');f.tight_layout(rect=[0,0,1,.97]);name='paired_effects_'+contrast.replace(' ','_');f.savefig(O/'figures'/f'{name}.pdf',transparent=True,bbox_inches='tight');f.savefig(O/'figures'/f'{name}.png',dpi=300,transparent=True,bbox_inches='tight');plt.close(f)
# Human-readable data dictionary accompanies the cross-session interface.
(O/'interfaces/INTERFACE_README.md').write_text('''# 只读位置接口

主键：dataset_id + section_id + position_id（保留完整barcode后缀）。265,002行，每个有效观测位置一行；外侧不按多个region重复计数。

- region_id：最近有效界面所属Epi-rich region。position_region_id：位置实际所属区域，外侧为Outside；无界面时相应键缺失或保留位置区域。
- region_label / region_evidence / region_confidence：组成支持区域，不能解释为病理确诊肿瘤。
- boundary_id：对应boundary_summary.tsv.gz与outline GeoJSON中的单条有效线段。
- distance_to_boundary：非负欧氏最短距离；signed_distance：内部负、外侧正；signed_distance_available=False时距离NA。
- distance_unit / coordinate_unit：本片原图像像素；coordinate_reference为原图像绝对路径。不同图像像素不是同一物理单位。
- distance_dNN：signed_distance / 本片中位最近邻距，适合相对分辨率展示，不是微米。
- distance_band：inner_core、inner_edge、outer_neighborhood、far_tissue或not_evaluable。
- evaluation_status=evaluable：可进入主分带描述；censored_near_support_edge距离可读但主分析排除；no_valid_interface_in_component不能当作0。
- distance_to_support_edge：距有效观测支持边距离，包含QC缺失与组织/截断不明边。
- patient_id仅为来源可核实ID；record_id包含匿名/假定分析单位，不能统称独立患者。

边界表中eligible_interface=True才是区域内外界面；missing_QC_observation_edge为完整tissue网格中有位置但本次QC观察缺失；tissue_or_capture_edge_unresolved无法区分组织外缘和捕获取材截断。support_topology另分outer_envelope与interior_hole。参照methods.md使用，不把缺失当作阴性。
''')
report=(O/'report.md').read_text();summary_par='''## 阅读重点

- 主比较为同片“外侧邻域−内侧边缘”。混合spot MP10的dataset中位差在6/7个dataset为负，GSE171351为正；这是混合组成层面的描述，不能据此认定上皮细胞内状态降低。
- Inhouse的7/7张可比较切片均显示外侧混合MP10较低、Treg在全30类中的占比较高；但Treg在T/NK内的比例仅3/7为正。跨dataset也有类似分母依赖，不能概括成统一的Treg特异富集。
- 已有Epi-component分数仅覆盖Epi anchors，外侧邻域没有对应分数，故跨界面Epi-component比较全部不可评估；core/inner-edge内可用比较另表保留。空白不是零效应。
- Fibroblast全类型占比的外侧−内侧dataset中位差在7/7为正；EC为6/7，GSE238145为负。区域本身由Epi比例定义，其他谱系比例上升可能含组成闭合效应；模型丰度与谱系内分母表同时提供。

'''
report=report.replace('## 证据覆盖',summary_par+'## 证据覆盖',1)
lim='GSE238145_B13的4,881个位置均未达到冻结Epi门槛，因此没有本轮Epi-rich区域或界面；该片仍保留完整输入、CNV、MP10和niche记录。Inhousedata_T4384654和GSE171351_D1无可评估core；不为补齐分带改阈值。'
report=report.replace('## 不一致与限制','## 不一致与限制\n\n'+lim,1)
report=report.replace('## Niche 对照','## Niche 对照\n\n同位置local与R1/R2标签的归一化互信息范围为'+f'{ns.normalized_mutual_information.min():.3f}–{ns.normalized_mutual_information.max():.3f}，两者不能直接互换；这是已赋值共同anchors上的描述性重叠，不是独立样本验证。详见niche_system_overlap_summary.tsv。',1)
report=report.replace('## 读取接口',f'代表片1w与2w敏感性在{int(wide.both_evaluable.sum())}个可比较的section×contrast×指标组合中，{int(wide.direction_changed.sum())}个方向改变；结果见representative_bandwidth_direction.tsv。带宽结果未用于选择主设置。\n\n## 读取接口',1)
(O/'report.md').write_text(report)
methods=(O/'methods.md').read_text();methods+='\n\n最终独立复核：运行 scripts/independent_checks.py，用独立Shapely点到线距离抽查62张有界面切片共744个位置，核查符号与mask内外；全部265,002位置原x/y和raw MP10逐值相符。最后运行 scripts/finalize_review.py，补充明确缺失说明、niche重叠和1w/2w方向表并封装review。\n';(O/'methods.md').write_text(methods)
# Plot annotations audit, including native categorical codes in legacy overview.
(O/'figures/figure_notes.md').write_text('所有空间图使用原图像像素，黑线为Epi-rich有效界面，不是病理肿瘤边界。Overview中region二元码0=outside，1=Epi-rich；MP10状态0=non-anchor、1=relative low、2=ties、3=relative high。Niche显示数字码按各片display_codes表解释。Diagnostics中CNV码0=diploid、1=aneuploid、2=undefined；profile/eligible码0=no、1=yes。meta-niche永远是candidate。Fine目标全17类均为c2l模拟估计比例。\n')
for name in ['paired_effects_outer_neighborhood_minus_inner_edge','paired_effects_inner_edge_minus_inner_core']:
 im=Image.open(O/'figures'/f'{name}.png');im.thumbnail((1600,1200));bg=Image.new('RGB',im.size,'white');bg.paste(im,mask=im.getchannel('A'));bg.save(O/'logs'/f'{name}_preview.jpg',quality=90)
fi=read('figures/figure_index.tsv');assert len(fi)==26
for r in fi.itertuples():
 assert (O/r.pdf).stat().st_size>1000
 with Image.open(O/r.png) as im:assert im.mode=='RGBA' and im.info.get('dpi',(0,))[0]>=299
rawcheck=json.loads((O/'logs/source_parity_validation.json').read_text());assert rawcheck['source_xy_unchanged'] and rawcheck['source_raw_MP10_unchanged'];valid=json.loads((O/'logs/validation.json').read_text());valid.update(independent_checks=rawcheck,figure_pairs_verified=len(fi),all_PNG_RGBA_300dpi=True,visual_review=['CanDis_2026_Sample1 overlay','Inhousedata_T4384654 overlay','paired effects outer minus inner'],review_status='complete_with_explicit_evidence_limits');(O/'logs/validation.json').write_text(json.dumps(valid,indent=2))
# Hash scripts and small critical config, never expression mother objects.
tab([dict(path=str(p.relative_to(O)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in list((O/'scripts').glob('*.py'))+[O/'regions/frozen_rules.json']],'logs/code_checksums.tsv')
small=[O/'README.md',O/'report.md',O/'methods.md',O/'figures/figure_index.tsv',O/'figures/figure_notes.md',O/'interfaces/INTERFACE_README.md',O/'logs/validation.json',O/'logs/source_parity_validation.json',O/'input_review/section_coverage.tsv',O/'input_review/dataset_coverage.tsv',O/'input_review/coordinate_alignment.tsv',O/'input_review/additional_platform_coverage.tsv',O/'cnv_review/cnv_availability.tsv',O/'regions/region_parameters.tsv',O/'target_profiles/cross_dataset_summary.tsv',O/'target_profiles/verified_patient_effects.tsv',O/'niche_review/niche_system_overlap_summary.tsv',O/'interfaces/representative_bandwidth_direction.tsv',O/'interfaces/representative_2w_sensitivity.tsv']
small+=list((O/'logs').glob('*preview.jpg'));small+=[O/'figures'/f'{s}_spatial_overview.pdf' for s in ['Inhousedata_T4384654','GSE319536_BC01']];small += [O/'figures'/f'{s}.pdf' for s in ['all_sections_band_coverage','major_composition_by_dataset','continuous_distance_profiles','paired_effects_outer_neighborhood_minus_inner_edge','paired_effects_inner_edge_minus_inner_core']]
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in small:z.write(p,p.relative_to(O))
with zipfile.ZipFile(O/'review.zip') as z:assert z.testzip() is None
manifest=[dict(path=str(p.relative_to(O)),bytes=p.stat().st_size) for p in O.rglob('*') if p.is_file() and p.name!='output_manifest.tsv'];tab(manifest,'output_manifest.tsv');(O/'completion.json').write_text(json.dumps(dict(status='completed_with_evidence_limits',n_sections=63,n_positions=265002,n_regions=1916,n_signed_positions=256788,n_primary_positions=217538,n_figure_pairs=len(fi),review_zip_bytes=(O/'review.zip').stat().st_size,outputs_root=str(O),old_inputs_read_only=True,communications_not_run=True,cnv_not_reinferred=True),indent=2));print((O/'completion.json').read_text())
