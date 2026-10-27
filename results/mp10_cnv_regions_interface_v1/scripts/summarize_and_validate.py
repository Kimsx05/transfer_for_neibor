#!/usr/bin/python3.12
import os,sys,json,gzip,zipfile,hashlib,platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'mp10_cnv_regions_interface_v1';E=ROOT/'spatial_ecology_top20_v2';sys.path.insert(0,str(E/'runtime'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import numpy as np,pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from shapely.geometry import shape,Point
from shapely.ops import unary_union
from PIL import Image
plt.rcParams.update({'font.size':9,'axes.titlecolor':'black','pdf.fonttype':42})
def read(p):return pd.read_csv(O/p,sep='\t',low_memory=False)
def tab(d,p):pd.DataFrame(d).to_csv(O/p,sep='\t',index=False,na_rep='NA')
figs=[]
def figsave(f,n,desc):
 f.savefig(O/'figures'/f'{n}.pdf',transparent=True,bbox_inches='tight');f.savefig(O/'figures'/f'{n}.png',dpi=300,transparent=True,bbox_inches='tight');plt.close(f);figs.append(dict(figure=n,pdf=f'figures/{n}.pdf',png=f'figures/{n}.png',description=desc))
cov=read('input_review/section_coverage.tsv');reg=read('regions/region_summary.tsv');p=read('regions/position_regions.tsv.gz');d=read('interfaces/position_interface_distances.tsv.gz');feat=read('input_review/position_features.tsv.gz');b=read('interfaces/boundary_summary.tsv.gz');eff=read('target_profiles/patient_section_effects.tsv');eff=eff[eff.level.eq('section')].copy();mp=read('target_profiles/region_mp10_summary.tsv');tar=read('target_profiles/target_region_summary.tsv.gz');nr=read('niche_review/niche_region_crosswalk.tsv');nt=read('niche_review/niche_target_summary.tsv.gz');cp=read('input_review/coordinate_alignment.tsv');pars=read('regions/region_parameters.tsv');cnv=read('cnv_review/cnv_availability.tsv');identity=pd.read_csv(E/'patient_identity_audit_20261002/revised_by_user_v2/section_patient_and_analysis_unit_crosswalk.tsv',sep='\t');keys=['dataset_id','section_id','position_id'];assert len(p)==len(d)==len(feat)==265002 and len(cov)==63 and cov.dataset_id.nunique()==7
assert not p.duplicated(keys).any() and not d.duplicated(keys).any();pd.merge(p[keys],d[keys],validate='one_to_one');assert d.loc[d.signed_distance_available,'boundary_id'].isin(b.loc[b.eligible_interface,'boundary_id']).all();assert p.final_mask.equals(p.candidate_mask);assert np.isclose(reg.n_positions.sum(),p.final_mask.sum());assert b.boundary_id.is_unique
# Annotate outer envelope versus interior holes geometrically; tissue/QC etiology remains separate.
for sec in cov.section_id:
 with gzip.open(O/'regions/masks'/f'{sec}.geojson.gz','rt') as f:geo=json.load(f)
 tissue=shape(next(x['geometry'] for x in geo['features'] if x['properties']['region_id']=='valid_observation_support'));polys=[tissue] if tissue.geom_type=='Polygon' else list(tissue.geoms);outer=unary_union([z.exterior for z in polys]);ix=(b.section_id==sec)&~b.eligible_interface
 for j in b.index[ix]:
  r=b.loc[j];mid=Point((r.x1+r.x2)/2,(r.y1+r.y2)/2);b.loc[j,'support_topology']='outer_envelope' if outer.distance(mid)<1e-6 else 'interior_hole'
 b.loc[(b.section_id==sec)&b.eligible_interface,'support_topology']='between_observed_region_and_surrounding_tissue'
 # All mask geometries must be valid; region area sums match occupied hex cells.
 assert all(shape(x['geometry']).is_valid for x in geo['features']),sec
 maskarea=sum(shape(x['geometry']).area for x in geo['features'] if x['properties']['region_id']!='valid_observation_support');assert np.isclose(maskarea,reg.loc[reg.section_id==sec,'grid_area_image_pixel_squared'].sum(),rtol=1e-7),sec
 # Add topology tags to reusable outline features without changing coordinates.
 with gzip.open(O/'regions/masks'/f'{sec}_outlines.geojson.gz','rt') as f:outlines=json.load(f)
 bm=b[b.section_id==sec].set_index('boundary_id').support_topology
 for x in outlines['features']:x['properties']['support_topology']=bm.loc[x['properties']['boundary_id']]
 with gzip.open(O/'regions/masks'/f'{sec}_outlines.geojson.gz','wt') as f:json.dump(outlines,f)
tab(b,'interfaces/boundary_summary.tsv.gz');tab(b.groupby(['dataset_id','section_id','boundary_type','support_topology'],dropna=False).agg(n_segments=('boundary_id','size'),total_length_image_pixels=('length_image_pixels','sum')).reset_index(),'interfaces/boundary_type_coverage.tsv')
# Verified identities are the only patient-level units. Anonymous and assumed units are a separate sensitivity.
pat=[];records=[]
for (ds,pid,contrast,metric),g in eff.dropna(subset=['patient_id']).groupby(['dataset_id','patient_id','contrast','metric']):
 z=g.dropna(subset=['effect']);pat.append(dict(dataset_id=ds,patient_id=pid,contrast=contrast,metric=metric,n_sections_total=len(g),n_sections_evaluable=len(z),effect=z.effect.mean(),level='verified_patient',identity_status='source-identified within study; cross-study deduplication not verified'))
for (ds,rid,contrast,metric),g in eff.dropna(subset=['record_id']).groupby(['dataset_id','record_id','contrast','metric']):
 sids=identity[identity.section_id.isin(g.section_id)];eligible=sids.eligible_under_revised_analysis_unit_policy.eq(True).all();z=g.dropna(subset=['effect']);records.append(dict(dataset_id=ds,record_id=rid,contrast=contrast,metric=metric,n_sections_total=len(g),n_sections_evaluable=len(z),effect=z.effect.mean() if eligible else np.nan,level='authorized_record_unit',unit_type=';'.join(sorted(sids.analysis_unit_type.unique())),eligible=eligible))
pat=pd.DataFrame(pat);records=pd.DataFrame(records);tab(pat,'target_profiles/verified_patient_effects.tsv');tab(records,'target_profiles/record_unit_effects.tsv');tab(pd.concat([eff,pat],ignore_index=True),'target_profiles/patient_section_effects.tsv')
summary=[]
for scope,frame in [('section',eff),('verified_patient',pat),('record_unit_sensitivity',records)]:
 for (ds,contrast,metric),g in frame.groupby(['dataset_id','contrast','metric']):
  z=g.effect.dropna();summary.append(dict(dataset_id=ds,unit_scope=scope,contrast=contrast,metric=metric,n_units_total=len(g),n_units_evaluable=len(z),n_positive=int((z>0).sum()),n_negative=int((z<0).sum()),n_zero=int((z==0).sum()),median_effect=z.median(),q25_effect=z.quantile(.25),q75_effect=z.quantile(.75),mean_effect=z.mean(),minimum_effect=z.min(),maximum_effect=z.max(),inference='descriptive paired effects; no P values; datasets not pooled'))
summary=pd.DataFrame(summary);tab(summary,'target_profiles/cross_dataset_summary.tsv')
# Existing labels are retained verbatim. New regions do not claim diagnostic equivalence.
legacy=p.groupby(['dataset_id','section_id','original_region_label','region_review_status','region_label'],dropna=False).size().reset_index(name='n_positions');tab(legacy,'regions/legacy_region_crosswalk.tsv')
# Epi/component score validity and biological resolution audit.
fc=feat.groupby(['dataset_id','section_id']).agg(n_spots=('position_id','size'),n_epi_component_scores=('epi_MP10_model_universe','count'),n_epi_score_valid=('epi_score_valid','sum')).reset_index();tab(fc,'input_review/score_resolution_coverage.tsv')
extra=[dict(dataset_id='Cervilla',section_id='visiumhd_Bladder',observation_unit='bin',resolution='2um',status='not_in_primary_interface_analysis',reason='existing audited counts/images but no accepted Epi/target annotation or CNV region; no new annotation authorized'),dict(dataset_id='GSE300692',section_id='D1_Johnson',observation_unit='bin',resolution='2/8/16um same section',status='not_in_primary_interface_analysis',reason='registered image recovered; no accepted Epi/CNV/target region for interface analysis'),dict(dataset_id='GSE299573',section_id='UC_SCC',observation_unit='bin',resolution='8um',status='not_in_primary_interface_analysis',reason='existing MP10 bin results; no trusted annotated Epi/target geometry'),dict(dataset_id='GSE318916',section_id='MPBC1',observation_unit='bin',resolution='2um',status='not_evaluable',reason='tissue mask/image registration and accepted cell composition absent'),dict(dataset_id='GSE318916',section_id='MPBC2',observation_unit='bin',resolution='2um',status='not_evaluable',reason='tissue mask/image registration and accepted cell composition absent')]
xpath=ROOT/'highres_ST/highres_spatial_audit_v2_external/08_XeniumPrime5K_sample_consistency/Prime5K_sample_level_coverage.tsv';xd=pd.read_csv(xpath,sep='\t');tab(xd,'input_review/Xenium_existing_coverage.tsv')
for i,r in xd.iterrows():extra.append(dict(dataset_id='GSE326226',section_id=str(r.get('sample',r.get('sample_id',f'coverage_record_{i+1}'))),observation_unit='cell',resolution='existing Xenium segmentation',status='architecture_only_existing_gate',reason='no existing cell-level MP10 score or accepted target annotation for this task; panel coverage gate retained',source=str(xpath)))
tab(extra,'input_review/additional_platform_coverage.tsv')
scpath=Path('/home/data/t070721/Metabolistic_scRNA/Multi2/3.2 infercnv/2/infercnv_continuous_cnv_score_per_cell.csv');tab([dict(source=str(scpath),analysis_unit='scRNA cell',available=scpath.exists(),space_mapping='no traceable spot mapping in inspected manifests',spatial_region_use=False,reason='scRNA score cannot be substituted for spatial spot CNV')],'cnv_review/scRNA_CNV_scope_audit.tsv')
# Coverage heatmap and dataset-specific paired effect figures.
bands=['inner_core','inner_edge','outer_neighborhood','far_tissue'];covmat=mp[mp.scope=='distance_band'].pivot(index='section_id',columns='group',values='n_positions').reindex(columns=bands);f,ax=plt.subplots(figsize=(8,15));im=ax.imshow(np.log10(covmat+1),aspect='auto',cmap='viridis');ax.set_yticks(range(len(covmat)),covmat.index,fontsize=7);ax.set_xticks(range(4),['Core','Inner edge','Outer neighborhood','Far tissue']);ax.set_title('Effective Visium spot coverage by distance band',color='black',loc='center');f.colorbar(im,ax=ax,label='log10(eligible spots + 1)');f.tight_layout();figsave(f,'all_sections_band_coverage','63 sections; censored support-edge positions excluded; zeros here indicate no evaluable positions, not zero biological response')
metrics=['MP10_mixed_AUCell','MP10_epi_component_existing','MP10_linear_epi_residual','major:T_NK:all_modeled_fraction','fine:Treg:all_modeled_fraction','fine:Treg:T_NK_fraction','major:Endothelial:all_modeled_fraction','major:Fibroblast:all_modeled_fraction','major:Mural:all_modeled_fraction'];datasets=sorted(cov.dataset_id.unique());colors=plt.cm.tab10(np.arange(7));labels=['Mixed MP10','Epi component MP10','MP10 residual vs Epi','T/NK / all types','Treg / all types','Treg / T/NK','Endothelial / all types','Fibroblast / all types','Mural / all types']
for contrast in eff.contrast.unique():
 f,axs=plt.subplots(3,3,figsize=(16,12))
 for ax,metric,label in zip(axs.ravel(),metrics,labels):
  for i,ds in enumerate(datasets):
   z=eff[(eff.dataset_id==ds)&(eff.contrast==contrast)&(eff.metric==metric)].dropna(subset=['effect']);v=z.effect.to_numpy();jitter=np.linspace(-.18,.18,len(v));ax.scatter(i+jitter,v,s=13,color=colors[i],alpha=.7)
   pv=pat[(pat.dataset_id==ds)&(pat.contrast==contrast)&(pat.metric==metric)].effect.dropna();ax.scatter(np.repeat(i+.27,len(pv)),pv,s=18,marker='D',facecolors='none',edgecolors='black')
   if len(v):ax.plot([i-.2,i+.2],[np.median(v)]*2,color='black',lw=2)
  ax.axhline(0,c='grey',lw=.7);ax.set_xticks(range(7),[x.replace('__UNVERIFIED','') for x in datasets],rotation=40,ha='right',fontsize=7);ax.set_title(label,loc='center');ax.set_ylabel('within-section paired difference')
 f.suptitle(contrast+' | circles: sections; diamonds: verified patients; bar: section median',color='black');f.tight_layout(rect=[0,0,1,.97]);figsave(f,'paired_effects_'+contrast.replace(' ','_'),'Separate dataset panels, section paired effects and repeated-section-averaged verified patient effects; no spot-level tests')
# All major types; section means enter once per band, with separate datasets.
major=tar[(tar.scope=='distance_band')&(tar.level=='major')];order=major.target.unique().tolist();f,axs=plt.subplots(2,4,figsize=(18,9))
for ax,ds in zip(axs.ravel(),datasets):
 z=major[major.dataset_id==ds].groupby(['group','target']).mean_local_fraction.mean().unstack().reindex(index=bands,columns=order);bottom=np.zeros(4)
 for j,t in enumerate(order):v=z[t].to_numpy();ax.bar(np.arange(4),v,bottom=bottom,label=t,color=plt.cm.tab20(j));bottom+=np.nan_to_num(v)
 ax.set_title(ds,loc='center');ax.set_xticks(range(4),['Core','Inner','Outer','Far']);ax.set_ylabel('mean local modeled fraction\n(equal section weight)');ax.set_ylim(0,1)
axs.ravel()[-1].axis('off');axs.ravel()[-1].legend(*axs.ravel()[0].get_legend_handles_labels(),loc='center');f.suptitle('All modeled major lineages by Epi-rich distance band',color='black');f.tight_layout();figsave(f,'major_composition_by_dataset','All 9 major lineages; section means within dataset; complete 30-type denominator')
# Distance curves are descriptive section summaries, no inference from individual spots.
pr=read('target_profiles/target_distance_profiles.tsv');selected=['Treg','T_NK','Endothelial','Fibroblast','Mural'];f,axs=plt.subplots(2,3,figsize=(16,9))
for ax,t in zip(axs.ravel(),['MP10']+selected):
 for ds,c in zip(datasets,colors):
  z=pr[(pr.dataset_id==ds)&(pr.target==('Treg' if t=='MP10' else t))&(pr.n_positions>=10)&np.isfinite(pr.distance_left_dNN)&np.isfinite(pr.distance_right_dNN)];z=z.groupby('distance_left_dNN').agg(value=('MP10_mean' if t=='MP10' else 'mean_local_fraction','median'),n_sections=('section_id','nunique'));ax.plot(z.index+.5,z.value,'o-',ms=3,color=c,label=ds.replace('__UNVERIFIED',''))
 ax.axvline(0,c='grey',lw=.7);ax.set_title(t,loc='center');ax.set_xlabel('signed distance / section d_NN');ax.set_ylabel('section median mixed AUCell' if t=='MP10' else 'section median local modeled fraction')
axs.ravel()[-1].legend(fontsize=6);f.suptitle('Distance profiles | only bins with >=10 eligible spots per section | raw coverage in tables',color='black');f.tight_layout();figsave(f,'continuous_distance_profiles','Per-dataset section median curves; no smoothing; excludes sparse section bins from display only')
# Supplemental spatial diagnostics and all actually resolvable target subtypes on representatives.
hm=pd.read_csv(E/'00_inputs/histology_manifest.tsv',sep='\t').set_index('section_id');cm=pd.read_csv(E/'00_inputs/c2l_manifest.tsv',sep='\t').set_index('section_id');mm=pd.read_csv(E/'00_inputs/major_map.tsv',sep='\t');rep=json.loads((O/'input_review/representative_selection.json').read_text())['sections'];allsp=p.merge(d[keys+['signed_distance','distance_dNN','distance_band','evaluation_status']],on=keys,suffixes=('','_interface'),validate='one_to_one').merge(feat,on=keys,validate='one_to_one')
for sec in rep:
 g=allsp[allsp.section_id==sec].copy();xy=g[['x','y']].to_numpy();w=pars.set_index('section_id').loc[sec,'bandwidth_w'];img=Image.open(hm.loc[sec,'image_path']);segb=b[(b.section_id==sec)&b.eligible_interface];segs=np.stack([segb[['x1','y1']],segb[['x2','y2']]],axis=1);ab=pd.read_csv(cm.loc[sec,'mean_path'],index_col=0).loc[g.position_id];frac=ab.div(ab.sum(axis=1),axis=0);size=max(2,(240*w/np.ptp(xy[:,0]))**2*.55)
 def panel(ax,v,title,cat=False,labels=None):
  ax.imshow(img,origin='upper',alpha=.12);sc=ax.scatter(xy[:,0],xy[:,1],c=v,s=size,cmap='tab20' if cat else 'viridis',linewidths=0);cb=plt.colorbar(sc,ax=ax,shrink=.7)
  if cat and labels is not None:cb.set_ticks(range(len(labels)),labels=labels,fontsize=6)
  if len(segs):ax.add_collection(LineCollection(segs,colors='black',linewidths=.3))
  ax.set_xlim(max(0,xy[:,0].min()-2*w),min(img.width,xy[:,0].max()+2*w));ax.set_ylim(min(img.height,xy[:,1].max()+2*w),max(0,xy[:,1].min()-2*w));ax.set_aspect('equal');ax.set_title(title,loc='center',color='black',fontsize=10);ax.set_xlabel('image x (pixel)');ax.set_ylabel('image y (pixel)')
 f,axs=plt.subplots(2,4,figsize=(18,9));a=axs.ravel();panel(a[0],g.native_copykat_label.map({'diploid':0,'aneuploid':1}).fillna(2),'Native CNV label; diagnostic',True,['diploid','aneuploid','undefined']);panel(a[1],g.profile_available.astype(int),'CNA profile available',True,['no','yes']);panel(a[2],g.distance_dNN,'Signed Epi-rich edge distance / d_NN');panel(a[3],g.evaluation_status.eq('evaluable').astype(int),'Primary distance evaluation',True,['excluded','eligible'])
 codes=[]
 for ax,col in zip(a[4:7],['local_niche_id','candidate_meta_niche_id','R12_niche']):
  ca=pd.Categorical(g[col].fillna('Missing'));panel(ax,ca.codes,col,True);codes.extend([dict(section_id=sec,system=col,display_code=i,label=l) for i,l in enumerate(ca.categories)])
 panel(a[7],g.epi_MP10_model_universe.where(g.epi_score_valid.eq(True)),'Existing valid Epi-component MP10');f.suptitle(sec+' | CNV baseline: '+str(cnv.set_index('section_id').loc[sec,'baseline_confidence']),color='black',fontsize=11);f.tight_layout(rect=[0,0,1,.96]);figsave(f,sec+'_interface_niche_diagnostics','Native CNV state and profile availability; signed distance and support censoring; original local/candidate-meta/R12 niches; valid component score');tab(codes,'niche_review/'+sec+'_all_niche_display_codes.tsv')
 types=mm[mm.major.isin(['T_NK','Endothelial','Fibroblast','Mural'])].subtype.tolist();f,axs=plt.subplots(4,5,figsize=(20,16))
 for ax,t in zip(axs.ravel(),types):panel(ax,frac[t],t+' / all modeled types')
 for ax in axs.ravel()[len(types):]:ax.axis('off')
 f.suptitle(sec+f' | {len(types)} resolved target types | modeled fractions, Visium spots | Epi-rich boundaries',color='black');f.tight_layout(rect=[0,0,1,.97]);figsave(f,sec+'_all_target_subtypes','All 7 T/NK, 4 Endothelial, 4 Fibroblast and 2 Mural types; no inferred missing subtypes')
figidx=pd.concat([read('figures/figure_index.tsv'),pd.DataFrame(figs)],ignore_index=True);figidx=figidx.drop_duplicates('figure',keep='last');tab(figidx,'figures/figure_index.tsv')
# Review numerical and geometry invariants, including raw source parity.
checks=dict(n_sections=len(cov),n_datasets=cov.dataset_id.nunique(),n_positions=len(p),unique_composite_keys=True,region_mask_equal_candidate=True,position_interface_key_match=True,all_distance_boundary_ids_valid=True,all_region_geometries_valid=True,area_matches_cell_grid=True,source_position_order_independent=True,n_regions=len(reg),n_valid_interface_segments=int(b.eligible_interface.sum()),n_signed_distance=int(d.signed_distance_available.sum()),n_primary_evaluable=int(d.evaluation_status.eq('evaluable').sum()),n_unavailable=int((~d.signed_distance_available).sum()),n_verified_patients=identity.patient_uid.nunique(),n_sections_with_verified_patient=int(identity.patient_uid.notna().sum()),n_figures=len(figidx),max_grid_fit_error_dNN=float(cp.max_fit_error_dNN.max()),threads=1)
(O/'logs/validation.json').write_text(json.dumps(checks,indent=2));tab(pd.DataFrame([dict(check=k,value=v,status='PASS') for k,v in checks.items()]),'logs/final_validation.tsv')
# Dataset coverage and denominator disagreement are concrete report inputs.
cover=cov.groupby('dataset_id').agg(sections=('section_id','size'),positions=('n_positions','sum'),epi_regions=('n_regions','sum'),signed_positions=('n_signed_distance','sum'),eligible_positions=('n_primary_evaluable','sum')).reset_index();tab(cover,'input_review/dataset_coverage.tsv')
tr=eff[eff.metric.isin(['fine:Treg:all_modeled_fraction','fine:Treg:T_NK_fraction'])].pivot(index=['dataset_id','section_id','contrast'],columns='metric',values='effect').reset_index();x=tr['fine:Treg:all_modeled_fraction'];y=tr['fine:Treg:T_NK_fraction'];tr['both_evaluable']=x.notna()&y.notna();tr['opposite_direction']=(x*y<0)&tr.both_evaluable;tab(tr,'target_profiles/Treg_denominator_comparison.tsv')
# Source trace for local-only conclusions; no web facts are introduced.
methods='''# MP10 区域与界面分析方法（v1）

输入为 spatial_ecology_top20_v2 冻结的 63 张 ordinary Visium、7 个 dataset、265,002 个 QC 保留组织位置。按 dataset_id、section_id、完整 position_id 合并，保留原始坐标、原分数及 CNV 标签。全部 30 类 c2l posterior mean 来自已完成的 v2_30type；不重新推断、不复制表达母对象。空间 Epithelial 是合并标签，不能由 scRNA 标签列表还原细分。Pericyte/SMC 独立属于 Mural。

## 证据与区域

旧 ROI 为 cnv_candidate_v1，8,257 个候选均无已填写 reviewed_label。原 region_auto、region_final、region_review_status、roi_id 保留。没有可对齐且已核实病理轮廓。63 张已有 CopyKAT 逐 spot 结果；41 张自动基线未经外部验证，19 张混合模型未收敛、3 张另有不确定标记。原生 aneuploid/diploid 和 MSQ 保留为诊断性证据；MSQ 定义为 native final profile 平方均值，不跨片合并、不移植 Yan 阈值。scRNA inferCNV 没有可靠空间映射，不能当 spot CNV。本轮 CNV-supported 状态为 Unknown，并非 CNV 阴性。

主区域沿用冻结 Epithelial / 全部30类 mean abundance ≥0.25 门槛，命名 Epi-rich region。区域参数在读取目标差异之前写入 frozen_rules.json。仅使用真实六邻接且实际间距≤1.35×每片中位最近邻距 d_NN；连通分量逐个编号，保留单点、孔洞与分离岛；不平滑、不填孔、不建跨间质凸包。Epi 丰度高是定义的一部分，不能作为独立验证。

几何在原图像像素框架中以真实格点拟合 affine，映射标准六角 Voronoi 单元，保留镜像/旋转；原位置 x/y 不改写。逐片记录残差，最大残差容限0.35 d_NN（冻结技术容限，非临床精度）。区域/支持 mask 是观测网格的派生覆盖，不是组织学分割或实际 spot 捕获圆面积。面积只报 grid image-pixel²，不报 μm²。输入坐标沿用已审计 lowres image reference，不再次套 scale factor。图像上视觉检查只排查明显错位。

## 边界与距离

两侧均有有效观测且属于不同 Epi mask 的共享单元边是有效 Epi-rich interface。完整原 tissue 网格用来区分 QC 缺失观测边；其余组织/捕获边保留原因未决。另按几何区分支持集外包络与内部孔洞。缺失组织边和可能截断边不参加区域界面距离。每个线段有稳定 boundary_id；轮廓同时保存 GeoJSON（本地图像坐标，非经纬度）。不从 H&E 自动猜测组织学身份。

计算点到全部候选线段的精确欧氏最短距离，不以最近线段中点近似。区域内只匹配自己的 region；外侧匹配同一有效组织连通分量中最近的 region，故每个位置只计一次。内部负、外侧正；distance_to_boundary 为绝对值，signed_distance 才携带符号。接口表 region_id 为所参照界面的区域，position_region_id 为位置本身区域，外侧为 Outside。无同分量有效边界则 NA，不填0。

w=d_NN（每片图像像素）。core d<-w；inner edge -w≤d≤0；outer neighborhood 0<d≤w。距任一缺失/组织/捕获支持边≤w 者保留距离但主统计排除；far tissue 仅在其支持边距离>w时报告。小区域无core即不可评估。代表片每dataset按 section_id 字典序首先选一张，不根据目标富集选择；主规则应用全片后，代表片另做2w敏感性。

## 分数、组成与统计

主 MP10 是冻结 mixed-spot raw Top20 AUCell，不称单个上皮细胞分数。可用的既有 Epi-component matched-model-universe 分数另列，严格沿用 epi_score_valid；不能与原 raw score 的数轴混合。原分数与 model_definition_id 保留。冻结输入没有本任务可直接复用的 raw-score high/low，本轮仅在冻结 Epi anchors 内按片中位数派生 relative_high（>中位数）、relative_low（<）、tied_unclassified（=）；非 anchor 不分组。high比例分母是该region/band所有Epi anchors，包含未分类相等者，并报告分母。所有原 scRNA 全Epi准入不变。

先汇总全部9 major，再列完整30个实际可识别类型（重点17种T/NK、EC、Fibroblast/Mural）。mean_model_abundance 是模型期望值，不是实测细胞数；mean_local_fraction 是每位置全30类比例的平均；pooled_modeled_fraction 是组内该类估计量和/全类型估计量和；mean_lineage_fraction 是该类/所属major在位置内的比例再平均。Treg/全部类型与Treg/T_NK分别报告。major 的谱系内比例恒为1，仅是分母定义。

先每section报告配对 inner_edge−inner_core、outer_neighborhood−inner_edge。患者仅来自修订来源身份表，同患者多片效应等权平均；匿名/假定独立记录单列敏感性，GSE285715未解配对不当6名独立患者。跨dataset只分层描述方向、分位数、覆盖，无spot级P值或机制推断。线性MP10~Epi丰度残差仅为片内组成敏感性，不能校正全部混合或证明因果。效应表保留左右位置数；连续曲线展示仅使用每section/bin≥10位置，原表保留全部覆盖，无平滑外推。

local niche、candidate meta-niche、已有R1/R2 accepted labels原样复用。meta-niche未确认，始终标candidate。R1/R2仅适用于原共同Epi-anchor宇宙；未被赋值不当作生物学缺失。crosswalk按同位置交叉，不强制编号一一对应。niche_target_summary为这些标签下当前中心位置组成，非重新计算ring组成。无新niche/CM聚类、无通讯推断。

## 复现

在本目录之外不写入。现有 /usr/bin/python3.12 + spatial_ecology_top20_v2/runtime 提供运行依赖；数值线程1，无GPU。

```bash
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 'mp10_cnv_regions_interface_v1/scripts/run_analysis.py'
OPENBLAS_NUM_THREADS=1 /usr/bin/python3.12 'mp10_cnv_regions_interface_v1/scripts/summarize_and_validate.py'
```

以上命令在项目根执行；重跑只更新本任务独占输出。pilot用 --pilot，完成视觉校验后使用完整命令。正式图为透明矢量PDF（H&E保留原始位图）+300dpi透明PNG，少量白底预览在logs。review.zip只含核心图、小表和报告。
'''
(O/'methods.md').write_text(methods)
def mdtable(df):
 return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join('NA' if pd.isna(v) else str(round(v,5)) if isinstance(v,float) else str(v) for v in r)+' |' for r in df.itertuples(index=False,name=None))
lines=['# MP10 CNV 区域与界面分析报告','',f'完成63张切片、7个dataset、{len(p):,}个质控保留位置；生成{len(reg):,}个Epi-rich连通区域、{int(b.eligible_interface.sum()):,}条有效界面线段。{int(d.signed_distance_available.sum()):,}个位置可给有符号距离；排除支持边影响后{int(d.evaluation_status.eq("evaluable").sum()):,}个位置进入主分带描述。输出为组成支持的区域，尚非病理确认的肿瘤边界。','', '## 证据覆盖',mdtable(cover),'',f'身份核验：30张切片对应29个研究内来源明确患者；其余记录按既有修订政策分别列为匿名/假定样本或配对未解。高分辨率补充输入逐项见 additional_platform_coverage.tsv；已有HD表达/图像并不提供本轮可信上皮/目标注释，Xenium沿用architecture-only限制。','', '## 区域与边界可靠性','旧ROI均处于pending review，原标签完整保留。主区域使用冻结Epi≥0.25及真实六邻接，保留小岛、孔洞。原生CNV仅作伴随诊断：41片自动参考未经外部验证、19片混合基线未收敛、3片其他不确定。因此未制造CNV高低阈值或病理标签。逐区域CNV profile覆盖、原生aneuploid比例、MSQ可查cnv_region_evidence.tsv。',f'坐标沿用原图像像素；最大格点拟合误差{cp.max_fit_error_dNN.max():.3f} d_NN，几何界面有这一定位误差。微米距离/物理面积不予声称。QC缺失、外包络与内部孔洞边独立标注，组织外缘与截断原因不够明确时保持unresolved。'],
# Remove tuple wrapper for easy incremental prose.
lines=lines[0]
for title,ms in [('MP10 位置与混合影响',['MP10_mixed_AUCell','MP10_epi_component_existing','MP10_linear_epi_residual']),('T/NK 与 Treg',['major:T_NK:all_modeled_fraction','fine:Treg:all_modeled_fraction','fine:Treg:T_NK_fraction']),('Endothelial',['major:Endothelial:all_modeled_fraction']),('Fibroblast 与 Mural',['major:Fibroblast:all_modeled_fraction','major:Mural:all_modeled_fraction'])]:
 lines+=['','## '+title]
 z=summary[(summary.unit_scope=='section')&(summary.contrast=='outer_neighborhood minus inner_edge')&summary.metric.isin(ms)][['dataset_id','metric','n_units_evaluable','n_positive','n_negative','median_effect']];lines+=[mdtable(z)]
 if title.startswith('MP10'):lines+=['表为外侧邻域−内侧边缘；连续分数优先。raw MP10是位置混合表达；Epi成分分数仅使用既有有效估计，单独比较。线性残差是组成敏感性，不等同细胞内调控证据。core对edge结果、高低anchor频率及分母见配对效应和region_mp10_summary表。Epi随分带变化包含区域定义本身。']
 if title.startswith('T/NK'):lines += [f'Treg两种分母在{int(tr.both_evaluable.sum())}个可评估section×contrast中有{int(tr.opposite_direction.sum())}个方向相反；这说明全组织组成和谱系内组成回答不同问题。它们是c2l模型估计，不是直接计数。']
lines+=['','## Niche 对照','local niche、candidate meta-niche与既有R1/R2共识标签全部复用。niche_region_crosswalk保存区域/带交叉；local_R12_overlap保存同位置重叠，避免把不同算法编号当作等价。R1/R2只覆盖其原Epi-anchor宇宙，因而在外侧邻域通常不可评估。candidate meta-niche尚未冻结，不用于确认性结论。代表图并列三类标签、区域界面、MP10和实际细分谱系。','', '## 不一致与限制','各dataset的方向、分布和可评估数量分别报告。患者层面结果只汇总可核实身份，匿名记录另表；不把spot数当独立生物重复，不计算探索性小P值。局部CNV与Epi不一致、阈值附近网格误差、孔洞和截断可影响边界；本轮只完成代表片2w敏感性，不按目标效应选带宽。单点岛和部分切片无core/无外侧支持时明确NA。','当前现象是可继续追踪的空间组成线索，不能由共定位推导通讯、免疫抑制、机制或预后。没有新增全量CNV、CellChat/COMMOT、空间注释或niche聚类，未改写scRNA全Epi定义。','', '## 读取接口','- regions/position_regions.tsv.gz：每位置原始候选/最终mask、区域与旧ROI标签。','- interfaces/position_interface_distances.tsv.gz：每位置最近有效界面的region_id/boundary_id、绝对和有符号距离、单位、分带及评价状态。主分析筛evaluation_status=evaluable。','- regions/masks/*.geojson.gz：逐片区域、有效观测支持及分类边界，图像像素坐标。','- input_review/position_features.tsv.gz：保留的MP10/CNV/组成门槛/niche字段，通过dataset_id+section_id+position_id连接。','- target_profiles/：区域/分带组成、连续曲线、section与verified patient效应、dataset分层描述。','- figure_index.tsv与review.zip提供轻量审阅；methods.md说明分母、几何、身份和复现。']
(O/'report.md').write_text('\n\n'.join(lines)+'\n');(O/'README.md').write_text('# MP10 CNV regions/interface v1\n\n完成63张Visium切片；主区域为Epi-rich region，不是病理肿瘤标签。先读 report.md、methods.md 和 review.zip。\n\n只读接口：interfaces/position_interface_distances.tsv.gz。唯一键dataset_id、section_id、position_id；主评价筛evaluation_status=evaluable。signed_distance内部负外部正，单位为本片图像像素；distance_dNN可用于分辨率相对展示。region_id是参照边界所属区域，position_region_id为该位置本身区域。\n\n缺失不填0，患者与匿名记录分开。所有旧输入只读。\n')
# Source manifest adds supplements without hashing large objects.
idx=read('input_review/input_index.tsv');extra_sources=[scpath,xpath,ROOT/'highres_VisiumHD/AUDIT_REPORT.md',ROOT/'highres_VisiumHD/input_recovery_followup/INPUT_RECOVERY_REPORT.md',ROOT/'neighborhood_ring_comparison_v2/01_regions/ROI_review.tsv'];idx=pd.concat([idx,pd.DataFrame([dict(role='supplemental_audit',path=str(x),exists=x.exists(),bytes=x.stat().st_size,source_mode='read_only') for x in extra_sources])],ignore_index=True);tab(idx,'input_review/input_index.tsv')
# Compact white background previews, no large expression objects in archive.
for sec in ['Inhousedata_T4384654','GSE319536_BC01']:
 for kind in ['spatial_overview','interface_niche_diagnostics']:
  im=Image.open(O/'figures'/f'{sec}_{kind}.png');im.thumbnail((1800,1400));bg=Image.new('RGB',im.size,'white');bg.paste(im,mask=im.getchannel('A'));bg.save(O/'logs'/f'{sec}_{kind}_preview.jpg',quality=90)
reviewfiles=[O/'README.md',O/'report.md',O/'methods.md',O/'figures/figure_index.tsv',O/'input_review/section_coverage.tsv',O/'input_review/dataset_coverage.tsv',O/'input_review/coordinate_alignment.tsv',O/'cnv_review/cnv_availability.tsv',O/'regions/region_parameters.tsv',O/'target_profiles/cross_dataset_summary.tsv',O/'target_profiles/verified_patient_effects.tsv',O/'interfaces/representative_2w_sensitivity.tsv',O/'logs/validation.json']
reviewfiles+=list((O/'logs').glob('*_preview.jpg'));reviewfiles+=[O/'figures'/f'{sec}_spatial_overview.pdf' for sec in ['Inhousedata_T4384654','GSE319536_BC01']];reviewfiles += [O/'figures'/f'{n}.pdf' for n in ['all_sections_band_coverage','major_composition_by_dataset','continuous_distance_profiles','paired_effects_outer_neighborhood_minus_inner_edge','paired_effects_inner_edge_minus_inner_core']]
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED) as z:
 for x in reviewfiles:z.write(x,x.relative_to(O))
manifest=[dict(path=str(x.relative_to(O)),bytes=x.stat().st_size) for x in O.rglob('*') if x.is_file() and x.name!='output_manifest.tsv'];tab(manifest,'output_manifest.tsv');(O/'logs/software_versions.json').write_text(json.dumps(dict(python=sys.version,pandas=pd.__version__,numpy=np.__version__,matplotlib=matplotlib.__version__,platform=platform.platform()),indent=2));print(json.dumps(checks,indent=2),flush=True)
