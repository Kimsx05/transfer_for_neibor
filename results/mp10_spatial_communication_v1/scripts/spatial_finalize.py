from runtime import *
import numpy as np,pandas as pd,json,re,shutil,zipfile
from scipy import sparse
Cdir=O/'comparison';Cdir.mkdir(exist_ok=True)
read=lambda p:pd.read_csv(p,sep='\t')
pilot=read(O/'input_review/pilot_sections.tsv');runs=[];contrasts=[];figures=[];params=[];coverages=[];checks=[]
def check(name,value,passed=True):checks.append(dict(check=name,value=value,pass_check=bool(passed)))
for method,folder,primary in [('COMMOT','commot','primary2dNN'),('SpatialCellChat','spatialcellchat','main')]:
 for sec in pilot.section_id:
  out=O/folder/sec/primary
  if not (out/'epi_state_type_contrasts.tsv').exists():raise RuntimeError(f'Missing completed primary {out}')
  c=read(out/'epi_state_type_contrasts.tsv');contrasts.append(c);figures.append(read(out/'figure_index.tsv'))
  if method=='COMMOT':
   status=json.loads((out/'COMPLETED.json').read_text());nLR=status['n_lr'];lrstat=read(out/'position_lr_index.tsv');activeLR=int(((lrstat.total_flow_with_diagonal-lrstat.same_spot_flow)>0).sum());elapsed=status['elapsed_seconds'];rss=status['peak_rss_GiB'];rawlog=(O/'logs'/f'commot_{sec}_primary.log').read_text(errors='replace');rsses=re.findall(r'Maximum resident set size \(kbytes\): (\d+)',rawlog);rss=max([rss]+[float(v)/1024**2 for v in rsses]);pr=json.loads((out/'parameters.json').read_text());cov=read(out/'gene_lr_coverage.tsv');cov['included']=cov.passes_5pct_all_subunits
   for k,v in pr.items():params.append(dict(method=method,section_id=sec,parameter=k,value=json.dumps(v) if isinstance(v,(list,dict)) else v))
  else:
   status=read(out/'run_status.tsv').iloc[0];nLR=status.n_LR;lrstat=read(out/'position_LR_support.tsv');activeLR=int((lrstat.cross_spot_links>0).sum());elapsed=status.elapsed_seconds;res=read(O/'spatialcellchat/resource_summary.tsv');res=res[(res.section==sec)&(res.run=='main')].iloc[0];rss=res.max_process_rss_kib/1024**2;pr=read(out/'parameters.tsv');cov=read(out/'gene_lr_coverage.tsv');cov['included']=cov.selected_for_inference
   for r in pr.itertuples():params.append(dict(method=method,section_id=sec,parameter=r.parameter,value=r.value))
  cov['method']=method;cov['section_id']=sec;coverages.append(cov)
  runs.append(dict(method=method,section_id=sec,n_spots=int(pilot.set_index('section_id').loc[sec,'n_spots']),n_LR_selected=nLR,n_LR_cross_supported=activeLR,elapsed_seconds=elapsed,peak_single_process_rss_GiB=rss,status='complete'))
  for lev in ['major','fine']:
   e=read(out/f'{lev}_edges.tsv');val=e['strength' if method=='COMMOT' else 'total_mean_strength'];check(f'{method}_{sec}_{lev}_nonnegative',val.min(),(val.dropna()>=0).all());check(f'{method}_{sec}_{lev}_finite_or_NE',len(val),not np.isinf(val).any())
  for (lev,lin,d),g in c.groupby(['network_level','lineage','direction']):
   for st in ['low','high']:
    v=g['share_within_lineage_'+st].sum(min_count=1);check(f'{method}_{sec}_{lev}_{lin}_{d}_{st}_share_sum',v,np.isnan(v) or abs(v-1)<1e-10)
  wide=O/folder/sec/('wide3dNN' if method=='COMMOT' else 'range3dNN')
  if (wide/'epi_state_type_contrasts.tsv').exists():
   w=read(wide/'epi_state_type_contrasts.tsv');keys=['network_level','lineage','node','direction'];s=c.merge(w,on=keys,suffixes=('_primary','_wide'))
   s['direction_same']=np.sign(s.delta_strength_primary)==np.sign(s.delta_strength_wide);s.to_csv(Cdir/f'{method}_{sec}_range_sensitivity.tsv',sep='\t',index=False,na_rep='NA');figures.append(read(wide/'figure_index.tsv'))
cc=pd.concat(contrasts,ignore_index=True);ff=pd.concat(figures,ignore_index=True);run=pd.DataFrame(runs);par=pd.DataFrame(params);gc=pd.concat(coverages,ignore_index=True)
shared=[]
for sec in pilot.section_id:
 cl=read(O/'commot'/sec/'primary2dNN/included_lr.tsv');sl=read(O/'spatialcellchat'/sec/'main/gene_lr_coverage.tsv');sl=sl[sl.selected_for_inference]
 canon=lambda a,b,sep: ('_'.join(sorted(str(a).split(sep))),'_'.join(sorted(str(b).split(sep))))
 ck={canon(r.ligand,r.receptor,'_') for r in cl.itertuples()};sk={canon(r.ligand_subunits,r.receptor_subunits,';') for r in sl.itertuples()}
 shared.append(dict(section_id=sec,COMMOT_unique_inferred_pairs=len(ck),SCC_unique_expanded_selected_pairs=len(sk),shared_expanded_pairs=len(ck&sk),COMMOT_only=len(ck-sk),SCC_only=len(sk-ck),scope='selected/inferred measurable LR before final SCC position support; canonical all complex subunits'))
pd.DataFrame(shared).to_csv(Cdir/'cross_method_LR_scope.tsv',sep='\t',index=False)

cc.to_csv(Cdir/'epi_state_type_contrasts.tsv',sep='\t',index=False,na_rep='NA');cc[[c for c in cc if 'share' in c or c in ['method','section_id','network_level','lineage','node','direction','lineage_labels','lineage_n_computed','computed']]].to_csv(Cdir/'lineage_shares.tsv',sep='\t',index=False,na_rep='NA')
summary=[]
for key,g in cc.groupby(['method','network_level','lineage','node','direction']):
 r=dict(zip(['method','network_level','lineage','node','direction'],key));r.update(n_sections=len(g),n_computed=int(g.computed.sum()),n_high_greater=int((g.delta_strength>0).sum()),n_low_greater=int((g.delta_strength<0).sum()),n_zero=int((g.delta_strength==0).sum()),n_NE=int((~g.computed).sum()),n_verified_patient_ids=g.loc[g.patient_id_status=='VERIFIED_WITHIN_STUDY','patient_uid'].nunique(),n_unverified_record_units=g.loc[g.patient_id_status!='VERIFIED_WITHIN_STUDY','analysis_unit_uid'].nunique(),sections=';'.join(g.section_id))
 for sh in ['all_TME','within_lineage']:
  r['share_'+sh+'_increases']=int((g['delta_share_'+sh]>0).sum());r['share_'+sh+'_decreases']=int((g['delta_share_'+sh]<0).sum())
 summary.append(r)
ss=pd.DataFrame(summary);ss.to_csv(Cdir/'cross_section_summary.tsv',sep='\t',index=False)
unavailable=read(O/'input_review/spatial_label_crosswalk.tsv');unavailable=unavailable[unavailable.identifiability_status.eq('NE_NOT_IN_SPATIAL_REFERENCE')].copy()
unavailable['computed']=False;unavailable['strength_low']=np.nan;unavailable['strength_high']=np.nan;unavailable['reason_NE']='No independent composition column in accepted spatial reference; not a calculated zero'
unavailable.to_csv(Cdir/'unavailable_types.tsv',sep='\t',index=False,na_rep='NA')
# Patient-key records are kept distinct from unverified analysis units; no pooling coordinates, no spot-level tests.
cc.to_csv(Cdir/'patient_or_record_paired_results.tsv',sep='\t',index=False,na_rep='NA')
keys=['section_id','network_level','lineage','node','direction'];agree=cc[cc.method=='COMMOT'].merge(cc[cc.method=='SpatialCellChat'],on=keys,suffixes=('_COMMOT','_SCC'))
agree['paired_computed']=agree.computed_COMMOT & agree.computed_SCC & agree.delta_strength_COMMOT.notna() & agree.delta_strength_SCC.notna();agree['direction_agrees']=(np.sign(agree.delta_strength_COMMOT)==np.sign(agree.delta_strength_SCC)).astype('boolean');agree.loc[~agree.paired_computed,'direction_agrees']=pd.NA;agree.to_csv(Cdir/'cross_method_direction_comparison.tsv',sep='\t',index=False,na_rep='NA')
# Compare only homologous named labels with scRNA; do not merge unavailable/ambiguous labels.
a=read(P/'mp10_subtype_screen_round3_v2/tables/focus_lineage_contrasts.tsv');a=a[(a.scope=='ALL')&(a.network_level=='fine')]
a=a[['node','lineage','direction','low_W','high_W','delta_W','delta_share_all_TME','delta_share_within_lineage','n_matched_specimens','n_matched_recorded_IDs']].rename(columns={col:'scRNA_'+col for col in ['lineage','low_W','high_W','delta_W','delta_share_all_TME','delta_share_within_lineage','n_matched_specimens','n_matched_recorded_IDs']})
ref=cc[cc.network_level=='fine'].merge(a,on=['node','direction'],how='left')
ref['matches_scRNA_ALL_direction']=np.where(ref.scRNA_delta_W.notna(),np.sign(ref.delta_strength)==np.sign(ref.scRNA_delta_W),None)
ref.to_csv(Cdir/'scRNA_ALL_direction_reference.tsv',sep='\t',index=False,na_rep='NA')
# Table remains raw per method and links scRNA reference rather than mixing units.
cc[cc.node.isin(['Treg','CD8 Exhausted','Tfh/CXCL13-like','NK','Capillary EC','Arterial EC','Venous EC','Lymphatic EC','Matrix_CAF','iCAF','vCAF','rCAF','Pericyte','SMC'])].to_csv(Cdir/'focus_named_subtypes.tsv',sep='\t',index=False,na_rep='NA')
ff['source_data']=ff.apply(lambda r:r.get('source_matrix') if pd.notna(r.get('source_matrix')) else str(Path(r.pdf).parent.parent/('epi_state_type_contrasts.tsv' if r.figure_type=='comparison' else ('spot_communication.tsv' if r.method=='COMMOT' else 'spot_sender_receiver.tsv'))),axis=1)
ff['scale_definition']=np.where(ff.figure_type=='circle','linear strength; common maximum across High/Low/Both in this lineage/section/method','axis-specific quantitative scales; see plot labels')
ff['direction']=ff['direction'].fillna('both directions for comparison; sender/receiver for spatial panels')
ff['relative_pdf']=ff.pdf.map(lambda p:str(Path(p).relative_to(O)));ff['relative_png']=ff.png.map(lambda p:str(Path(p).relative_to(O)));ff.to_csv(O/'figure_index.tsv',sep='\t',index=False);run.to_csv(Cdir/'method_run_summary.tsv',sep='\t',index=False);par.to_csv(O/'parameters.tsv',sep='\t',index=False);gc.to_csv(O/'gene_lr_coverage.tsv',sep='\t',index=False,na_rep='NA')
for name in ['input_index.tsv','spatial_label_crosswalk.tsv','pilot_sections.tsv','type_support.tsv']:shutil.copyfile(O/'input_review'/name,O/name)
for name in ['epi_state_type_contrasts.tsv','lineage_shares.tsv','cross_section_summary.tsv']:shutil.copyfile(Cdir/name,O/name)
for r in ff.itertuples():
 check('figure_exists_'+Path(r.pdf).name,Path(r.pdf).stat().st_size,Path(r.pdf).exists())
 if r.figure_type=='circle':check('circle_mask_'+r.method+'_'+r.section_id+'_'+Path(r.pdf).stem,r.n_positive_edges,(not r.has_self_loops) and r.hidden_positive_edges==0)
circ=ff[ff.figure_type=='circle'];check('n_primary_circles',len(circ),len(circ)==108)
for key,g in circ.groupby(['method','section_id','network_level','lineage']):check('triplet_common_scale_'+str(key),g.edge_max.nunique(),len(g)==3 and g.edge_max.nunique()==1 and set(g.view)=={'High','Low','Both'})
pd.DataFrame(checks).to_csv(Cdir/'numerical_and_figure_checks.tsv',sep='\t',index=False);assert all(r['pass_check'] for r in checks)
# Lightweight markdown formatter avoids optional tabulate dependency.
def fmt(x):
 if pd.isna(x):return 'NE'
 if isinstance(x,(float,np.floating)):return f'{x:.5g}'
 return str(x).replace('|','/')
def table(df):return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x) for x in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
section_table=pilot[['section_id','dataset','n_spots','n_epi_high','n_epi_low','patient_id_status']]
report=['# MP10 空间细胞通讯试点报告','完成日期：2026-10-27（北京时间）。本轮在同一批两张完整Visium切片上完成SpatialCellChat与COMMOT，另对首片完成一次较宽作用距离敏感性。旧scRNA与空间输入均只读。主要输出为全部纳入LR的总强度，NECTIN2–TIGIT未被拿来替代总网络。','## 试点与可辨识范围',table(section_table),'两片分别为768和3917个混合spot，均20/20 MP10基因可测。GSE171351有来源核验患者标识；Inhouse仅为未核实身份的记录单位，因此不能称为两个已核验患者，也不作患者级显著性检验。选片在通讯之前依据完整性和每数据集最小计算规模；不按Treg结果选择。','每片保留完整背景。使用accepted cell2location 30类型后验均值除以全部建模类型总量的Q，未做0.1截断或dominant硬标签。三个原参考低可辨识T/NK列合并为T_NK_unresolved，保留其质量。Treg仍单列；Pericyte/SMC遵循空间原映射保留在Mural，未改成CAF。Tip EC和Cycling EC未独立建模，主圈图明确NE。其余细群多数仍属于INTERPRET_WITH_CAUTION，不等于已确认单细胞身份。','H/L只在既有Epi比例≥0.25 anchors内按冻结MP10中位数严格划分，等值/其他spot上皮成分留作EpiOther背景。每片H/L共用一个位置网络和一份TME Q。Epi_H/L准确含义是MP10高/低关联spot中的上皮成分，未分辨同spot内两种上皮细胞。','## 模型、总强度与距离',table(run),'SpatialCellChat主指标total_mean_strength为各通过原生位置QC的LR跨spot矩阵按官方Visium Avg helper计算F/K后求和；K为该LR有链接且两类型有支持的spot-pair数。没有group P二次筛选。COMMOT在完整片中对全部纳入LR联合求竞争性最优传输，主G=QᵀTcrossQ/[sumQa·sumQb−sum(QaQb)]，分母包含跨spot零边。每种方法均另存total_flow和明确分母；major与fine独立聚合，绝不相加fine均值冒充major。','实测counts仅做每spot总量1万+log1p一次，未用c2l比例或MP10替代LR表达。human库覆盖Secreted、ECM、Contact三类蛋白信号，方法实际库版本和筛选后LR数不同；SCC按原DB交互记录计数，COMMOT精确ligand/receptor对去重计数，展开复合物后的交集见cross_method_LR_scope.tsv。主分泌/ECM距离为2×dNN=200μm，接触为1.05×dNN=105μm；采用可靠Visium阵列索引及标称100μm间距，非误用图像倍率。额外1e−6μm仅处理浮点边界。邻spot接触信号不证明膜接触。','所有模型先保留位置对角项推断，汇总副本中再去除同spot项并单列原贡献。这不等于重新求解禁止同位置通讯的模型，也不能把Visium同spot项称为准确的自分泌。','两算法的绝对数值、分母和内部筛选不同，也与scRNA W不同，不能直接相减、平均或比较大小。共享数据库/去卷积时，方向一致并非独立实验验证。','## Major结果']
for method in ['SpatialCellChat','COMMOT']:
 g=cc[(cc.method==method)&(cc.network_level=='major')];report += [f'### {method}',table(g[['section_id','node','direction','strength_low','strength_high','delta_strength','share_all_TME_low','share_all_TME_high']])]
report+=['## T/NK内的Treg']
treg=cc[(cc.network_level=='fine')&(cc.node=='Treg')];report +=[table(treg[['method','section_id','direction','strength_low','strength_high','delta_strength','rank_strength_high','rank_delta_strength','share_all_TME_low','share_all_TME_high','share_within_lineage_low','share_within_lineage_high']])]
for r in treg.itertuples():report.append(f'{r.method} / {r.section_id} / {r.direction}：High强度在T/NK内排名{int(r.rank_strength_high)}/{r.lineage_n_computed}，增量排名{int(r.rank_delta_strength)}/{r.lineage_n_computed}；谱系内份额{100*r.share_within_lineage_low:.2f}%→{100*r.share_within_lineage_high:.2f}%，全TME份额{100*r.share_all_TME_low:.2f}%→{100*r.share_all_TME_high:.2f}%。这分别描述绝对水平和组成份额，不能直接称为相对特异机制。')
for r in treg.itertuples():
 fd=r.total_flow_high-r.total_flow_low
 if np.sign(fd)!=np.sign(r.delta_strength):report.append(f'{r.method} / {r.section_id} / {r.direction}的主强度与flow方向相反：主强度{r.strength_low:.6g}→{r.strength_high:.6g}，flow {r.total_flow_low:.6g}→{r.total_flow_high:.6g}。这是归一化分母与流量口径不同，不能把主mean下降写成所有通讯总量减少。')
report.append('SpatialCellChat总mean为sum_l(F_l/K_l)，一般不等于sum_l F_l / sum_l K_l；表中summed_LR_K仅作链接覆盖汇总，不是主mean的统一分母。COMMOT的G则有其单一跨位置权重分母D，两者不能混用。SCC的K计有正Q支持的保留链接数，不直接除以目标类型Q总量；COMMOT的D则按Q加权的全部跨位置配对定义。因此稀有类型在两方法中的绝对排名可以不同，不能把某算法较高名次直接解释成真实细胞通讯更强。')
report += ['T/NK对照保留CD8 Exhausted、Tfh/CXCL13-like、NK及原低可辨识质量合成的T_NK_unresolved；后者是背景集合，不应与单一生物亚型等同。完整双方向强度、增量、两个分母份额见各片comparison_T_NK面板。Treg即使从未占某spot第一也没有被删掉。',table(cc[(cc.network_level=='fine')&(cc.lineage=='T/NK')][['method','section_id','node','direction','strength_high','delta_strength','rank_strength_high','delta_share_within_lineage']])]
for lin,title in [('Endothelial','Endothelial全部可辨亚群'),('Fibroblast','Fibroblast与独立Mural'),('Mural','Pericyte/SMC原映射')]:
 report += ['## '+title,table(cc[(cc.network_level=='fine')&(cc.lineage==lin)][['method','section_id','node','direction','strength_low','strength_high','delta_strength','rank_strength_high','delta_share_all_TME','delta_share_within_lineage']])]
 for method in ['SpatialCellChat','COMMOT']:
  for sec in pilot.section_id:
   g=cc[(cc.method==method)&(cc.section_id==sec)&(cc.lineage==lin)&(cc.network_level=='fine')]
   for d,h in g.groupby('direction'):
    winner=h.sort_values('strength_high',ascending=False).iloc[0];inc=h.sort_values('delta_strength',ascending=False).iloc[0];sh=h.sort_values('delta_share_within_lineage',ascending=False).iloc[0]
    report.append(f'{method} / {sec} / {d}：High绝对强度最高为{winner.node}；最大有符号High−Low为{inc.node}（{inc.delta_strength:.5g}）；最大有符号谱系份额变化为{sh.node}（{100*sh.delta_share_within_lineage:.3f}个百分点）。若最大值仍为负，则表示下降最少，不是增强。')
report+=['## 跨片重复、算法分歧与scRNA对照',table(ss[(ss.network_level=='fine')&ss.lineage.isin(['T/NK','Endothelial','Fibroblast','Mural'])][['method','lineage','node','direction','n_high_greater','n_low_greater','n_zero','share_within_lineage_increases','share_within_lineage_decreases']]),'只有两片，逐片效应方向比试点合并排名更有解释性。表中次数是切片/记录次数，不是已核验患者重复。相同方向也可能来自Q组成与混合表达共享因素，不能由总强度推出免疫抑制、差预后或因果机制。']
for sec,g in agree[agree.paired_computed & (agree.network_level=='fine')&agree.lineage.isin(['T/NK','Endothelial','Fibroblast','Mural'])].groupby('section_id'):report.append(f'{sec}：重点细群双向比较共{len(g)}项，COMMOT与SpatialCellChat有{int(g.direction_agrees.sum())}项High−Low方向一致；其余{int((~g.direction_agrees).sum())}项保留为算法分歧，详见cross_method_direction_comparison.tsv。')
report+=['scRNA重绘结果作为独立尺度对照：ALL Epi→Treg强度和增量均为T/NK第10/11，全TME份额下降而T/NK内部略升；五队列内部份额却仅1增4减。Tip/Cycling EC在scRNA pooled强度突出但可计算队列少；空间参考缺少它们，故不能以空间NE否定scRNA信号。Capillary在scRNA两方向4/4增加；Matrix_CAF→Epi pooled强度和增量领先但三个队列内部份额都下降；Pericyte两方向4/4增加。空间本表需要按自身方向与可测标签判断是否相符，不复制scRNA名次。','## 首片距离敏感性']
for path in sorted(Cdir.glob('*_range_sensitivity.tsv')):
 s=read(path);g=s[(s.network_level=='fine')&s.lineage.isin(['T/NK','Endothelial','Fibroblast','Mural'])];report +=[f'{path.stem}：重点双向方向{int(g.direction_same.sum())}/{len(g)}保持一致。',table(g[g.node=='Treg'][['node','direction','delta_strength_primary','delta_strength_wide','direction_same','delta_share_within_lineage_primary','delta_share_within_lineage_wide']])]
report+=['敏感性只有一次较宽范围：COMMOT仅Secreted改到3dNN；SpatialCellChat原生同一个interaction.range同时作用Secreted与ECM，故两者一起扩到3dNN，Contact不变。没有把两次扰动冒称完全相同的参数试验，也未大量扫描阈值。SpatialCellChat的K会随保留链接范围变化，因此mean下降不等同于total_flow下降；两种量分别查表。','## 覆盖、局限与交付','所有Q列达到原验证类型覆盖规则（≥10个spot的mean abundance≥0.01且总abundance>1），逐位置比例不截断。独立组成列存在不等于可确定每个spot中某个LR的细胞来源。原参考只涵盖已建模类型；Q是在其内部归一化。实测spot背景与Q软归属保留，其总强度仅为此模型与采样尺度下的通讯归属。','位置内部筛选、LR数量、spot数量都不当患者H/L检验。未运行额外空间位置置换、tradeSeq、impact、新niche、CNV或预后模型。现成可靠的完整细胞级注释ROI不足，未把HD bin冒充单细胞开展可选ROI。',f'交付{len(circ)}张High/Low/Both主圈图，全部只保留Epi↔目标边、无群体自环，三图固定位置和共同线宽；另有重点比较、空间分布与敏感性面板。共{len(checks)}项数值/图形检查通过。原推断位置矩阵保留在各方法目录；review.zip只含PDF、报告和重点小表，不含表达、qs/RDS/npz或PNG。','入口：[图索引](FIGURE_INDEX.md)、[方法](methods.md)、[完整逐片表](comparison/epi_state_type_contrasts.tsv)、[跨片方向](comparison/cross_section_summary.tsv)。']
lead=['## 主要结果（探索性）']
for d,label in [('epi_to_TME','Epi→Treg'),('TME_to_epi','Treg→Epi')]:
 z=treg[treg.direction==d]
 lead.append(f'{label}：在{len(z)}个方法×切片组合中，主强度High−Low有{int((z.delta_strength>0).sum())}项为正、{int((z.delta_strength<0).sum())}项为负。排名及谱系份额逐片变化，不能把总强度或份额称为免疫抑制强度，也不能由其推断功能和因果。')
focus=cc[(cc.network_level=='fine')&cc.lineage.isin(['Endothelial','Fibroblast','Mural'])]
for d,label in [('epi_to_TME','Epi→EC/基质'),('TME_to_epi','EC/基质→Epi')]:
 g=focus[focus.direction==d];pos=[];neg=[];mixed=[]
 for n,h in g.groupby('node'):
  if len(h)!=4 or not h.computed.all():continue
  if (h.delta_strength>0).all():pos.append(n)
  elif (h.delta_strength<0).all():neg.append(n)
  else:mixed.append(n)
 lead.append(f'{label}：两片与两算法方向均增加的可辨标签为{", ".join(pos) or "无"}；均下降的为{", ".join(neg) or "无"}；存在片间或方法间不一致的为{", ".join(mixed) or "无"}。这四项比较共享输入/数据库因素，不是四次独立实验。')
lead.append('主强度与流量总和可能方向相反，尤其应看Treg的逐片F与分母。只有一位来源核验患者加一份未核验Inhouse记录，结果是有限切片探索；没有患者级显著性结论。')
report=report[:2]+lead+report[2:]
(O/'pilot_report.md').write_text('\n\n'.join(report)+'\n')
methods='''# 空间通讯方法与复现

输入详见[input_review/README.md](input_review/README.md)，涵盖counts、Q、坐标、冻结MP10与身份来源。先运行scripts/prepare_inputs.py（实际脚本名见scripts），再依次运行commot_run.py和scc_batch.sh，最后spatial_postprocess.py及spatial_finalize.py。复现脚本run_all.sh固定顺序和环境。

## COMMOT

官方源码0.0.3，commit d117445bc07eaa19109c7609e97d9e35b26e99ca；源文件未改。项目隔离runtime补POT/plotly/narwhals，复用现有科学Python库，仅按命名空间导入官方推断/database模块，避免无关tradeSeq/pysal/karateclub导入。每片完整背景，human CellChatDB、三类蛋白LR。多亚基在过滤和推断均heteromeric=True、min、下划线分隔。每亚基必须实测且在≥5% spots表达。精确ligand/receptor对去重保留一次，所有数据库记录另存。

所有LR一次联合求解，默认四种竞争COT权重各0.25，eps=0.1、rho=10、最多10000迭代、无平滑。距离字典覆盖当前源码会访问的ligand×receptor笛卡尔积；非真实LR不会形成竞争边。官方toy例与小矩阵F/D显式求和均通过。Inhouse为运行效率采用项目并发调度封装：四个完整原生COT组件分别接收相同全LR输入，并按原0.25权重与相同运算次序合并；没有拆分LR集合或修改组件算法。多亚基、多对多、零表达和边界小矩阵的逐LR与总矩阵parity最大误差均0，见commot/parallel_scheduling_parity.tsv。主/子进程整组资源另由日志抽样记录。API未暴露收敛诊断，不把有限/非负检查当作收敛证明。参数和实测资源详见逐片COMPLETED.json。

原T_total与逐LR网络保存；Tcross只为汇总将对角项置零，原求解不变。F=QᵀTcrossQ；D=outer(sumQ,sumQ)−QᵀQ；G=F/D，D=0为NE。D包括跨位置实际零边，不是有边条件均值。Q不是类型特异实测表达，纯软归属无法确定受体究竟来自同spot哪个细胞。

## SpatialCellChat

独立包0.1.0，commit be88f300464cf970b36a6165ea7d2e41a47a7511，隔离Rlib；旧CellChat1.6.1对象未改变。完整设置见[算法说明](spatialcellchat/METHOD_NOTES.md)和parameters.tsv。实测log1pCP10k，无ALRA补值；原生MERINGUE空间特征筛选和原生位置级概率/QC。没有group permutation P，主指标不做group P二次筛选。

每LR跨位置P经drop0清显式零，F=QᵀPQ，K=I(Q>0)ᵀ I(P>0) I(Q>0)。官方computeAvgCommunProb_LR_Visium_Avg给B=F/K，无链接贡献0，不可评估类型NE。总mean_strength=sum_l B，另存sum_l F和逐LR K。这不是每个真实细胞平均。major/fine分别从同一P聚合；固定LR矩阵不会因展示哪个谱系而变化。

## 距离、状态和呈现

使用真实Visium阵列坐标，名义dNN100μm。Secreted/ECM200μm，Contact105μm；1e−6μm仅数值边界容差。首片COMMOT Secreted300而ECM保持200；SpatialCellChat受原生接口限制Secreted+ECM同时300，Contact保持105。两算法逐片共背景/共LR/共状态阈值，绝不拼接切片坐标。

每个major和所有可辨识fine谱系均High/Low/Both独立圈图。显示掩膜只作用于矩阵副本，保留双向箭头，隐藏TME↔TME、Epi↔Epi、群体自环。等大节点；三图目标固定位置、颜色和共同线性映射，所有正边显示，无top-N。不同谱系、片和算法线粗不可直接比较。NE与真实0严格区分；Tip/Cycling EC不存在的空间组成列以NE注明。

两种份额分母分别为同对象共同可计算TME集及同fine谱系集；low/high使用同一集合。排名min-tie，单一子类份额1只是结构事实。所有表同时保留低/高、差、flow、分母、覆盖和身份状态。spot/LR和内部QC不当患者；这两片仅1个来源核验患者+1个未核验记录单位，不进行患者差异检验。

## 官方来源

[SpatialCellChat代码与教程](https://github.com/jinworks/SpatialCellChat/tree/be88f300464cf970b36a6165ea7d2e41a47a7511)、[官方Visium聚合实现](https://github.com/jinworks/SpatialCellChat/blob/be88f300464cf970b36a6165ea7d2e41a47a7511/R/modeling.R)、[COMMOT空间API](https://commot.readthedocs.io/en/latest/api/commot.tl.spatial_communication.html)、[COMMOT数据库API](https://commot.readthedocs.io/en/latest/api/commot.pp.ligand_receptor_database.html)、[COMMOT原论文](https://www.nature.com/articles/s41592-022-01728-4)。
'''
(O/'methods.md').write_text(methods)
index=['# 图索引','所有主圈图为High/Low/Both三张，共尺度；PDF透明矢量，PNG透明300dpi，*_white.png为少量白底预览。','']
for r in ff.itertuples():index.append(f'- [{r.method} / {r.section_id} / {r.network_level} / {r.lineage} / {r.figure_type} / {'' if pd.isna(getattr(r,'view','')) else r.view}]({Path(r.pdf).relative_to(O).as_posix()})')
(O/'FIGURE_INDEX.md').write_text('\n'.join(index)+'\n')
(O/'README.md').write_text('# MP10 空间通讯 v1\n\n已完成两片×两算法完整试点、首片距离敏感性与全部主圈图。\n\n- [中文报告](pilot_report.md)\n- [图索引](FIGURE_INDEX.md)\n- [方法](methods.md)\n- [输入与身份](input_review/README.md)\n- [逐片结果](comparison/epi_state_type_contrasts.tsv)\n\nreview.zip为轻量审核包；原矩阵留在各算法输出目录。scRNA另见 ../mp10_subtype_screen_round3_v2。\n')
# Archive contains portable reports, all PDF figures, compact required tables and methods.
files=[O/x for x in ['README.md','pilot_report.md','methods.md','FIGURE_INDEX.md','figure_index.tsv','input_index.tsv','spatial_label_crosswalk.tsv','pilot_sections.tsv','type_support.tsv','parameters.tsv','gene_lr_coverage.tsv','epi_state_type_contrasts.tsv','lineage_shares.tsv','cross_section_summary.tsv']]
files += list((O/'spatialcellchat').glob('*.tsv'))+list((O/'commot').glob('*.tsv'))+list(Cdir.glob('*.tsv'))+[O/'input_review/README.md',O/'spatialcellchat/METHOD_NOTES.md',O/'spatialcellchat/SUMMARY_NOTES.md']+[Path(p) for p in ff.pdf]
with zipfile.ZipFile(O/'review.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in dict.fromkeys(files):z.write(f,str(f.relative_to(O)))
with zipfile.ZipFile(O/'review.zip') as z:assert z.testzip() is None;assert not any(n.endswith(('.npz','.rds','.qs','.png','.mtx.gz')) for n in z.namelist())
(O/'logs/COMPLETED.txt').write_text(f'Two sections, two methods; 108 main circles; {len(ff)} total PDFs; {len(checks)} checks PASS; review {((O/"review.zip").stat().st_size/1024**2):.2f} MiB\n')
print((O/'logs/COMPLETED.txt').read_text())
