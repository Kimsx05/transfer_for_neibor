#!/usr/bin/env python3
"""Generate an honest Chinese report and portable light review archive from existing tables."""
from summary_common import *
import argparse,datetime,zipfile,shutil,hashlib
p=argparse.ArgumentParser();p.add_argument('--allow-partial',action='store_true');a=p.parse_args();comp=O/'comparison';Z=json.loads((comp/'SUMMARY_STATUS.json').read_text());partial=Z['delivery_status']!='final'
if partial and not a.allow_partial:raise RuntimeError('Report refused: SUMMARY_STATUS is partial; full inference and all postprocessing must finish first.')
I=read(comp/'section_input_identity_coverage.tsv');S=read(comp/'section_method_status.tsv');C=read(comp/'all_section_contrasts.tsv');U=read(comp/'patient_and_unverified_record_summary.tsv');D=read(comp/'dataset_stratified_summary.tsv');E=read(comp/'dataset_equal_summary.tsv');L=read(comp/'leave_one_dataset_out.tsv');T=read(comp/'exploratory_patient_tests.tsv');A=read(comp/'method_concordance_summary.tsv');FI=read(O/'figure_index.tsv');now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')
for k in ['computed','result_available']:C[k]=C[k].map(truth)
if not partial:
 live=read(O/'logs/task_status.tsv');assert len(live)==126 and live.status.isin(['reused','completed','failed','blocked']).all()
 for r in live[live.status.isin(['reused','completed'])].itertuples():assert (job_dirs(r.method,r.section_id)[0]/'POSTPROCESS_COMPLETED.json').exists()
 assert not Z.get('common_LR_is_QA_override',False),'A QA common-LR table cannot be packaged as final.'
 assert (comp/'common_LR_sensitivity.tsv').exists(),'Native-cache common-LR sensitivity missing.'
 assert (O/'figures/summary/figure_index.tsv').exists(),'Summary figures missing.'
 assert len(I)==63 and I.section_id.nunique()==63
 tag='全量任务已终态；成功任务均已完成表图' if not S.status.isin(['failed','blocked']).any() else '全量任务已终态；存在有原因及恢复入口的失败/受阻任务'
else:tag='运行中快照：不是全量完成结果'
def fmt(v,d=3):
 if pd.isna(v):return 'NE'
 if isinstance(v,(int,np.integer)):return str(v)
 if isinstance(v,(float,np.floating)):return f'{v:.{d}g}'
 return str(v).replace('|','/').replace('\n',' ')
def mdtable(df,columns=None,maxrows=None):
 if columns:df=df[columns]
 if maxrows:df=df.head(maxrows)
 if not len(df):return '暂无可列出的记录。\n'
 lines=['| '+' | '.join(map(str,df.columns))+' |','| '+' | '.join(['---']*len(df.columns))+' |']
 lines+=['| '+' | '.join(fmt(v) for v in r)+' |' for r in df.itertuples(index=False,name=None)];return '\n'.join(lines)+'\n'
def link(p,label=None):return f'[{label or Path(p).name}]({p})'
coverage=[]
for ds,g in I.groupby('dataset',sort=True):
 v=g[g.identity_stratum.eq('verified_patient')];js=S[S.dataset.eq(ds)];coverage.append({'dataset':ds,'切片':len(g),'核验患者':v.patient_uid.nunique(),'未核验记录':int(g.identity_stratum.eq('unverified_record').sum()),'H/L可评估输入':int(g.hl_status.eq('evaluable').sum()),'network_only':int(g.hl_status.eq('network_only').sum()),'完整20/20':int(g.n_mp10_genes_measured.eq(20).sum()),'tumor/normal/unknown':'/'.join(str(int(g.tissue_condition.eq(x).sum())) for x in ['tumor','normal','unknown']),'算法成功/失败受阻/进行中':f'{int(js.status.isin(["completed","reused"]).sum())}/{int(js.status.isin(["failed","blocked"]).sum())}/{int(js.status.isin(["queued","running"]).sum())}'})
# Parameter inventory records actual cache settings rather than inferred defaults.
par=[]
for r in S.itertuples():
 out,source,var=job_dirs(r.method,r.section_id)
 if r.status not in ['completed','reused']:continue
 if (source/'parameters.json').exists():
  j=json.loads((source/'parameters.json').read_text());par.extend({'method':r.method,'section_id':r.section_id,'parameter':k,'value':json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v,'source':str(source/'parameters.json')} for k,v in j.items())
 elif (source/'parameters.tsv').exists():
  t=read(source/'parameters.tsv');par.extend({'method':r.method,'section_id':r.section_id,'parameter':row.parameter,'value':row.value,'source':str(source/'parameters.tsv')} for row in t.itertuples())
save(par,O/'parameters.tsv')
methods='''# 方法与解释边界

本轮在只读 v1 的已执行定义、隔离环境和缓存上扩展全部 63 张普通 Visium。v1 两片输入、参数与聚合口径一致的四个方法缓存以 CACHE_REFERENCE.json 指向原路径；所有新表图在 v2。母清单 selected 仅为上一轮试点选择，不是本轮准入条件。不复制大型 RNA 对象或作大文件全量 hash，不重新训练 MP10、cell2location、CNV 或区域模型。

## 输入、组成与状态

采用 hard-QC 后完整组织内位置、非负实测 counts；按每个 spot 全部实测基因的总量归一化至 10,000，再 log1p 一次。不只在 LR 基因上算 library size，不用 Harmony/scVI latent、integrated scale、负残差或 ALRA 代替实测 RNA。

accepted cell2location v2_30type 的全部 30 原始参考类型 posterior mean 以全类型总和归一化为 Q。没有 q05 替换、dominant-label 硬分类、按展示类型重归一化或逐位置丰度截断。保留 T_NK_unresolved 的不可辨识质量，Mural 包括 Pericyte/SMC；未建模 Tip/Cycling EC 明确 NE。

冻结 EpiMP10_Top20_AUCell，Epi anchor 为 Epi posterior-mean Q≥0.25；各片合格 anchors 内严格高于/低于 MP10 中位数为 high/low，等于中位数、不可评分、非 anchor 为 other。Q_EpiH+Q_EpiL+Q_EpiOther=原 Q_Epi；TME 共享同一份 Q。没有 malignant-like/CNV 过滤。EpiH/L 是高/低关联 spot 内的上皮成分，不是在同一混合 spot 分辨出两种上皮细胞。

类型支持：至少 10 个 spot 的原始 posterior mean abundance≥0.01 且全片 abundance 总和>1。Epi 状态用分配到状态的原始 Epi mean 逐一检查。这个条件只决定可评估性，不截断 Q。H/L 任一状态不足时仍推断完整位置网络，状态为 network_only，H/L 为 NE，不画空的“有效对照”圈图。单一类型不足不删除该片其他类型。19/20 MP10 仍运行并单列补充评分覆盖层，不与完整20/20主层合并。

## 原生网络与空间距离

每片每算法仅拟合一个共同背景位置网络，不分 H/L 两套网络，不按样本总体 MP10 分组，不跨片拼接坐标。真实普通 Visium 阵列 nominal pitch=100μm；Secreted Signaling/ECM-Receptor=200μm，Cell-Cell Contact=105μm，边界容差1e-6μm。QC 留下的最近邻间距不改变平台距离尺度。Contact 是数据库分类与空间范围，不声称 spot 分辨率确认细胞膜接触。

COMMOT 0.0.3，commit d117445bc07eaa19109c7609e97d9e35b26e99ca；human bundled CellChatDB 三类蛋白 LR，多亚基min规则，原表达检测比例0.05，cot_nitermax=10000、eps_p=.1、rho=10、四组件权重各.25、不smooth。每个完整原生组件接收相同完整 LR 集，不按不竞争小批次拆分 LR。官方API不直接暴露停止诊断；保留有限、非负和原实现等价性检查，不把“函数返回”称为收敛证明。

SpatialCellChat 0.1.0，commit be88f300464cf970b36a6165ea7d2e41a47a7511；原生MERINGUE feature_adjP=.05、variable_both=FALSE、use_AGAN=TRUE、Kh=.5、Hill n=1、scale_distance=.0101；native filterProbability 的位置抽样观察概率分位QC（nboot100、seed666、quantile.95，minlinks10、minsource10、mintarget10）。没有 group P 二次筛选；位置过滤不是患者 H/L 差异检验。实际逐片参数与 LR 清单以 parameters.tsv 和各原生缓存为准。

全量扩展包含一项稀疏矩阵尺寸兼容修复：SpatialCellChat 原生 computeCellDistance 的 sparseMatrix 未显式指定维度，末尾孤立位置可使矩阵缩短。v2 仅在进程内函数副本中指定 dims=c(NC,NC)，保留全部位置、既定距离以及原生对角规则；未改安装库或 v1。CanDis_2026_Sample2 的 10061×10061 恢复为 10062×10062，原有子矩阵误差为0；两张 v1 试点逐值误差为0，额外概率等价检查误差为0。MERINGUE/LR 小缓存经精确匹配604特征、1841 LR后恢复；旧失败尝试保留审计。验证表见 logs/scc_compat_qa/。


## 主聚合及分母

原生推断保存含同位置项的网络。主汇总去掉位置对角线，仅改变跨位置归属汇总，不是重新拟合“禁止同位置”的网络。major/fine 分别从同一位置网络聚合，不能相加 fine 均值冒充 major。

COMMOT：T_cross=T_total 去掉对角线；F=QᵀT_crossQ；D=outer(colSums(Q),colSums(Q))−QᵀQ；G=F/D。主指标 abundance_weighted_mean_cross_spot_strength；F 是 total_flow。D 包含空间距离约束导致的零配对，不是只统计可达配对。D=0 为 NE。

SpatialCellChat：逐 LR 去对角位置概率 P_l；F_l=QᵀP_lQ；K_l=I(Q>0)ᵀI(P_l>0)I(Q>0)；B_l=F_l/K_l；主指标 total_mean_strength=Σ_l B_l（按官方 Visium Avg 的零支持处理），并保存 ΣF_l 和各LR支持分母。主指标不是 ΣF_l/ΣK_l。类型总支持符合条件、已计算而无正支持时为0；缺类型/状态/可定义分母为NE。两方法原始数值单位不同，不能直接相减或平均。

## 对照、份额和图

每片先 H/L 再汇总；保存 mean、flow、H−L、仅当两者均正的log2(H/L)、零类别、支持、覆盖和身份。全TME份额分母取同片、同方法、同level、同方向 H/L共同可评估目标；fine谱系份额只来自相应fine和共同可评估目标。分母标签全部保存；零总分母的份额为NE，无任意伪计数。

每方法每H/L有效片：major全部8群及fine全部8谱系各画High/Low/Both三张，主边用全部实际纳入LR总mean，仅显示Epi→目标和目标→Epi，方向分弧；TME↔TME、EpiH↔EpiL、自环、EpiOther展示边隐藏但原生背景不改。三图共享目标位置、顺序、颜色、等节点尺寸和同一edge max映射，不按top-N隐藏正边。不同片、谱系和方法分别标尺度；0无边、NE明示。透明矢量PDF、300dpi PNG及少量白底预览。每片另有所有major/fine的效应与双份额比较图。

## 患者、数据集与统计

核验患者包括SOURCE_VERIFIED_PATIENT及用户接受的来源优先冲突映射，保留证据标记。同一核验患者、同组织条件和治疗状态先以切片有限log2(H/L)中位数汇总；未核验记录单列，不自动视为独立患者。tumor/normal/unknown和不同治疗状态绝不合并；未知治疗仍标unknown。多片单元方向若同时有增强和减弱记mixed；零类别和NE按全部片保留。spots、LR、两算法都不增加生物学重复数。

dataset分别报告核验患者/记录、可评估量、增强/减弱/相等/NE方向比例及零类别。dataset主效应为患者或记录有限log2(H/L)中位数；跨dataset描述给每个有有限效应的dataset同等权重，取其median-logratio的算术均值，并逐一leave-one-dataset-out。身份、组织、治疗和评分覆盖层分开。H−L仅逐片/辅助展示，跨片原始G受组织几何和完整配对分母影响，不能直接解释成batch或生物学幅度。

探索性检验预先写于EXPLORATORY_TEST_SPEC.md：同dataset/组织/治疗/complete20、每对照至少8名核验患者有有限比值才做双侧Wilcoxon和2000次患者bootstrap中位数区间；按方法×dataset×组织×治疗×方向×major/fine分族BH。不足量明确NE，不放宽条件或伪计数。主要结论仍为描述，不以显著性筛选展示类型。

## 共同 LR 与诊断

从已完成native逐LR缓存，按method×dataset×组织×治疗×评分覆盖层取实际精确多亚基复合物LR交集，逐片事后重聚合；保存交集、原清单、保留数/比例和分层片数。单片层交集等于原集，不能作为跨片稳健证据。COMMOT仍处原完整LR竞争背景，不是共同集重拟合。SCC同精确复合物可能有多个原生ID，保留其原计数而不擅自去重。方法间数据库/覆盖差异不强行等同。

诊断包括全63输入的spot UMI/nFeature、Epi/目标Q、状态支持、基因/LR覆盖、组织规模、105/200μm可达配对及Q加权可达质量占完整D比例。结果只是质量、组成和几何的描述性关联，不是批次校正，也不保证片内H/L消除局部质量混杂。保持19/20、未核验身份、unknown组织、均值与flow方向不同、dataset方向相反、覆盖不足结果。

## 恢复、资源和可追溯性

scripts/resume.sh以锁保护单调度器，按切片×算法独立缓存/日志，先2任务测量、按余量升至3；按用户后续要求，在保留原worker PID和checkpoint的情况下提升为最多5个并行任务。新SCC worker为3线程，已有6线程worker完成前保持原设置；模型输入与科学参数未变。本会话预算24CPU/128GiB不变，模型加后处理共享22CPU affinity、预留2核审查，任务预算21CPU；按实际进程树CPU/PSS和服务器余量动态准入，项目总限48CPU/200GB。execution_plan、parallel_upgrade_verified和resource_admission_latest日志记录实际执行变更。其后按用户资源要求进一步保守准入：模型内存预留合计≤80GiB，观测PSS>64GiB不启动新任务，session硬预算仍128GiB。磁盘按本账户quota（soft500GiB/hard505GiB）而非文件系统df空闲量判断，每30秒审查，给运行任务按max(0.5,spots×0.0002)GiB预留输出并另留8GiB安全余量；不足时任务保持queued并写admission_waits，不删除既有结果。此策略保留所有模型参数和完整位置背景，可能等待配额/持久存储调整，不会以降采样规避资源限制。COMMOT组件算入预算，BLAS/OpenMP禁止嵌套扩张。失败记录具体原因，降低并发/复用checkpoint合理重试后继续其他片。不抽样、不切ROI、不删LR缩规模。logs/task_status.tsv是运行状态唯一写入源；summarize/report只读它。仅所有126任务终态且成功任务后处理完成才允许最终报告。完成后依次common_lr.py、summarize.py、report.py；轻量review.zip不包含巨型原生网络或RNA副本，完整缓存路径、图索引和恢复脚本保留。
'''
(O/'methods.md').write_text(methods)
lines=[f'# MP10 全部普通 Visium 空间细胞通讯：中文报告\n\n**{tag}**\n\n生成时间：{now}。范围固定63片×2方法。以下数值来自同一次摘要快照；运行中的任务不会补零或写成完成。\n', '## 一、覆盖情况\n',mdtable(pd.DataFrame(coverage)),f'输入已准备 {len(I)}/63 片；H/L可评估 {int(I.hl_status.eq("evaluable").sum())} 片，其中完整20/20主层 {int((I.hl_status.eq("evaluable")&I.n_mp10_genes_measured.eq(20)).sum())} 片。网络成功/复用任务 {Z["n_success_jobs"]}/126；后处理完成 {Z["n_postprocessed_jobs"]}/126；失败或受阻 {Z["n_failed_or_blocked_jobs"]}。已索引 {len(FI)} 组PDF/PNG。当前可评估片×算法×类型×方向对照 {int(C.computed.sum())} 条。\n']
if partial:lines.append('**本报告只用于运行审查，不作全队列结论，也不继承两张试点的候选排序。剩余任务继续由已启动调度器执行；queued/running 均保留原状态。**\n')
lines+=['\nnetwork_only 明细（完整网络仍计算）：\n',mdtable(I[I.hl_status.eq('network_only')][['section_id','n_epi_anchors','n_epi_high','n_epi_low','threshold','hl_reason']]),'其中中位数为0的切片可有较多anchors，但严格小于中位数的low为空；这是冻结评分及严格切分的结果，不能降低门槛制造对照。D1只有3个anchors，B13无anchors。\n', '评分不足切片：\n',mdtable(I[I.n_mp10_genes_measured.lt(20)][['section_id','n_mp10_genes_measured','hl_status','tissue_condition','patient_id_status']]),'组织和身份来自可追溯源表。unknown没有自动当tumor；患者别名未核验的记录没有自动当独立患者。治疗状态缺失保持unknown。详情见 [输入/身份/覆盖总表](comparison/section_input_identity_coverage.tsv)、[方法运行状态](comparison/section_method_status.tsv) 和 [证据](input_review/tissue_treatment_evidence.tsv)。\n']
# Keep every type; state-stratified equal-dataset effects, never a selected candidate ranking.
def effect_tables(level,lineage):
 es=E[(E.network_level==level)&(E.lineage.eq(lineage) if level=='fine' else True)]
 out=[]
 for key,g in es.groupby(['tissue_condition','treatment_status','score_coverage_status','identity_stratum'],sort=True,dropna=False):
  if not g.n_datasets_finite.gt(0).any():continue
  out.append('### '+' / '.join(key)+'\n')
  tab=g[['node','method','direction','n_datasets_finite','dataset_equal_mean_log2_HL','min_dataset_log2_HL','max_dataset_log2_HL','dataset_equal_higher_fraction_evaluable','dataset_equal_lower_fraction_evaluable','dataset_opposite_directions']].copy();tab.columns=['类型','方法','方向','有限dataset数','等权log2效应','dataset最小','dataset最大','等权增强比例','等权减弱比例','dataset反向'];out.append(mdtable(tab))
 if not out:out=['暂无当前层可报告的有限dataset效应；见完整覆盖与零值表。\n']
 return out
if (O/'input_review/workbook_annotations.tsv').exists():
 lines.append('用户补充样本汇总工作簿已用于更新27片注释：NPJ 7片为tumor、NMIBC/pT1；CanDis 13片按队列一致肿瘤分类标tumor，逐片MIBC/NMIBC及分期仍unknown（不按行顺序推断对应）；Inhouse 7片仅有MIBC疾病背景，取材位置及肿瘤累及对应未解，组织仍unknown。患者身份、治疗、Q和原生推断均未改变。当前全队列为53 tumor、3 normal、7 unknown。证据层级和限制见 [workbook_annotations.tsv](input_review/workbook_annotations.tsv)。\n')
lines+=['## 二、全部 major\n','八大群均用相同处理和汇总规则。表中“等权效应”是各dataset内部患者/记录median-log2(H/L)的等权平均；不是两方法原始强度平均。正负表示片内MP10高关联上皮成分相对于低的通讯归属方向。单侧零/双侧零没有被加伪计数塞入有限比值；其方向数量在逐dataset表保留。\n']+effect_tables('major','all_major')
lines+=['## 三、各 major 内全部 fine\n','所有可辨识亚群独立计算，没有以major显著或增强为前提。T_NK_unresolved保留未分辨质量；Mural独立。Tip/Cycling EC未单独建模，全部NE。单一亚群谱系的内部份额等于1是结构结果，不能当富集证据。\n']
for lineage in LINEAGES:
 lines+=['### '+lineage+'\n']+effect_tables('fine',lineage)
lines+=['## 四、跨 dataset 一致性、异质性与方法对照\n']
valid=C[C.computed];discord=int((valid.mean_flow_direction_agrees.map(truth)==False).sum()) if 'mean_flow_direction_agrees' in valid else 0
lines.append(f'当前 {len(valid)} 条可评估片内对照中，mean 与 total_flow 的方向不一致有 {discord} 条；两者度量不同，均保留，不将flow增强替代主mean增强。全部逐片数值、双份额及明确分母标签见 [all_section_contrasts.tsv](comparison/all_section_contrasts.tsv)。\n')
opp=E[E.dataset_opposite_directions.map(truth)];lines.append(f'分层的类型×方向×方法摘要中有 {len(opp)} 个出现dataset中位效应正负相反；这是异质性提示，不自动归因为batch。\n')
if len(opp):lines.append(mdtable(opp[['method','tissue_condition','score_coverage_status','identity_stratum','network_level','lineage','node','direction','n_datasets_finite','min_dataset_log2_HL','max_dataset_log2_HL']]))
lines.append('同患者多片先在同组织/治疗内汇总。dataset表保留有效片、核验患者、未核验记录、方向与零类别计数。等权汇总与leave-one-dataset-out均按身份、组织、治疗和评分覆盖分层；少数/单dataset层不能被解释为跨dataset一致性。\n')
AA=A[A.n_both_evaluable.gt(0)]
if len(AA):
 at=AA.groupby(['tissue_condition','score_coverage_status','identity_stratum'],dropna=False).agg(共同可评估片类型方向条数=('n_both_evaluable','sum'),方向一致条数=('n_direction_agrees','sum')).reset_index();at['描述性一致比例']=at['方向一致条数']/at['共同可评估片类型方向条数'];lines.append(mdtable(at))
lines.append('上表将类型/方向展开只是描述方法一致性，条数不是生物学重复数。方法具有不同原生模型、LR覆盖和分母；方向不一致不删去，也不直接比较两方法原始数值。细分见 [方法一致性](comparison/method_concordance_summary.tsv)。\n')
common_exists=Z.get('common_LR_sensitivity_available',False) and (comp/'common_LR_sensitivity_summary.tsv').exists()
if common_exists and Z.get('common_LR_is_QA_override',False):lines.append('当前共同LR表来自明确标记的两个v1缓存QA，用来验证事后汇总链条；不能作为全量敏感性完成证据。最终调度器将从全量真实缓存重新建立各层交集。\n')
elif common_exists:
 cs=read(comp/'common_LR_sensitivity_summary.tsv');cr=read(comp/'common_LR_sensitivity.tsv');multi=cr[cr.n_stratum_sections.gt(1)] if 'n_stratum_sections' in cr else cr
 lines.append(f'共同LR事后汇总已输出 {len(cr)} 条对照，其中所属分层含多片的记录 {len(multi)} 条。原始/共同集均有限时比较方向及幅度，覆盖不足不补零。保留集合和层大小见 [共同LR表](comparison/common_LR_sensitivity.tsv)；COMMOT仍在原完整LR竞争背景，不能称作共同LR重拟合。单片层交集等于原集不构成跨片稳健性。\n')
else:lines.append('共同LR事后敏感性尚待队列完成后生成，目前不声称该项已完成。\n')
nt=int(T.p_value.notna().sum());lines.append(f'按事前条件（同dataset/组织/治疗/完整20基因、≥8名有有限比值核验患者）当前可进行探索性检验 {nt} 项。所有资格、患者量、零值排除覆盖和BH结果见 [探索性表](comparison/exploratory_patient_tests.tsv)。不足条件保持描述，不放宽条件来筛候选。\n')
lines+=['## 五、质量、组成、几何与方法限制\n','全部输入已有UMI/nFeature、Epi与目标Q和空间可达配对诊断：[状态质量](input_review/spot_quality_diagnostics.tsv)、[Q加权可达分母](input_review/geometry_denominator_diagnostics.tsv)。方法缓存还记录LR覆盖与全片组织规模。描述性关联表和图用于发现覆盖/局部质量/几何伴随变化，不是batch校正或机制证明。COMMOT原始G的全配对分母随组织规模和形状变化；片内比值也不保证消除局部质量与组成差异。\n','这里是混合Visium表达上的、细胞类型比例加权跨spot通讯归属。不能据总通讯升降直接推断免疫抑制、差预后、细胞膜接触或已验证机制。unknown组织/治疗和未核验身份限制生物学可比性；19/20评分层单列。完整背景、原同spot网络和EpiOther仍在缓存中，主图去对角只是跨位置汇总选择。\n','## 六、交付、失败与恢复\n',f'全部图： [figure_index.tsv](figure_index.tsv)，透明PDF与300dpiPNG配对。汇总图在 [figures/summary](figures/summary/)，每片图在对应method/section子目录。逐LR与原网络通过各目录 CACHE_REFERENCE.json 或原生文件定位；v1只读，未复制大型网络到review包。\n']
failed=S[S.status.isin(['failed','blocked'])]
if len(failed):lines+=['失败/受阻必须逐项恢复，不能宣称全部推断成功：\n',mdtable(failed[[c for c in ['method','section_id','status','reason','cache_dir','output_dir'] if c in failed]])]
else:lines.append('摘要快照中暂无failed/blocked任务。queued/running并不属于成功。\n')
lines.append('资源准入同时检查真实账户quota。运行时采用模型预留≤80GiB、观测PSS>64GiB不新启；磁盘保留运行输出预算和8GiB安全边际。空间不足时queued是等待资源，不能写成完成；不删除已有文件。见 [disk_quota_latest.json](logs/disk_quota_latest.json)、[admission_waits.json](logs/admission_waits.json)。\n')
lines+=['恢复：运行 `bash scripts/resume.sh`；失败重试按日志修复后用 `bash scripts/resume.sh --retry-failed`。锁保护单调度器；如果遗留worker仍活跃，先核实它们，避免重复求解。单片后处理可运行 `python3 scripts/postprocess.py METHOD SECTION --force`，不重推断；全量完成后依次运行 `python3 scripts/common_lr.py`、`python3 scripts/summarize.py`、`python3 scripts/report.py`。运行审查可加 `--allow-partial`，输出明示partial。\n','轻量 [review.zip](review.zip) 包含报告、方法、实际参数、脚本、摘要/对照/支持表、全部汇总PDF及按dataset中位组织规模预定选择的代表图；完整figure_index标明全部图路径。代表图只为便于审查，选择不依赖效应方向或显著性。完整结果仍在服务器新目录。\n']
(O/'中文报告.md').write_text('\n'.join(lines))
readme=f'''# MP10 spatial communication v2 — all 63 ordinary Visium sections

**{tag}**。生成时间 {now}。请首先阅读 [中文报告](中文报告.md)。状态基于明确的任务/后处理完成标记；运行中不是全量完成。

- [逐片输入、身份、组织、评分及状态](comparison/section_input_identity_coverage.tsv)
- [126任务运行与缓存位置](comparison/section_method_status.tsv)；实时权威表 [logs/task_status.tsv](logs/task_status.tsv)
- [全部逐片H/L、flow、双份额和支持](comparison/all_section_contrasts.tsv)
- [患者与未核验记录](comparison/patient_and_unverified_record_summary.tsv)、[dataset](comparison/dataset_stratified_summary.tsv)、[dataset等权](comparison/dataset_equal_summary.tsv)、[LODO](comparison/leave_one_dataset_out.tsv)
- [方法一致性](comparison/method_concordance_summary.tsv)、[共同LR敏感性](comparison/common_LR_sensitivity.tsv)、[质量/几何诊断](comparison/batch_geometry_diagnostics.tsv)
- [全图索引](figure_index.tsv)、[方法](methods.md)、[实际参数](parameters.tsv)、[脚本](scripts/)、[轻量review包](review.zip)

两种方法单位不同；先片内H/L，再同组织/治疗患者内汇总，再dataset。unknown不当tumor、未核验记录不当核验患者，19/20单列。0、NE与failed/blocked严格区分。所有谱系相同规则，候选判断不继承试点排序。

新目录只保存v2结果；v1缓存只读以source_dir引用。恢复入口 `bash scripts/resume.sh`；全量自动最终链 common_lr.py → summarize.py → report.py。输出中有queued/running时不得称为完成。
''';(O/'README.md').write_text(readme)
# Portable light review: numerical data + scripts + all summary PDFs and predefined representative circle subset.
review=O/'review';review.mkdir(exist_ok=True);include=set()
for fn in ['README.md','中文报告.md','methods.md','parameters.tsv','EXPLORATORY_TEST_SPEC.md','figure_index.tsv']:include.add(O/fn)
for f in comp.glob('*.tsv'):include.add(f)
for f in comp.glob('*.json'):include.add(f)
for f in (O/'scripts').glob('*'):
 if f.is_file() and f.suffix in ['.py','.R','.sh']:include.add(f)
for name in ['candidate_sections.tsv','section_status.tsv','section_identity_crosswalk.tsv','tissue_treatment_evidence.tsv','input_index_revision_audit.tsv','type_support.tsv','spatial_label_crosswalk.tsv','spot_quality_diagnostics.tsv','geometry_denominator_diagnostics.tsv','INPUT_PREPARATION_COMPLETED.json','workbook_annotations.tsv','WORKBOOK_ANNOTATIONS_APPLIED.json','README.md']:
 f=O/'input_review'/name
 if f.exists():include.add(f)
for name in ['task_status.tsv','task_status.json','INFERENCE_QUEUE_FINISHED.json','FINALIZATION_FAILED.json','execution_plan.json','parallel_upgrade_verified.json','resource_admission_latest.json','admission_waits.json','disk_quota_latest.json','quota_guard_handoff.json','input_resources_summary.json']:
 f=O/'logs'/name
 if f.exists():include.add(f)
for f in (O/'logs/scc_compat_qa').glob('*.tsv'):include.add(f)
for name in ['SCC_COMPATIBILITY_FIX.json','scc_repair_handoff.json']:
 f=O/'logs'/name
 if f.exists():include.add(f)
# Aggregate per-job QA, actual cache pointers and native LR coverage, not large LR/network objects.
qa=[];selected=[]
for r in S.itertuples():
 out,source,var=job_dirs(r.method,r.section_id)
 for fn in ['POSTPROCESS_COMPLETED.json','postprocess_checks.tsv','CACHE_REFERENCE.json','HL_STATUS.md','unavailable_types.tsv']:
  if (out/fn).exists():include.add(out/fn)
 if (out/'postprocess_checks.tsv').exists():
  q=read(out/'postprocess_checks.tsv');q['method']=r.method;q['section_id']=r.section_id;qa.append(q)
 if r.status in ['completed','reused']:
  selected.append({'method':r.method,'section_id':r.section_id,'cache_dir':str(source),'new_output_dir':str(out),'native_LR_list':str(source/('LR_selected.tsv' if r.method=='SpatialCellChat' else 'included_lr.tsv')),'native_network':str(source/('position_network_native.rds' if r.method=='SpatialCellChat' else 'position_total_with_diagonal.npz'))})
if qa:save(pd.concat(qa,ignore_index=True),review/'postprocess_QA_all_jobs.tsv');include.add(review/'postprocess_QA_all_jobs.tsv')
save(selected,review/'native_cache_inventory.tsv');include.add(review/'native_cache_inventory.tsv')
selection=[]
for ds,g in I[(I.hl_status=='evaluable')&(I.n_mp10_genes_measured==20)].groupby('dataset',sort=True):
 available=set(C[C.result_available].section_id);g=g[g.section_id.isin(available)]
 if len(g):
  median=g.n_spots.median();chosen=g.assign(delta=(g.n_spots-median).abs()).sort_values(['delta','section_id']).iloc[0];selection.append({'dataset':ds,'section_id':chosen.section_id,'n_spots':chosen.n_spots,'selection':'closest to dataset median available complete20 H/L-evaluable section size; tie section_id; no effect selection'})
SEL=pd.DataFrame(selection);save(SEL,review/'representative_selection.tsv');include.add(review/'representative_selection.tsv');selset=set(SEL.section_id) if len(SEL) else set()
FI['included_in_light_review']=False
for idx,r in FI.iterrows():
 keep=str(r.get('relative_pdf','')).startswith('figures/summary/') or r.figure_type=='network_only_status' or (r.get('section_id') in selset and r.figure_type=='circle' and (r.get('network_level')=='major' or r.get('view')=='Both'))
 if keep and Path(r.pdf).exists():include.add(Path(r.pdf));FI.loc[idx,'included_in_light_review']=True
save(FI,review/'figure_index_review.tsv');include.add(review/'figure_index_review.tsv')
(review/'README.md').write_text('''# 轻量审查包

根目录中文报告明示当前是partial还是final；不能从zip存在推断全量完成。

包含数值表、方法、参数、脚本、全部汇总PDF，以及每dataset按可用完整评分且H/L有效切片的中位组织规模选代表：major三联和所有fine的Both图。完整High/Low/Both和300dpiPNG仍在服务器，见根figure_index及review/figure_index_review.tsv；included_in_light_review标出包内PDF。不含巨大原生位置网络/RNA，review/native_cache_inventory列其路径。

代表图不依效应方向/显著性选择。脚本引用服务器绝对缓存路径，本包便于阅读和审查而不是脱离原输入独立重算所有网络。
''');include.add(review/'README.md')
# Archive manifest hashes only small delivered files, never large RNA/native network caches.
manifest=[]
for path in sorted(include):
 if path.is_file():manifest.append({'relative_path':str(path.relative_to(O)),'bytes':path.stat().st_size,'sha256_small_delivered_file':hashlib.sha256(path.read_bytes()).hexdigest()})
save(manifest,review/'manifest.tsv');include.add(review/'manifest.tsv');tmp=O/'review.zip.tmp'
with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for path in sorted(include):
  if path.is_file():z.write(path,str(path.relative_to(O)))
with zipfile.ZipFile(tmp) as z:assert z.testzip() is None;nfiles=len(z.namelist())
tmp.replace(O/'review.zip')
result={'status':'partial' if partial else 'completed_with_documented_failures' if len(failed) else 'completed','created_at':now,'n_expected_jobs':126,'n_terminal_jobs':Z['n_terminal_jobs'],'n_postprocessed_jobs':Z['n_postprocessed_jobs'],'n_figures':len(FI),'review_files':nfiles,'review_bytes':(O/'review.zip').stat().st_size,'review_integrity':'zip CRC verified; small deliverable SHA256 manifest only','summary_snapshot':str(comp/'SUMMARY_STATUS.json')}
atom_json(result,O/('DELIVERY_PARTIAL.json' if partial else 'DELIVERY_COMPLETED.json'))
if not partial and (O/'DELIVERY_PARTIAL.json').exists():(O/'DELIVERY_PARTIAL.json').unlink()
print(json.dumps(result,ensure_ascii=False),flush=True)
