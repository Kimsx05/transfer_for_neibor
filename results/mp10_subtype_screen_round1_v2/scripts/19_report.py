from common import *
import shutil,re
s=read('samples_v2.tsv');r=read('proportion_models_v2.tsv');axes=read('selected_review_axes.tsv');d=read('selected_axes_all_network_context.tsv');q=read('candidate_dataset_evidence.tsv');pat=read('candidate_patient_differences.tsv');e=read('candidate_sample_expression.tsv.gz');old=read('old_vs_new_samples.tsv');m=read('cells_v2.tsv.gz')
def md(df):
 def fmt(x):
  if pd.isna(x):return 'NA'
  if isinstance(x,(float,np.floating)):return f'{x:.4g}'
  return str(x).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(map(str,df.columns))+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def model(ct,lev,den,comp,var='primary_asin',scope='ALL'):
 z=r[r.subtype.eq(ct)&r.annotation_level.eq(lev)&r.denominator_type.eq(den)&r.comparison.eq(comp)&r.analysis_variant.eq(var)&r.scope.eq(scope)]
 return z.iloc[0] if len(z) else None
# Complete evidence for each selected axis with actual models, per-dataset differences, patient support and drivers.
rows=[]
for ar in axes.itertuples():
 lev='major' if ar.network_level=='major' else 'subtype';den='TME' if lev=='major' else 'parent';a=r[r.subtype.eq(ar.TME_type)&r.annotation_level.eq(lev)&r.denominator_type.eq(den)&r.comparison.str.contains('MP10')];ctx=d[d.axis_id.eq(ar.axis_id)];eq=q[q.network_level.eq(ar.network_level)&q.TME_type.eq(ar.TME_type)&q.direction.eq(ar.direction)&q.LR_id.eq(ar.LR_id)];pp=pat[pat.network_level.eq(ar.network_level)&pat.TME_type.eq(ar.TME_type)&pat.direction.eq(ar.direction)&pat.LR_id.eq(ar.LR_id)];xx=e[e.network_level.eq(ar.network_level)&e.TME_type.eq(ar.TME_type)&e.direction.eq(ar.direction)&e.LR_id.eq(ar.LR_id)&e.coverage_eligible]
 for ds in ['ALL']+sorted(s.dataset.unique()):
  cc=ctx[ctx.scope.eq(ds)];ee=eq if ds=='ALL' else eq[eq.dataset.eq(ds)];pt=pp if ds=='ALL' else pp[pp.dataset.eq(ds)];sp=xx if ds=='ALL' else xx[xx.dataset.eq(ds)];z=dict(axis_id=ar.axis_id,scope=ds,TME_type=ar.TME_type,parent='Plasma' if ar.TME_type=='Plasma' else 'B',network_level=ar.network_level,direction=ar.direction,LR_id=ar.LR_id,evidence_tier=ar.evidence_tier,n_paired_specimens=sp.sample_uid.nunique(),n_paired_patient_ids=len(pt),n_positive_patients=int(pt.aux_product_delta.gt(0).sum()),n_negative_patients=int(pt.aux_product_delta.lt(0).sum()),n_zero_patients=int(pt.aux_product_delta.eq(0).sum()),largest_absolute_expression_specimens=';'.join(sp.assign(abs_delta=sp.aux_product_delta.abs()).nlargest(3,'abs_delta').sample_uid),spatial_question='Compare candidate gene expression at high/low MP10 epithelial boundaries adjacent to B/Plasma niches; test proximity and avoidance independently',MP10_gene_overlap=ar.MP10_gene_overlap,interpretation_limit='descriptive inferred network; patient IDs unverified; no activation/contact inference; same-data expression support')
  for co in ['epi_state_within_specimen','sample_MP10_group_allEpi','Normal_Tumor','T1_T2']:
   v=cc[cc.comparison_type.eq(co)];z[co+'_raw_delta']=v.raw_delta.iloc[0] if len(v) else np.nan;z[co+'_status_low']=v.status_low.iloc[0] if len(v) else 'not_available';z[co+'_status_high']=v.status_high.iloc[0] if len(v) else 'not_available'
  for var in ['primary_asin','parent30_asin','parent30_logit','one_specimen_asin']:
   for co in ['sample_MP10_group_allEpi','MP10_continuous','MP10_stage_complete_unadjusted','MP10_stage_adjusted']:
    v=model(ar.TME_type,lev,den,co,var,ds);prefix=var+'_'+co
    for field in ['coefficient','CI_low','CI_high','P','q','n_specimens','status']:z[prefix+'_'+field]=v[field] if v is not None else np.nan
  rows.append(z)
save(rows,'main_candidate_evidence_v2.tsv')
# Full annotation-level evidence stays available irrespective of B shortlist.
other=read('all_TME_state_LR_descriptive.tsv');rr=r[r.analysis_variant.eq('primary_asin')];summ=[]
for (lev,ct,major),a in rr.groupby(['annotation_level','subtype','major']):
 levnet='major' if lev=='major' else 'fine';nn=other[other.network_level.eq(levnet)&other.TME_type.eq(ct)&other.scope.ne('ALL')];mm=a[a.scope.eq('ALL')&a.comparison.str.contains('MP10')]
 summ.append(dict(annotation_level=lev,subtype=ct,parent=major,n_pooled_MP10_models_estimated=int(mm.status.eq('estimated').sum()),min_pooled_MP10_q=mm.q.min(),n_datasets_with_supported_state_LR=nn.scope.nunique(),example_descriptive_LRs=';'.join(nn.LR_id.drop_duplicates().head(8)),selection_note='All applicable proportion models and full numeric networks retained; no matched-patient candidate validation outside B/Plasma focus'))
save(summ,'all_groups_evidence_summary.tsv')
# Explicit coverage for B labels, including nonrankable GC B.
cnt=read('expression_fine_counts.tsv');state_ids=set(s.loc[s.epi_state_eligible,'sample_uid']);bc=cnt[cnt.major.isin(['B','Plasma'])].copy();bc['target_ge30']=bc.n_cells.ge(30);bc['common_epi_specimen']=bc.sample_uid.isin(state_ids);save(bc,'B_Plasma_target_coverage.tsv')
for name in ['stacked_bar_summary','stacked_bar_boundaries','stacked_bar_unavailable']:
 shutil.copy2(T/(name+'.tsv'),T/(name+'_v2.tsv'))
# Add specimen coverage and frozen thresholds to the lightweight epi metadata only.
epi=m[m.major.eq('Epithelial')].merge(s[['sample_uid','score_coverage','score_coverage_trusted','MP10_eligible','sample_MP10_group_allEpi','MP10_threshold_allEpi','MP10_median_allEpi','n_Epi_MP10_low','n_Epi_MP10_high','epi_state_eligible']],on='sample_uid');save(epi,'epi_MP10_cells_v2.tsv.gz')
coverage=s.groupby('dataset').agg(标本数=('sample_uid','size'),记录患者ID数=('patient_uid_recorded','nunique'),肿瘤标本=('tissue_status',lambda x:x.eq('Tumor').sum()),全上皮MP10合格=('MP10_eligible','sum'),双epi状态合格=('epi_state_eligible','sum'),补算细胞=('n_score_filled_same_definition','sum')).reset_index()
bulk=r[r.scope.eq('ALL')&r.analysis_variant.eq('primary_asin')&r.annotation_level.eq('major')&r.denominator_type.eq('TME')&r.subtype.isin(['B','Plasma'])&r.comparison.str.contains('MP10')][['subtype','comparison','coefficient','CI_low','CI_high','q','n_specimens']]
bs=r[r.scope.eq('ALL')&r.analysis_variant.isin(['primary_asin','parent30_asin','parent30_logit'])&r.annotation_level.eq('subtype')&r.denominator_type.eq('parent')&r.major.eq('B')&r.comparison.eq('MP10_stage_adjusted')][['subtype','analysis_variant','coefficient','CI_low','CI_high','q','n_specimens']]
bcounts=m[m.major.eq('B')].groupby(['dataset','subtype']).size().unstack(fill_value=0).reset_index()
axsummary=axes[['axis_id','network_level','TME_type','direction','LR_id','n_datasets_same_direction','n_supported_patients','support_datasets']]
paired=q.merge(axes[['network_level','TME_type','direction','LR_id','axis_id']],on=['network_level','TME_type','direction','LR_id']);pbrief=paired[paired.axis_id.isin(['AX01','AX02','AX03','AX04'])][['axis_id','dataset','n_patients','n_supported_patients','n_positive','n_negative','raw_delta','supported_median_delta','coverage_concordant']]
mat=read('matched_B_Plasma_LR_sensitivity.tsv');niche=read('existing_niche_B_Plasma_context.tsv');niche=niche[niche.niche.isin(['niche1','niche2'])]
# Per-dataset primary vs sample-group directions must be shown even when inconsistent.
co=d[d.axis_id.isin(['AX01','AX02','AX03','AX04'])&d.scope.ne('ALL')&d.comparison_type.isin(['epi_state_within_specimen','sample_MP10_group_allEpi'])].pivot(index=['axis_id','scope'],columns='comparison_type',values='raw_delta').reset_index()
olds=old[old.reentered][['dataset','orig.ident','old_n_malignant','n_epi_total','MP10_sample_allEpi','sample_MP10_group_allEpi']]
flips=old[old.group_flipped][['dataset','orig.ident','old_MP10_sample','MP10_sample_allEpi','old_MP10_group','sample_MP10_group_allEpi']]
text=f'''# MP10 round1 v2 中文分析报告

## 本轮结论

本轮已按全上皮定义完成模块0–3；v1保持只读。**标本级 MP10 与 B/Plasma 丰度没有得到经多重校正的关联证据，但标本内 epi 状态分析给出可以进一步验证的 LR 线索。** 目前优先考虑 epi→B/Plasma 的 **APP–CD74**；B/Memory B→epi 的 **ADGRE5–CD55** 是方向相反的通讯路线候选，Plasma→epi 的 **CD46–JAG1** 保留为更有限的候选。它们不是已证实的接触、激活或抑制机制。

APP–CD74 主状态结果在四队列偏 epi-high，HRA003620 明确偏 epi-low；共同 B/Plasma 标本复核保留 CNP、Inhouse、PRJNA 的 high 方向及 HRA 的反向，但 GSE 只剩一个记录 ID且方向改变。不能把四队列初始同向直接表述为稳健四队列验证。Memory/Naive B 的公共队列覆盖比 GC B 更适合下一轮具体空间问题；Inhouse IFN-B 是保留标签的单队列探索对象，不能等同于公共队列 B 亚群。

主状态与标本分组是不同层次：前者比较同一标本内两类上皮与共享 TME；后者比较不同标本的全上皮平均 MP10。两者不一致不自动否定任一结果，也不能相互冒充独立验证。

## 1. 实际输入、样本和评分变化

有效细胞 **374,329**，全部上皮 **141,881**；58个标本、42个记录患者ID、5个dataset。15个Normal、43个Tumor。所有患者ID仍是derived/provisional，未经生物学身份核实。Inhouse有13个肿瘤标本、12个记录ID，其中12标本符合MP10和双epi状态条件；同一记录ID的两个标本在网络保留、表达复核先平均患者差，比例主模型采用固定一标本回退。

{md(coverage)}

123,699 个上皮复用历史分数，18,182 个按相同90,011-feature counts秩公式补齐。当前全部上皮均有可信有限分数，**无缺失补0**。60个跨队列重叠细胞 counts 完全一致，分数最大绝对误差小于1e-12。三份未进入标本 MP10 模型的肿瘤标本不是缺分数，而是epi少于30：GSE326225_BC5为22、GSE326225_BC7为19、T4466585为18；其细胞与其他适用描述/条件分析仍保留。

v1有27个MP10合格肿瘤标本，v2有40个；**新增13个，原合格标本6个二分标签翻转**。阈值是在各dataset新的合格肿瘤集合上一次固定，后续不随覆盖或分期重划。全部40个合格肿瘤标本也满足两类epi各≥30。

新增资格：

{md(olds)}

原合格标本的分组翻转：

{md(flips)}

原始均值、描述中位数、旧新阈值、计数、状态覆盖和排除原因均在 `tables/old_vs_new_samples.tsv` 与 `samples_v2.tsv`。图见 `figures/overall/old_vs_new_eligibility.pdf`。

## 2. B/Plasma 丰度与 B 内部组成

全局主分析使用asin√原比例；ALL调整dataset。40个标本经重复记录ID固定回退后实际拟合39标本；分期完整模型调整前后同为25标本。B parent>0 对应37标本，分期完整25；B≥30分期完整对应19标本，两种转换使用完全相同的这19标本。

以下为major/TME，连续效应按每1原始MP10单位报告，区间为95% CI：

{md(bulk)}

B/TME的sample high−low系数约−0.0102、q=0.832；Plasma/TME约+0.0102、q=0.832。连续、分期调整和all_cells分母均未支持确定性关联；全体谱系的ALL主MP10模型也无q<0.05结果。非显著不等于证明无关系；当前样本量和队列差异限制精度。

**Inhouse与公共队列没有统一丰度方向。** Inhouse B的二分与连续系数均正；CNP、GSE、HRA均负，PRJNA为正。Plasma在Inhouse二分/连续均负；公共队列二分有正有负，连续只有PRJNA为正，其余负。逐队列结果均未提供B/Plasma的校正显著关联。原始点、低B标记及连续散点已逐队列输出。

B内部组成的分期调整敏感性：

{md(bs)}

主asin中Memory B负、Naive B正；限制B≥30后仍如此，但换成相同标本的logit，Memory B、GC B和IFN-B等部分方向改变。所有区间均跨0且q不支持确定性结论。这说明低/稀疏组成和转换口径影响解释，不能择取更好看的方向。Naive B正向也只能是待验证趋势。二分、连续及分期完整调整前/后的三种分析完整见 `B_sensitivity_effects.pdf` 和模型表，不能用敏感性替代主分析。

### B标签的可比性核查

{md(bcounts)}

Inhouse 2,785个B中2,753个标为IFN-stimulated B；公共队列该标签仅GSE1个、HRA2个。Inhouse Memory B仅3个、Naive B仅5个、GC B24个。标签本身高度队列相关，无法靠dataset调整完全恢复同一亚群的可比生物学分布。

冻结marker面板保留B身份支持；在Inhouse IFN-B中可见B身份表达，但ISG15/STAT1/MX1相对其他标签/队列的模式并不构成充分、独立的跨队列一致IFN定义。公共队列那1–2个IFN-B不能当作验证集。marker点图已标细胞和标本数；基因是否在母RNA存在、逐队列是否有非零证据另表保存。此核查既不证明标签错误，也不证明其比例差是生物学差异；沿用冻结注释并保留原技术警告。

## 3. 主状态通讯与患者表达复核

主状态有12个独立推断对象（5dataset＋ALL，major/fine各一），次级条件有56个，共68个；另有6个共同B/Plasma覆盖敏感性对象。大群从原始大群抽样重推，细群独立抽样重推。两类epi逐标本等量、共享TME只出现一次。ALL仅概览；未支持的已算概率仍保存原值，LR筛掉、节点不足或测量不明均不能当0。

下面是用于阅读的10条代表轴，不是“显著top10”。同向dataset数要求每队列≥2个表达支持记录ID、网络原始差与配对表达中位方向一致，且该队列至少一状态有内部支持。该数仍只是覆盖描述，不是差异P。

{md(axsummary)}

GC B未列为候选：合格同标本覆盖只见HRA和PRJNA各1个记录ID，没有任何dataset达到2个；两者方向还相反，**不能据此排序**。其完整状态和表达表仍保留。

四条大群轴的逐队列患者方向与网络差如下。患者正/负计数基于全部可测、足覆盖的配对表达，`n_supported_patients`另计≥10%检测支持；不把表达辅助乘积当CellChat概率。

{md(pbrief)}

**APP–CD74：** 比较的是同标本epi APP表达差及共享B/Plasma CD74端。B、Plasma和B亚群共享同一个epi APP差，因此这些多个接收群不是多次独立验证。HRA的低状态增强在配对表达中也出现，必须保留。Memory B与Naive B的主要公共队列支持来自CNP、GSE、PRJNA；HRA Memory B反向。Inhouse IFN-B只在该队列有足够覆盖。

**ADGRE5–CD55：** 反方向比较B/Memory B的共享ADGRE5端与epi CD55差。B大群在CNP、HRA、Inhouse具较完整同向支持；GSE网络无状态差，PRJNA网络正而患者表达中位方向负，不应称为跨五队列一致。Naive B反方向只有CNP满足同向覆盖，等级更低。

**Plasma→epi：** CD46–JAG1在CNP和PRJNA具high方向覆盖，HRA则偏low；它作为独立于APP接收路线的有限候选保留。Plasma的ADGRE5–CD55即使某些原始概率差为正，也不能在该队列内部支持不足时算作强支持；完整表保留它，而不强行与B大群共用最佳反向轴。

IFN-B的COL4A3–SDC4仅为Inhouse单队列探索；标签与测量可比性限制依然适用。所有候选都须回到具体发送/接收基因及空间位置验证，不能由通讯方向推断激活、抑制或直接接触。

### 一次共同B/Plasma标本限制复核

两类epi、B和Plasma各≥30后，共21标本、20记录ID。该复核比单目标覆盖更严格，GSE只剩1个ID，结果只能描述。大群代表轴结果：

{md(mat[['axis_id','scope','raw_delta','n_specimens','n_patients']])}

APP–CD74在CNP、Inhouse、PRJNA仍偏high、HRA仍偏low；GSE变为low方向且仅一个ID，提示不能忽视未匹配epi对原pool的贡献。B反方向ADGRE5–CD55总体仍偏high，但是否有内部支持和样本表达一致性需一起读。Plasma反方向CD46–JAG1在共同覆盖限制后，仅PRJNA保留非零状态差，其余队列差为已计算的0，因此不宜升级为稳健多队列主线。不同敏感性对象重新归一化，不能把绝对差大小变化解释为效应放大倍数。

## 4. 与次级标本分组、组织状态和T分期的关系

{md(co)}

表内两列分别是标本内epi状态差和不同标本MP10分组差，不能混称high/low。例如APP–CD74的B大群ALL主状态差为+0.0102，而标本MP10组差为−0.0205；CNP和PRJNA也由状态high增强转为标本组high较弱，Inhouse与GSE同为正、HRA同为负。APP–CD74可在同标本状态中偏high，而在不同标本均值分组时反向或不一致；这反映比较单位、细胞组成、受体背景及独立对象尺度均不同，不能声称“同一效应被验证/推翻”。每个队列的热图和独立对象circle还给出Normal/Tumor和T1/T2，未计算的LR为NE。

只有CNP同时提供T1和T2的队列内条件；其他队列缺一条件或缺明确T，不能制造组间比较。ALL的T比较也受dataset和样本构成限制。Normal/Tumor有配对记录ID补充，但Normal来源和治疗资料缺失，解释应保留这一限制。

Treg、FOLR2/C1QC/SPP1 TAM、CAF、内皮及其他T/NK/髓系均保留原比例模型、独立网络、IO与完整数值边。`all_groups_evidence_summary.tsv`与`all_TME_state_LR_descriptive.tsv`提供适用汇总；后者只是每对象内数值差较大的LR，不是完成了同等患者表达验证的新生物学优先排名。本轮不继承v1的Treg/FOLR2 shortlist，也不因某一pooled强边覆盖B/Plasma研究重点。

## 5. 下一轮可检验的空间问题

本轮只读取现有accepted SAMPLE_R12_LEIDEN标签与轻量组成摘要，未做空间重分析。既有粗粒度标签是GC B、Naive/Memory B、Plasma，不能直接等同于本轮细分Memory、Naive或IFN-B。读取到的相关摘要为：

{md(niche)}

因此可把 **niche2的Naive/Memory B与Plasma相关环境** 作为APP–CD74路线的对接入口；niche1虽有GC B相关组成，但本轮scRNA GC B覆盖不足，不宜作为已获支持的首选GC轴。这里的B/Plasma富集来自既有组成定义，不是新的独立发现；需要下一轮确认标签版本、切片/患者单位及映射。

建议下一轮具体检验：在控制B/Plasma总量、切片与患者后，MP10高/低上皮与这些既有niche是邻接、重叠还是回避；APP高的epi是否位于CD74阳性B/Plasma邻近边界；反方向候选的ADGRE5阳性B与CD55阳性epi、CD46阳性Plasma与JAG1阳性epi是否在相应位置匹配。HRA应预设为反向待解释队列，GSE不能当共同覆盖稳健重复；IFN-B单独作为Inhouse探索层，不与公共Naive/Memory标签直接合并。

## 6. 交付与阅读顺序

先读本报告及`methods.md`，再看`review/FIGURE_INDEX.md`。每个dataset都包含主状态和原条件图、B/Plasma覆盖/比例图、真实LR圈图和配对表达图，含负向、无支持及NE，不只收录ALL或阳性图。透明PDF与300dpi PNG留在完整目录；轻量`review.zip`只装核心PDF、关键表、方法和报告，不含原始表达、母对象或大型输入。完整数值网络和重现脚本留在v2目录，v1未改写。

本轮结论止于探索性关联及空间可检验假设，没有患者层面CellChat差异P、空间接触证据或功能因果结论；不自动开始下一轮空间模块。
'''
(O/'round1_v2_report.md').write_text(text)
# Small tables specifically split primary and secondary to prevent ambiguous labels.
allc=read('B_Plasma_LR_contrasts.tsv.gz');save(allc[allc.comparison_type.eq('epi_state_within_specimen')],'epi_state_LR_contrasts_v2.tsv.gz');save(allc[allc.comparison_type.ne('epi_state_within_specimen')],'secondary_condition_LR_contrasts_v2.tsv.gz')
print('Report and integrated evidence written')
