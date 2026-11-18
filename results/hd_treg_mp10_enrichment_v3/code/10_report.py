from common import *
import shutil,re

def rd(f):return pd.read_csv(R/'tables'/f,sep='\t')
def fmt(x):
 if pd.isna(x):return 'NA'
 if isinstance(x,(float,np.floating)):return 'Inf'if np.isinf(x)else f'{x:.4g}'
 return str(x)
def md(d):
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(x)for x in row)+' |'for row in d.itertuples(index=False,name=None))
counts=rd('label_counts_and_overlap.tsv');ef=rd('enrichment_effects.tsv');ne=rd('neighborhood_enrichment.tsv');co=rd('continuous_associations.tsv');el=rd('measurement_eligibility_by_section_radius.tsv');hot=rd('spatial_robustness.tsv');th=rd('MP10_high_thresholds.tsv');cal=rd('calibration_performance_pooled.tsv');diag=rd('label_measurement_diagnostics.tsv');ff=rd('exclusion_funnel.tsv');cur=rd('cross_distance_results.tsv');bo=rd('bootstrap_diagnostics.tsv');k=rd('standard_cross_KL.tsv');roi=rd('fixed_ROI_v3_evidence.tsv');er=rd('continuous_and_ROI_robustness.tsv')
main=ef[ef.family.eq('primary_FE')];h=hot[hot.radius_um.eq(200)&hot.label_definition.eq('FOXP3_supported')&hot.high_definition.eq('H05')&hot.removed_centers.gt(0)].groupby('section').agg(FE_min=('FE','min'),FE_max=('FE','max'),max_removed=('removed_centers','max'),hotspots=('hotspot','count')).reindex(SECTIONS).reset_index();put(h,'hotspot_influence_summary.tsv')
cr=cur[cur.label_definition.eq('FOXP3_supported')&cur.high_definition.eq('H05')&cur.distance_um.eq(200)];cc=co[co.readout.eq('Epi_MP10')&co.block_um.eq(400)];br=bo[bo.high_definition.eq('H05')&bo.label_definition.eq('FOXP3_supported')&bo.radius_um.eq(200)]
# Existing raw marker/rank evidence is linked and a small table copied for portable inspection.
shutil.copy2(F/'tables/detected_zero_UCell_examples.tsv',R/'tables/historical_zero_UCell_rank_examples.tsv')
roi_text=[]
for q in roi.itertuples():
 obs=f'FOXP3/IL2RA/CTLA4/IKZF2 原始UMI={q.FOXP3_UMI}/{q.IL2RA_UMI}/{q.CTLA4_UMI}/{q.IKZF2_UMI}; FOXP3-supported={q.FOXP3_supported}, strict={q.strict}。'
 obs+=f'200μm独立Epi64={q.r200_Epi_n}, H05={q.r200_H05_n}; 50μm独立Epi64={q.r50_Epi_n}。'
 obs+='无独立上皮读数，不能以本ROI检验MP10邻近。'if q.r200_Epi_n==0 else'仅为区域表达邻近，不能确认单个Treg身份或通讯。'
 if q.section=='Cervilla'and q.roi.startswith('candidate'):obs+='两个候选是同一热点相距约32μm的邻接测量；周围两个可评分Epi64均非H05。'
 roi_text.append(dict(section=q.section,roi=q.roi,interpretation=obs,raw_evidence_atlas=str(F/'figures'/f'{q.section}_{q.roi}_512atlas.png'),new_overlay=f'figures/{q.section}_ROI_{q.roi}.png'))
put(pd.DataFrame(roi_text),'ROI_interpretation.tsv')
ix=rd('figure_index.tsv');ir=rd('figure_index_R.tsv');ix=pd.concat([ix[~ix.figure.isin(ir.figure)],ir],ignore_index=True);put(ix,'figure_index.tsv')
(R/'FIGURE_INDEX.md').write_text('# 图件索引\n\n白底Arial；PNG300dpi；PDF文字/坐标/图例矢量，密集空间点栅格化。带图例图另有`_panel_only.pdf`和`_legend_only.pdf`。每一图的统计单位均为区域。\n\n'+md(ix[['figure','description']])+'\n\n原固定24个ROI的原始多基因、H&E或表达覆盖图册继续引用v2，路径见`tables/ROI_interpretation.tsv`。本轮24个新图为新标签与H05配对几何覆盖叠图，未替换旧原始表达证据。\n')
report=f'''# hd_treg_mp10_enrichment_v3

完成日期：{datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y-%m-%d %H:%M %Z')}。这是查看过v1/v2后的探索性区域空间分析。五张切片全部纳入；T32是表达聚合区域，不是Treg细胞。没有运行通讯或上游重分析。

## 主要结论

简化标签增加了可描述的FOXP3支持T区域，但没有形成跨切片一致的MP10-high偏好。H05主FE计划10项，{main.P.notna().sum()}项可计算限制置换P，BH后均未检出明确富集(最小q={main.q.min():.4g})。50μm有严重的测量机会不足；MPBC1的阳性可交换支持不足，MPBC2无H05曝露，不能将NA解释为没有生物邻近。

200μm连续Treg程序在D1呈小幅负向调整效应(β≈−0.036)，400/800μm区块区间均偏负；其他四片区间跨零。它与FOXP3标签的FE接近1同时成立，不能把两个读数当作同一身份指标或互相独立的验证。本轮没有新增连续检验P；v1正式连续P/q另表原样保留。

## 定义与测量

T_eligible核验与历史一致：nFeature≥325；CD3D/CD3E/CD3G/TRAC/CD2中至少2个检出且原始UMI合计≥3。连续11基因UCell使用全部T_eligible，包括零值；完整背景、average ties、maxRank1300完全复用。新标签仅T_eligible且raw FOXP3>0，没有任何附加门槛。旧strict仍为原3基因/core/旧cutoff联合规则，旧cutoff=0.19955389201641047。

{md(counts[['section','T_eligible','FOXP3_supported','strict','intersection','FOXP3_only','strict_only']])}

差异不是重新评分或更多真实Treg被证实。旧strict不要求FOXP3，因此两标签不是简单包含关系。FOXP3-only大量未超过旧UCell门槛，一部分不足3个panel基因；strict-only则可以凭其他core/程序信号满足旧规则。支持marker、检测量和程序分布见`label_expression_distributions.tsv`及skill小提琴图。

{md(diag[diag.group.eq('FOXP3_supported')][['section','n','zero_program','median_program','median_nFeature','less3_signature_genes','below_old_cutoff']])}

上述FOXP3阳性零程序分别为12/39、594/1254、72/285、3/35、39/179。实际秩例子(Cervilla_T32_000033)FOXP3=1但平均秩1356.5，超出1300窗口，UCell=0。`historical_zero_UCell_rank_examples.tsv`保留原v2例子。FOXP3单个原始计数、程序窗口及宽度要求衡量不同事实。活化/Tfh可贡献共享程序是竞争解释，不能仅靠这些分布证明真实Treg身份；本轮没有增加新硬排除条件。

## 校准分母修复

复用v2的21600个冻结混合物及原训练/样本留出，未重生成或重划分。逐条件先应用完整候选规则，以**所有T_eligible阴性混合物**作为分母。最小可行非负`score > cutoff`为0.27019811718420683；此前v2有错误条件分母的候选0.319938没有被复用。本轮候选单独保存，既不替换旧strict，也不进入FOXP3定义。

{md(cal)}

表中池化率仅是摘要，阈值选择约束是每个非空训练条件。训练唯一FP在Epi90/Tfh/5000条件，为1/30；留出唯一FP在Epi40/Tfh/5000，为1/59。小条件区间很宽；无T_eligible阴性的条件率为NA。完整逐条件分子/分母、Wilson描述区间与最小性验证见calibration表。留出已在v2看过，不是新的独立验证；模拟重复也不是患者重复。修复后仍低检出，不能据此宣称区域身份问题已解决。

## MP10和有效空间总体

历史Top20、18085基因背景、平均并列秩、K=905、Epi64测量资格和分数均未变。H05/H10分别为每片合格Epi64的严格大于q95/q90(linear/type7)，没有拆分并列值。

{md(th[['section','definition','threshold','n_eligible','n_high','actual_fraction']])}

仅同片同输入组件的真实μm质心距离；50/200均累计圆。共享任何原始bin的整对T32–Epi64剔除，原分数不作分摊。U_r需要T_eligible、原supplied circle coverage≥0.5、至少1个独立可评分Epi64。所有观察/null/曲线沿用同一来源排除。

{md(el[['section','radius_um','T_eligible','U_r','FOXP3_in_U','strict_in_U']])}

50μm首拒绝原因：

{md(ff[ff.radius_um.eq(50)&ff.count_type.eq('first_exclusion')][['section','reason','n','denominator']])}

粒度与来源排除确实限制近距离分析：Cervilla有33个中心的全部候选对因共享来源而失去资格，D1为1355，UC_SCC433，MPBC1为7，MPBC2为18。另有大量位置本来就没有达到Epi64测量资格的邻居。不能将这些位置放入非high对照。

## 主FE及历史标签比较

FE分母按唯一T32中心计数，PR比较high与有测量但非high组的标签阳性率；OR不加0.5。完整a/b/c/d、0/Inf/NA及其分母在`contingency_tables.tsv`和`enrichment_effects.tsv`。

{md(main[['section','radius_um','N','label_positive','a','near_high_denominator','FE','FE_low','FE_high','P','q']])}

以上区间为400μm空间块重采样，NA包括样本支持不足、过多退化抽样。800μm敏感性中，全部五片200μm FOXP3 FE区间仍跨1。主FE的双侧P基于9999次同组件/检测量/T信号/上皮机会/coverage分层内标签置换。分层不使用MP10高低或结果选择；每条标签的可交换数和不可交换比例单独保存。该零模型不是完整空间自相关保持模型，因此区块/热点检查不可省略。

严格标签在Cervilla/UC_SCC的200μm有效总体仅各2个，MPBC1/2为0，推断不足。D1有17个，H05 FE≈{ef[(ef.section=='D1_Johnson')&(ef.radius_um==200)&(ef.label_definition=='strict')&(ef.high_definition=='H05')].FE.iloc[0]:.3f}；这也没有经辅助家族校正的明确证据。H10及strict共30个辅助FE计划行完整保留，没有依据方向替换H05主定义。

## 邻接边、距离和混合

{md(ne[ne.high_definition.eq('H05')&ne.label_definition.eq('FOXP3_supported')][['section','radius_um','O','null_mean','null_sd','Z','P','q']])}

边数对有更多上皮邻居的中心赋予更多配对机会，不能混称FE。全部40项辅助边检验共同BH；没有明确校正后信号。零方差或可交换支持不足时Z/P=NA，观察数保留。

来源排除后的累计交叉邻近曲线固定U200，0–200μm每10μm，报告每阳性中心的高上皮邻居数及同mask限制置换基线；不在曲线峰值追加新检验。200μm观察/随机均值边比为：

{md(cr[['section','N_fixed_U200','n_positive','observed_edges','null_mean_edges','observed_over_null']])}

曲线包络是逐距离点的条件随机包络，不是同时包络，不能据个别越界点宣称新的特异距离。主50μm总体与曲线固定U200不同，数字不应直接替换。

标准形式K/L仅对前三片已知mask、supplied in-tissue观察窗做同源全局不重叠子集的8μm像素border估计。D1因此从2884个Epi64减少至1001个，子集损失很大。保留各组件及不足点数的NA；这是自实现的描述性homogeneous border估计，未解决组织非均匀强度，不是未经修改的原生Kcross显著性检验。MPBC1/2缺真实组织窗，标准K/L未执行并有明确NA；没有用包围矩形替代。适配曲线五片均完成。

双向MS/NMS全部80个组合已计算，包含观察范围内无边点。A-A为无序对、A-B计一次；零A-A分母NA，合法零跨边才为0。MS/NMS有不同的总体(coverage合格T及其组件内全部合格Epi)，不冒充FE的U_r；999次同规则随机基线仅供描述，没有额外P。NMS=1不是非均匀组织中自动的无关联基准。见`mixing_scores.tsv`与双向图。

## 连续程序及稳健性

{md(cc[['section','radius_um','N','spearman','adjusted_rank_beta','low','high']])}

原始Spearman与固定秩标准化调整均保留。调整候选顺序为T-lineage、T侧UMI、Epi侧深度、独立Epi区域数、coverage和原T区域背景；顺序残差方差≤0.1的共线列按统一规则不纳入。机会数是本轮明确加入的诊断调整，因此不把系数冒充v1原系数。实际设计矩阵保存。D1的200μm负向效应在800μm区间约[−0.0697,−0.0033]，去除原固定候选ROI后β≈−0.034；只有片内、探索性支持，未跨片重复。

200μm逐个阳性热点排除后的FE范围：

{md(h)}

热点是预定64μm距离连接的同组件阳性区域组；仅删除组内中心，剩余表达/邻域不变。另有原固定候选ROI合并范围内全部T中心排除，以及每个400/800μm空间块leave-one-out。Cervilla FE可从0.72至1.36，MPBC2从0至1，说明局部正例少、影响较大；D1/UC_SCC范围更窄但总体仍接近1。没有按影响大小删除对象后选择最终模型。

## 固定ROI与局部验证建议

原24个ROI ID/中心/候选对照身份全部保留，以512μm视野显示50/200圆。新叠图同时展示有测量的Epi64及H05，统计表明确哪些配对被来源排除。原v2图册继续提供H&E/原始多基因/成员支持图，逐ROI链接和解释在`ROI_interpretation.tsv`，没有新挑选最漂亮的重叠区域。

Cervilla两个相邻候选不是独立重复：candidate_1有FOXP3=1、IL2RA=1且strict；candidate_2 FOXP3=0、IKZF2=1而strict。两者50μm均无独立可评分Epi64；200μm各有相同的2个非H05 Epi64。这支持核查局部身份与测量来源，不支持本轮H05邻近富集。

有限后续优先级(本轮未启动)：D1 candidate_1及control_2可用于区分“FOXP3+联合严格程序”与“FOXP3+但不strict”，有原高分辨率图与较多全片测量支持；Cervilla相邻候选与control_1适合保留旧cell来源差异的局部核查；UC_SCC candidate_3核心多基因可检测，但200μm没有独立合格Epi64，主要用于身份/测量核验。MPBC候选缺H&E和原组织mask，只支持表达区域描述，不宜先进行细胞位置归因。优先级依据身份/覆盖与材料可用性，不依据富集P。

## 解释边界与完成状态

本轮在可评价T区域作背景，并控制部分T信号、检测深度和上皮机会，因此检验的是标签/程序相对于一般T区域的偏好；仍不能完整区分活化T、稀疏FOXP3和真实Treg细胞。连续分数不是比例，FOXP3标签不是细胞鉴定，空间混合不是通讯。

已完成：五片定义审核、冻结标签与MP10核验、校准分母修复、40组FE/PR/OR与40组边统计、双向80组MS/NMS、连续与空间稳健性、五片适配距离曲线、已知窗的描述性K/L尝试、24固定ROI新叠图、轻量加载及重读。标准K/L在MPBC无窗及小类型/内部点不足组合不可估计；strict多数分支推断不足；50μm的部分主效应/P/区间不可估计，均保留计划行和原因。没有因缺输入停止独立模块。

没有重分割、重聚合counts、重跑BANKSY或改变历史MP10/Treg分数；不需要新增计数守恒试验，复用原已验证聚合及来源映射。原始文件mtime/尺寸、逐ID数值、来源排除、置换分母、无序配对因子2、退化抽样和输出重读均核验。没有复制全表达矩阵或上传发布。

切片患者关系仍未知/未建立独立性，MPBC tissue_original=NA；不作五患者总体P。结论为：**未检出跨片一致的正向富集；D1有小幅负向连续线索；多个近距离/严格标签分支信息不足。** 有理由进行有限局部身份核验，没有证据据此开展全片细胞重建或宣称通讯增强。

文件入口：`code/load_result.py`、`code/load_result.R`；复跑`rerun.sh`，仅从已存表重绘`replot.sh`。主要表格位于tables，分片紧凑区域表位于各section目录。原counts、T32/Epi64成员及配对依赖v1路径在manifest/loader中明确。资源和磁盘占用见`tables/resource_and_delivery.tsv`。
'''
(R/'REPORT.md').write_text(report)
# Resource accounting is bounded to this task's logs; no server-wide audit.
rr=[]
for p in (R/'logs').glob('*.log'):
 text=p.read_text(errors='replace');mem=re.search(r'Maximum resident set size \(kbytes\): (\d+)',text);elapsed=re.search(r'Elapsed \(wall clock\) time[^\n]*: ([^\n]+)',text)
 if mem:rr.append(dict(step=p.stem,peak_RSS_KiB=int(mem.group(1)),wall_clock=elapsed.group(1)if elapsed else'not recorded',threads_limit=4))
put(pd.DataFrame(rr),'runtime_resources.tsv')
# Pair table is reused, not copied; all input values must remain exactly as read.
checks=[]
for s in SECTIONS:
 t=pd.read_parquet(R/s/'T_regions.parquet').set_index('region_id');old=pd.read_parquet(V/s/'T32_regions.parquet').set_index('region_id');e=pd.read_parquet(R/s/'Epi_regions.parquet').set_index('region_id');oe=pd.read_parquet(V/s/'Epi64_regions.parquet').set_index('region_id');assert t.index.equals(old.index)
 for col in ['Treg_main','T_lineage','Tfh','Activated_T','raw_FOXP3','x_um','y_um','nCount','nFeature']:assert t[col].equals(old[col])
 assert t.strict.equals(old.evidence.eq('Treg_like_multigene'))
 assert np.array_equal(e.MP10,oe.loc[e.index,'MP10']);checks.append(dict(section=s,IDs_unique=True,old_scores_exact=True,strict_exact=True,raw_FOXP3_exact=True,MP10_exact=True,coordinates_exact=True))
audit=pd.read_csv(R/'definitions_audit.tsv',sep='\t');a=audit[audit.path.notna()]
for row in a.itertuples():p=Path(row.path);assert p.stat().st_size==row.size and abs(p.stat().st_mtime_ns-float(row.mtime_ns))<=512
put(pd.DataFrame(checks),'delivery_numeric_reread.tsv');print('REPORT_COMPLETE_AND_SOURCE_NUMERIC_VERIFIED')
